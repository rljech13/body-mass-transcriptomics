#!/usr/bin/env python3
"""Working sets and fits after BOTH recovery rounds.

Gene position priority: the species' own annotation, then the round-1 miniprot
locus (chicken protein, script 50), then the round-2 locus (four-bird proteins,
script 56). Measurement rule is the same for all seven species: peaks within
WINDOW bp of the 5' end on the gene's strand, 1 TPM per library, summed, averaged
over the tissue's libraries.

Working set = genes of the KEGG map detected in every one of the seven birds.
"""
from __future__ import annotations

import argparse
import importlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
s53 = importlib.import_module("53_three_checks")
from oxphos_thread_lib import (  # noqa: E402
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    load_geneset,
    load_geo_rnaseq,
    normalize_symbol,
    peak_to_symbol_de,
)

EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
R1 = RES / "gene_recovery"
R2 = RES / "recovery_round2"
PROTSRC = Path("/path/to/storage/kleiber_birds_gemcat/ref/chicken_prot/query_proteins_all.faa")
OUT = RES / "working_sets_round2"
SETS = {"hsa00190 OXPHOS": "kegg_00190", "hsa04714 Thermogenesis": "kegg_04714"}
THR = 1.0


def loci_from_gff(gff: Path, pid2sym: dict[str, str] | None) -> dict[str, dict]:
    best: dict[str, tuple[float, dict]] = {}
    if not gff.exists():
        return {}
    for line in gff.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        at = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
        tgt = at.get("Target", "").split()[0]
        sym = tgt.split("##")[1] if "##" in tgt else (pid2sym or {}).get(tgt)
        if not sym:
            continue
        seq = f[0].split()[0]
        if seq.startswith("ENA|"):
            seq = seq.split("|")[-1]
        sc = float(f[5]) if f[5] != "." else 0.0
        st = f[6]
        cand = dict(seqname=seq, strand=st, five=float(f[3] if st == "+" else f[4]))
        if sym not in best or sc > best[sym][0]:
            best[sym] = (sc, cand)
    return {k: v[1] for k, v in best.items()}


def measure(species: str, window: int) -> dict[str, pd.Series]:
    d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
    fp: dict[str, dict] = {}
    ann = d.dropna(subset=["geneStart", "geneEnd"])
    ann = ann[ann["symbol"].ne("")].copy()
    ann["sg"] = ann["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(ann["strand"])
    for sym, g in ann.groupby("symbol", observed=True):
        r = g.iloc[0]
        st = r["sg"]
        fp[sym] = dict(seqname=r["seqnames"], strand=st,
                       five=float(r["geneStart"] if st == "+" else r["geneEnd"]))
    pid2sym = dict(re.findall(r"^>(\S+)\s+(\S+)\s+\[organism=", PROTSRC.read_text(), re.M))
    # priority: shipped annotation, then the whole-proteome re-annotation
    # (script 58), then the two targeted recovery rounds
    REANNO = Path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies")
    for gff, p2s in (
        (REANNO / species / "reanno.gff", None),
        (R1 / f"{species}_miniprot.gff", pid2sym),
        (R2 / f"{species}_q2.gff", None),
    ):
        for sym, c in loci_from_gff(gff, p2s).items():
            fp.setdefault(sym, c)

    cols = {t: lib_cols(d, t) for t in ("fib", "epi")}
    pool = d.assign(mid=(d["start"] + d["end"]) / 2.0)
    by = {k: v for k, v in pool.groupby("seqnames", observed=True)}
    out = {t: {} for t in ("fib", "epi")}
    for sym, c in fp.items():
        sub = by.get(c["seqname"])
        if sub is None:
            continue
        near = sub[(sub["strand"] == c["strand"]) & ((sub["mid"] - c["five"]).abs() <= window)]
        if near.empty:
            continue
        for t in ("fib", "epi"):
            cc = cols[t]
            if not cc:
                continue
            m = near[cc].apply(pd.to_numeric, errors="coerce").fillna(0.0)
            m = m.where(m >= THR, 0.0)
            v = float(m.sum(axis=0).mean())
            if v > 0:
                out[t][sym] = v
    return {t: pd.Series(out[t], dtype=float) for t in ("fib", "epi")}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=1000)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    mats = {sp: measure(sp, a.window) for sp in SPECIES_ORDER}
    cage = {t: pd.DataFrame({sp: mats[sp][t] for sp in SPECIES_ORDER}) for t in ("fib", "epi")}
    idx = sorted(set(cage["fib"].index) | set(cage["epi"].index))
    ok = None
    for t in ("fib", "epi"):
        m = cage[t].reindex(idx)
        g = m.notna() & (m > 0)
        ok = g if ok is None else (ok & g)
    rna = pd.DataFrame({sp: (lambda g: g["gene_mat"].div(g["lib_total"], axis=1).mul(1e6).mean(axis=1))(
        load_geo_rnaseq(sp)) for sp in SPECIES_ORDER})

    rows = []
    for name, key in SETS.items():
        genes = sorted({normalize_symbol(g) for g in load_geneset(key)} - {""})
        det = ok.reindex(genes).fillna(False).infer_objects(copy=False)
        per = det.sum(axis=0)
        core = sorted(det.index[det.all(axis=1)])
        (OUT / f"working_{key}_w{a.window}.txt").write_text("\n".join(core) + "\n")
        print(f"\n{name}: карта {len(genes)}, найдено по видам "
              + ", ".join(f"{sp} {int(per[sp])}" for sp in SPECIES_ORDER)
              + f"; РАБОЧИЙ НАБОР {len(core)}")
        for layer, w in (("CAGE фибробласты", cage["fib"]), ("CAGE эпибласт", cage["epi"]),
                         ("RNA-seq фибробласты", rna)):
            present = [g for g in core if g in w.index]
            r, p = s53.fit(np.log10(w.reindex(present).mean(axis=0).to_numpy(float)))
            pn, nm = s53.matched_null_p(w, present, r)
            rows.append(dict(набор=name, окно=a.window, n_рабочий=len(core), слой=layer,
                             r=round(r, 3), p=round(p, 4),
                             p_против_нуля=None if np.isnan(pn) else round(pn, 4),
                             медиана_r_нуля=None if np.isnan(nm) else round(nm, 3)))
            print(f"  {layer:22s} r={r:+.3f} p={p:.4f}  против нуля "
                  f"{'—' if np.isnan(pn) else round(pn, 4)}")
    t = pd.DataFrame(rows)
    f = OUT / f"fits_w{a.window}.tsv"
    t.to_csv(f, sep="\t", index=False)
    print(f"\nwrote {f}")


if __name__ == "__main__":
    main()
