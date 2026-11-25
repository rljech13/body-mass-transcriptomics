#!/usr/bin/env python3
"""Gene x library matrices for WGCNA, plus the widened gene universe.

Everything in this project so far collapsed libraries to a per-tissue mean. WGCNA
needs the libraries themselves, so the same measurement rule is applied per
library here:

    gene 5' end  <- shipped annotation, then reanno.gff (whole-proteome
                    re-annotation, script 58), then the two targeted recovery
                    rounds (scripts 50 and 56)
    value(gene, library) = sum over peaks within WINDOW bp of that 5' end on the
                    gene's strand, each peak first thresholded at 1 TPM
    stored as log10(value + 1)   -- a log, not a VST: species are compared here

Universe, widened as planned: the old 2828 required non-zero signal in BOTH
tissues in all seven species. Here a gene enters the fibroblast universe if it is
non-zero in at least one fibroblast library in every species, and likewise for
epiblast; the combined universe is the union.

Networks are built per tissue across species, which is what removes the
tissue-contrast bias without any correction:
    M_fib  42 libraries   M_epi  23 libraries   M_all  65 libraries
plus one matrix per species (7-11 libraries) as a sensitivity.
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
s57 = importlib.import_module("57_working_sets_round2")
from oxphos_thread_lib import (  # noqa: E402
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    peak_to_symbol_de,
)

EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
OUT = RES / "wgcna_input"
THR = 1.0


def five_prime(species: str, d: pd.DataFrame) -> dict[str, dict]:
    """Same locus priority as script 57, exposed here so both stay in step."""
    fp: dict[str, dict] = {}
    ann = d.dropna(subset=["geneStart", "geneEnd"])
    ann = ann[ann["symbol"].ne("")].copy()
    ann["sg"] = ann["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(ann["strand"])
    for sym, g in ann.groupby("symbol", observed=True):
        r = g.iloc[0]
        st = r["sg"]
        fp[sym] = dict(seqname=r["seqnames"], strand=st,
                       five=float(r["geneStart"] if st == "+" else r["geneEnd"]))
    pid2sym = dict(re.findall(r"^>(\S+)\s+(\S+)\s+\[organism=", s57.PROTSRC.read_text(), re.M))
    reanno = Path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies")
    for gff, p2s in (
        (reanno / species / "reanno.gff", None),
        (s57.R1 / f"{species}_miniprot.gff", pid2sym),
        (s57.R2 / f"{species}_q2.gff", None),
    ):
        for sym, c in s57.loci_from_gff(gff, p2s).items():
            fp.setdefault(sym, c)
    return fp


def per_library(species: str, window: int) -> pd.DataFrame:
    d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
    fp = five_prime(species, d)
    cols = lib_cols(d, "fib") + lib_cols(d, "epi")
    pool = d.assign(mid=(d["start"] + d["end"]) / 2.0)
    num = pool[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    num = num.where(num >= THR, 0.0)
    pool = pool.assign(**{c: num[c] for c in cols})
    by = {k: v for k, v in pool.groupby("seqnames", observed=True)}

    rows: dict[str, np.ndarray] = {}
    for sym, c in fp.items():
        sub = by.get(c["seqname"])
        if sub is None:
            continue
        near = sub[(sub["strand"] == c["strand"]) & ((sub["mid"] - c["five"]).abs() <= window)]
        if near.empty:
            continue
        v = near[cols].to_numpy(float).sum(axis=0)
        if v.sum() <= 0:
            continue
        rows[sym] = v
    m = pd.DataFrame.from_dict(rows, orient="index", columns=cols)
    m.columns = [f"{species}|{'fib' if c.startswith('norm.C') else 'epi'}|{c.split('.')[1]}" for c in m.columns]
    return np.log10(m + 1.0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=1000)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    per: dict[str, pd.DataFrame] = {}
    for sp in SPECIES_ORDER:
        per[sp] = per_library(sp, a.window)
        f = OUT / f"M_{sp}_w{a.window}.tsv"
        per[sp].to_csv(f, sep="\t")
        nf = sum(1 for c in per[sp].columns if "|fib|" in c)
        ne = sum(1 for c in per[sp].columns if "|epi|" in c)
        print(f"[{sp:8s}] генов {len(per[sp]):>6}  библиотек {per[sp].shape[1]} (фибр {nf}, эпи {ne})")

    def universe(tissue: str) -> list[str]:
        sets = []
        for sp in SPECIES_ORDER:
            c = [x for x in per[sp].columns if f"|{tissue}|" in x]
            m = per[sp][c]
            sets.append(set(m.index[(m > 0).any(axis=1)]))
        return sorted(set.intersection(*sets))

    uf, ue = universe("fib"), universe("epi")
    ua = sorted(set(uf) | set(ue))
    print(f"\nвселенная: фибробласты {len(uf)}, эпибласт {len(ue)}, объединение {len(ua)}")
    for nm, u in (("fib", uf), ("epi", ue), ("all", ua)):
        (OUT / f"universe_{nm}_w{a.window}.txt").write_text("\n".join(u) + "\n")

    def combine(tissue: str | None, genes: list[str]) -> pd.DataFrame:
        parts = []
        for sp in SPECIES_ORDER:
            c = [x for x in per[sp].columns if tissue is None or f"|{tissue}|" in x]
            parts.append(per[sp].reindex(genes)[c])
        return pd.concat(parts, axis=1)

    mats = {"fib": combine("fib", uf), "epi": combine("epi", ue), "all": combine(None, ua)}
    for nm, m in mats.items():
        m = m.dropna(how="all")
        m.to_csv(OUT / f"M_{nm}_w{a.window}.tsv", sep="\t")
        print(f"M_{nm}: {m.shape[0]} генов x {m.shape[1]} библиотек, пропусков {int(m.isna().sum().sum())}")

    meta = []
    for c in mats["all"].columns:
        sp, tis, lib = c.split("|")
        meta.append(dict(sample=c, species=sp, tissue=tis, library=lib,
                         mass_g=SPECIES[sp]["mass_g"], log10_mass=np.log10(SPECIES[sp]["mass_g"])))
    pd.DataFrame(meta).to_csv(OUT / f"samples_w{a.window}.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT}/M_*_w{a.window}.tsv, universe_*, samples_w{a.window}.tsv")


if __name__ == "__main__":
    main()
