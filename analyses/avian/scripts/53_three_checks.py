#!/usr/bin/env python3
"""The three checks the recovered result needs before it can be stated.

A. Matched null. Is the level decline of the recovered core stronger than that of
   random gene sets matched on post-recovery detection breadth and on chicken
   abundance decile? Without this, -0.880 has nothing to be compared against.

B. Symmetric measurement. In script 51 the three Augustus species got their genes
   measured by our window rule while the four native species kept the old
   annotation-based rule. Same gene, two procedures, and the procedures differ
   exactly along the axis under test. Here every gene in every species is measured
   by ONE rule: take the gene's 5' end (from the species' own annotation where it
   exists, from the miniprot locus where it does not), collect peaks within
   1000 bp on the gene's strand, threshold each peak at 1 TPM per library, sum,
   average over the tissue's libraries.

C. The promoter-class filter. Recovered genes bypassed it because that class was
   assigned by the annotation being corrected. Reported both ways: any class, and
   promoter class only, under the symmetric rule.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    load_geneset,
    normalize_symbol,
    peak_to_symbol_de,
)

EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
GFFDIR = RES / "gene_recovery"
PROTSRC = Path("/path/to/storage/kleiber_birds_gemcat/ref/chicken_prot/query_proteins_all.faa")
OUT = RES / "three_checks"
AUG = ["Turkey", "Peacock", "Emu"]
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], float))
WINDOW = 1000
THR = 1.0
N_NULL = 2000
RNG = np.random.default_rng(20260920)
LIST_NAME = "его OXPHOS 133"


def fit(y) -> tuple[float, float]:
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan
    r = stats.linregress(LOGM[ok], y[ok])
    return float(r.rvalue), float(r.pvalue)


def five_prime_table(species: str, peaks: pd.DataFrame, use_miniprot: bool = True) -> pd.DataFrame:
    """symbol -> (seqname, 5' position, strand), one row per symbol.

    Annotation first: every peak row carries the coordinates of the model it was
    assigned to. Then the miniprot loci, for the genes the annotation lacks.
    """
    rows = []
    d = peaks.dropna(subset=["geneStart", "geneEnd"]).copy()
    d = d[d["symbol"].ne("")]
    d["strand_g"] = d["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(d["strand"])
    for sym, g in d.groupby("symbol", observed=True):
        r = g.iloc[0]
        st = r["strand_g"]
        rows.append(
            dict(symbol=sym, seqname=r["seqnames"], strand=st,
                 five=float(r["geneStart"] if st == "+" else r["geneEnd"]), src="annotation")
        )
    have = {r["symbol"] for r in rows}

    gff = GFFDIR / f"{species}_miniprot.gff"
    if use_miniprot and gff.exists():
        pid2sym = dict(re.findall(r"^>(\S+)\s+(\S+)\s+\[organism=", PROTSRC.read_text(), re.M))
        best: dict[str, tuple[float, dict]] = {}
        for line in gff.read_text().splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            at = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
            sym = pid2sym.get(at.get("Target", "").split()[0])
            if not sym or sym in have:
                continue
            seq = f[0].split()[0]
            if seq.startswith("ENA|"):
                seq = seq.split("|")[-1]
            sc = float(f[5]) if f[5] != "." else 0.0
            st = f[6]
            cand = dict(symbol=sym, seqname=seq, strand=st,
                        five=float(f[3] if st == "+" else f[4]), src="miniprot")
            if sym not in best or sc > best[sym][0]:
                best[sym] = (sc, cand)
        rows.extend(v[1] for v in best.values())
    return pd.DataFrame(rows)


def measure(species: str, promoter_only: bool, use_miniprot: bool = True) -> dict[str, pd.Series]:
    """One rule for every species: window around the 5' end."""
    d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
    fp = five_prime_table(species, d, use_miniprot)
    cols = {t: lib_cols(d, t) for t in ("fib", "epi")}
    pool = d[d["class"].astype(str).str.lower().eq("promoter")] if promoter_only else d
    pool = pool.assign(mid=(pool["start"] + pool["end"]) / 2.0)
    by_seq = {k: v for k, v in pool.groupby("seqnames", observed=True)}

    out = {t: {} for t in ("fib", "epi")}
    for _, g in fp.iterrows():
        sub = by_seq.get(g["seqname"])
        if sub is None:
            continue
        near = sub[(sub["strand"] == g["strand"]) & ((sub["mid"] - g["five"]).abs() <= WINDOW)]
        if near.empty:
            continue
        for t in ("fib", "epi"):
            c = cols[t]
            if not c:
                continue
            m = near[c].apply(pd.to_numeric, errors="coerce").fillna(0.0)
            m = m.where(m >= THR, 0.0)
            v = float(m.sum(axis=0).mean())
            if v > 0:
                out[t][g["symbol"]] = v
    return {t: pd.Series(out[t], dtype=float) for t in ("fib", "epi")}


def matched_null_p(wide: pd.DataFrame, core: list[str], obs_r: float) -> tuple[float, float]:
    det = wide.notna() & (wide > 0)
    breadth = det.sum(axis=1)
    chick = wide["Chicken"]
    pool = wide.index[breadth > 0]
    dec = pd.qcut(chick.reindex(pool).rank(method="first"), 10, labels=False, duplicates="drop")
    key = pd.DataFrame({"b": breadth.reindex(pool), "d": dec}, index=pool)
    hit = [g for g in core if g in key.index]
    if len(hit) < 5 or not np.isfinite(obs_r):
        return np.nan, np.nan
    want = key.loc[hit].groupby(["b", "d"]).size()
    buckets = {k: key.index[(key["b"] == k[0]) & (key["d"] == k[1])].to_numpy() for k in want.index}
    null = np.empty(N_NULL)
    for i in range(N_NULL):
        pick: list[str] = []
        for k, need in want.items():
            av = buckets[k]
            if len(av):
                pick.extend(RNG.choice(av, size=min(need, len(av)), replace=False))
        null[i] = fit(np.log10(wide.loc[pick].mean(axis=0).to_numpy(float)))[0]
    null = null[np.isfinite(null)]
    return float((np.abs(null) >= abs(obs_r)).mean()), float(np.median(null))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = sorted({normalize_symbol(g) for g in load_geneset("ruslan_oxphos_133")} - {""})

    rows = []
    for promoter_only in (False, True):
        tag = "только promoter" if promoter_only else "любой класс"
        mats = {sp: measure(sp, promoter_only) for sp in SPECIES_ORDER}
        wide = {t: pd.DataFrame({sp: mats[sp][t] for sp in SPECIES_ORDER}) for t in ("fib", "epi")}
        idx = sorted(set(wide["fib"].index) | set(wide["epi"].index))
        ok = None
        for t in ("fib", "epi"):
            m = wide[t].reindex(idx)
            g = m.notna() & (m > 0)
            ok = g if ok is None else (ok & g)
        breadth = ok.sum(axis=1)
        for k in (7, 6, 5, 4, 3):
            core = [g for g in genes if breadth.get(g, 0) >= k]
            for t in ("fib", "epi"):
                w = wide[t]
                r, p = fit(np.log10(w.reindex(core).mean(axis=0).to_numpy(float)))
                pn, nm = (np.nan, np.nan)
                if k in (7, 5):
                    pn, nm = matched_null_p(w, core, r)
                rows.append(dict(
                    правило=tag, k=k, n_core=len(core), tissue=t,
                    r=round(r, 3), p=round(p, 4),
                    p_matched_null=None if np.isnan(pn) else round(pn, 4),
                    null_median_r=None if np.isnan(nm) else round(nm, 3),
                    n_detected=";".join(str(int(v)) for v in w.reindex(core).notna().sum(axis=0)),
                ))
        print(f"[{tag}] генов измерено: " + ", ".join(f"{sp} {int((wide['fib'][sp] > 0).sum())}" for sp in SPECIES_ORDER))

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "three_checks.tsv", sep="\t", index=False)
    pd.set_option("display.width", 260)
    for tag in tab["правило"].unique():
        print("\n" + "=" * 118)
        print(f"ЕДИНОЕ ПРАВИЛО ДЛЯ ВСЕХ СЕМИ ВИДОВ — {tag}")
        print("=" * 118)
        print(tab[tab["правило"] == tag].drop(columns=["правило"]).to_string(index=False))
    print(f"\nwrote {OUT / 'three_checks.tsv'}")


if __name__ == "__main__":
    main()
