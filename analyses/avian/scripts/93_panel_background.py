#!/usr/bin/env python3
"""Where a set's mass slope sits against its own panel.

Each matrix is scored alone. The bird fibroblast matrix is the reference
for the same statistic, not a row in a pooled regression.

For a set: log10(sum) vs log10(mass). The tail is the share of 2000 random
sets of the same size whose correlation is at least as negative. Background
is the median per-gene correlation in that matrix.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import DATA, RES, SPECIES, SPECIES_ORDER, normalize_symbol  # noqa: E402

OUT = RES / "mammal_size_sets"
N_NULL = 2000
RNG = np.random.default_rng(1)
SIXTEEN = [
    "ATP5F1E", "ATP6V0E1", "ATP6V0E2", "COX6A1", "COX6C", "NDUFA5",
    "NDUFA7", "NDUFAB1", "NDUFB2", "NDUFB5", "NDUFB6", "NDUFB8",
    "NDUFC2", "SDHD", "UQCR10", "UQCRC1",
]
KEGG_MT = {
    "ATP6": "MT-ATP6", "ATP8": "MT-ATP8",
    "COX1": "MT-CO1", "COX2": "MT-CO2", "COX3": "MT-CO3",
    "CYTB": "MT-CYB",
    "ND1": "MT-ND1", "ND2": "MT-ND2", "ND3": "MT-ND3",
    "ND4": "MT-ND4", "ND4L": "MT-ND4L", "ND5": "MT-ND5", "ND6": "MT-ND6",
}


def resolve(raw_genes, index: set[str]) -> list[str]:
    hit: list[str] = []
    seen: set[str] = set()
    for g in raw_genes:
        raw = str(g).upper().strip()
        norm = normalize_symbol(raw)
        for a in (norm, raw, KEGG_MT.get(raw, ""), KEGG_MT.get(norm, "")):
            if a and a in index and a not in seen:
                hit.append(a)
                seen.add(a)
                break
    return hit


def gene_r(wide: pd.DataFrame, logm: np.ndarray) -> np.ndarray:
    y = np.log10(wide.to_numpy(float))
    yc = y - y.mean(axis=1, keepdims=True)
    xc = logm - logm.mean()
    return (yc @ xc) / np.sqrt((yc * yc).sum(axis=1) * float(xc @ xc))


def sum_fit(wide: pd.DataFrame, genes: list[str], logm: np.ndarray) -> tuple[float, float, float]:
    y = np.log10(wide.loc[genes].to_numpy(float).sum(axis=0))
    fit = stats.linregress(logm, y)
    return float(fit.rvalue), float(fit.pvalue), float(fit.slope)


def left_tail(wide: pd.DataFrame, n: int, logm: np.ndarray, obs: float) -> tuple[float, float]:
    mat = wide.to_numpy(float)
    null = np.empty(N_NULL)
    for i in range(N_NULL):
        take = RNG.choice(mat.shape[0], size=n, replace=False)
        y = np.log10(mat[take].sum(axis=0))
        null[i] = stats.linregress(logm, y).rvalue
    return float((null <= obs).mean()), float(np.median(null))


def main() -> None:
    kegg = pd.read_csv(DATA / "genesets" / "kegg_hsa_all_pathways.tsv", sep="\t")
    ox_raw = str(kegg.loc[kegg.pathway == "hsa00190", "genes"].iloc[0]).split(";")
    panels: list[tuple[str, pd.DataFrame, np.ndarray]] = []
    for name in ("ma2016", "wagner", "gse105160"):
        wide = pd.read_csv(RES / "mammal_honest" / f"{name}_M.tsv", sep="\t", index_col=0)
        mass = pd.read_csv(RES / "mammal_honest" / f"{name}_mass.tsv", sep="\t")
        wide = wide[mass.species.tolist()]
        logm = np.log10(mass.mass_g.to_numpy(float))
        panels.append((name, wide, logm))
    bird = pd.read_csv(RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv", sep="\t", index_col=0)
    bird = bird[SPECIES_ORDER]
    bird = bird.loc[bird.notna().all(axis=1) & (bird > 0).all(axis=1)]
    blogm = np.log10(np.array([SPECIES[s]["mass_g"] for s in SPECIES_ORDER], float))
    panels.append(("birds_fib_cage", bird, blogm))

    rows = []
    for name, wide, logm in panels:
        r = gene_r(wide, logm)
        index = set(wide.index)
        base = dict(
            panel=name, n_genes=int(wide.shape[0]), n_species=int(wide.shape[1]),
            background_median_r=float(np.median(r)),
            background_frac_neg=float((r < 0).mean()),
        )
        for set_name, raw in (("hsa00190", ox_raw), ("bird16", SIXTEEN)):
            genes = resolve(raw, index)
            rec = dict(base)
            rec["set"] = set_name
            rec["n_set"] = len(genes)
            if len(genes) < 8:
                rec.update(r=np.nan, p=np.nan, slope=np.nan, left_tail=np.nan, null_median_r=np.nan)
                rec["genes_present"] = ";".join(genes)
            else:
                rr, pp, sl = sum_fit(wide, genes, logm)
                tail, nmed = left_tail(wide, len(genes), logm, rr)
                rec.update(r=rr, p=pp, slope=sl, left_tail=tail, null_median_r=nmed, genes_present="")
            rows.append(rec)
            print(f"{name:16} {set_name:8} n={rec['n_set']:3d} r={rec['r']} tail={rec['left_tail']}", flush=True)
    tab = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    tab.to_csv(OUT / "background_vs_oxphos.tsv", sep="\t", index=False)
    print(f"wrote {OUT / 'background_vs_oxphos.tsv'}", flush=True)


if __name__ == "__main__":
    main()
