#!/usr/bin/env python3
"""Which gene SETS, as a sum, track body mass on the honest CAGE matrix.

The 16 OXPHOS subunits are one set. This asks whether any other a priori
KEGG map, or any WGCNA module, does the same thing: log10(mean TPM of the
genes present in all seven) vs log10 mass, against 2000 abundance-matched
nulls. Empty cells are not filled; a gene not in the 3566 is dropped, not
zeroed.

Universe = honest_gene_scan/M_fib_nearest2kb.tsv, CAGE in all 7.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    RES,
    SPECIES,
    SPECIES_ORDER,
    load_geneset,
    normalize_symbol,
)

OUT = RES / "honest_set_scan"
M_PATH = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
KEGG = DATA / "genesets" / "kegg_hsa_all_pathways.tsv"
MODS = RES / "wgcna_out_v2" / "gene_modules_fib.tsv"
LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in SPECIES_ORDER], float))
N_NULL = 2000
RNG = np.random.default_rng(20260921)
MIN_N = 8
SIXTEEN = [
    "ATP5F1E", "ATP6V0E1", "ATP6V0E2", "COX6A1", "COX6C", "NDUFA5",
    "NDUFA7", "NDUFAB1", "NDUFB2", "NDUFB5", "NDUFB6", "NDUFB8",
    "NDUFC2", "SDHD", "UQCR10", "UQCRC1",
]


YC = LOGM - LOGM.mean()
YSS = float(YC @ YC)


def r_of_y(y: np.ndarray) -> float:
    z = y - y.mean()
    den = np.sqrt(float(z @ z) * YSS)
    return float((z @ YC) / den) if den else np.nan


def fit(y) -> tuple[float, float]:
    r = stats.linregress(LOGM, np.asarray(y, float))
    return float(r.rvalue), float(r.pvalue)


class AbundanceNull:
    """Precompute chicken-abundance deciles; null draws are integer row indices."""

    def __init__(self, wide: pd.DataFrame):
        self.genes = np.asarray(wide.index)
        self.pos = {g: i for i, g in enumerate(self.genes)}
        self.mat = wide.to_numpy(float)
        chick = np.log10(self.mat[:, SPECIES_ORDER.index("Chicken")].clip(1e-12))
        self.dec = pd.qcut(pd.Series(chick).rank(method="first"), 10,
                           labels=False, duplicates="drop").to_numpy()
        self.buckets = [np.flatnonzero(self.dec == k) for k in range(10)]

    def rows_of(self, genes: list[str]) -> np.ndarray:
        return np.array([self.pos[g] for g in genes if g in self.pos], dtype=int)

    def r_mean(self, rows: np.ndarray) -> float:
        return r_of_y(np.log10(self.mat[rows].mean(axis=0)))

    def matched(self, rows: np.ndarray, obs_r: float) -> tuple[float, float]:
        if len(rows) < MIN_N or not np.isfinite(obs_r):
            return np.nan, np.nan
        want = np.bincount(self.dec[rows], minlength=10)
        null = np.empty(N_NULL)
        for i in range(N_NULL):
            parts = []
            for k, need in enumerate(want):
                if need == 0:
                    continue
                av = self.buckets[k]
                parts.append(RNG.choice(av, size=min(int(need), len(av)), replace=False))
            pick = np.concatenate(parts)
            null[i] = self.r_mean(pick)
        null = null[np.isfinite(null)]
        return float((np.abs(null) >= abs(obs_r)).mean()), float(np.median(null))


def score(null: AbundanceNull, genes: list[str]) -> tuple[float, float, int, np.ndarray]:
    rows = null.rows_of(genes)
    y = np.log10(null.mat[rows].mean(axis=0))
    r, p = fit(y)
    return r, p, int(len(rows)), y


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(M_PATH, sep="\t", index_col=0)[SPECIES_ORDER]
    det = raw.notna() & (raw > 0)
    wide = raw.loc[det.all(axis=1)]
    print(f"вселенная {len(wide)} генов, честный CAGE во всех семи", flush=True)
    null = AbundanceNull(wide)

    catalog: list[tuple[str, str, list[str]]] = []
    catalog.append(("set16", "16 packed OXPHOS (найденный набор)", list(SIXTEEN)))
    for name, file in (
        ("hsa00190", "kegg_00190"),
        ("hsa04714", "kegg_04714"),
        ("ruslan_133", "ruslan_oxphos_133"),
        ("TCA", "TCA"),
        ("GLYCOLYSIS", "GLYCOLYSIS"),
        ("MITO_BIOGENESIS", "MITO_BIOGENESIS"),
    ):
        catalog.append((name, name, sorted({normalize_symbol(g) for g in load_geneset(file)} - {""})))
    ox = {normalize_symbol(g) for g in load_geneset("kegg_00190")} - {""}
    th = {normalize_symbol(g) for g in load_geneset("kegg_04714")} - {""}
    catalog.append(("04714_minus_00190", "термогенез без OXPHOS", sorted(th - ox)))

    kegg = pd.read_csv(KEGG, sep="\t")
    for rec in kegg.itertuples(index=False):
        genes = [normalize_symbol(x) for x in str(rec.genes).split(";") if x.strip()]
        genes = [g for g in genes if g]
        catalog.append((rec.pathway, rec.name, genes))

    seen: set[str] = set()
    rows = []
    for sid, label, genes in catalog:
        if sid in seen:
            continue
        seen.add(sid)
        hit = [g for g in genes if g in wide.index]
        if len(hit) < MIN_N:
            continue
        r, p, n, _ = score(null, hit)
        pn, nm = null.matched(null.rows_of(hit), r)
        kind = "kegg" if str(sid).startswith("hsa") else "named"
        rows.append(dict(
            id=sid, набор=label, kind=kind, n=n, r=r, p=p,
            p_нуля=pn, медиана_нуля=nm,
        ))
        if sid in {"set16", "hsa00190", "hsa04714", "04714_minus_00190", "hsa03010", "hsa03050"}:
            print(f"  {sid:22s} n={n:3d} r={r:+.3f} p={p:.4f} нуль {pn:.4f}")

    t = pd.DataFrame(rows)
    kegg_t = t[t.kind == "kegg"].copy()
    kegg_t["bh"] = stats.false_discovery_control(kegg_t["p_нуля"].fillna(1.0).to_numpy(), method="bh")
    t = t.merge(kegg_t[["id", "bh"]], on="id", how="left")
    t = t.sort_values("r")
    t.to_csv(OUT / "sets_vs_mass.tsv", sep="\t", index=False)

    print(f"\nKEGG карт с ядром ≥{MIN_N}: {len(kegg_t)}")
    print(f"BH<0.05 по p нуля: {(kegg_t.bh < 0.05).sum() if 'bh' in kegg_t else 0}")
    print("\n10 самых отрицательных KEGG:")
    show = t[t.kind == "kegg"].nsmallest(10, "r")
    for r in show.itertuples():
        bh = "" if pd.isna(r.bh) else f" BH={r.bh:.3f}"
        print(f"  {r.id:12s} n={r.n:3d} r={r.r:+.3f} нуль {r.p_нуля:.4f}{bh}  {r.набор[:50]}")
    print("\n10 самых положительных KEGG:")
    show = t[t.kind == "kegg"].nlargest(10, "r")
    for r in show.itertuples():
        bh = "" if pd.isna(r.bh) else f" BH={r.bh:.3f}"
        print(f"  {r.id:12s} n={r.n:3d} r={r.r:+.3f} нуль {r.p_нуля:.4f}{bh}  {r.набор[:50]}")

    if MODS.exists():
        gm = pd.read_csv(MODS, sep="\t")
        mrows = []
        for mod, g in gm.groupby("module"):
            if mod == "grey":
                continue
            hit = [x for x in g["gene"] if x in wide.index]
            if len(hit) < MIN_N:
                continue
            r, p, n, _ = score(null, hit)
            pn, nm = null.matched(null.rows_of(hit), r)
            mrows.append(dict(модуль=mod, n=n, r=r, p=p, p_нуля=pn, медиана_нуля=nm))
        mt = pd.DataFrame(mrows).sort_values("r")
        mt["bh"] = stats.false_discovery_control(mt["p_нуля"].fillna(1.0).to_numpy(), method="bh")
        mt.to_csv(OUT / "modules_vs_mass.tsv", sep="\t", index=False)
        print("\nмодули (честный CAGE, членство со старой сети 3691):")
        for r in mt.itertuples():
            if r.p_нуля < 0.05 or abs(r.r) >= 0.7:
                print(f"  {r.модуль:16s} n={r.n:3d} r={r.r:+.3f} нуль {r.p_нуля:.4f} BH={r.bh:.3f}")

    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
