#!/usr/bin/env python3
"""KEGG set scan in Ruslan's metric: official list, missing gene×species = 0.

Does not freeze membership to genes detected in all seven. log10(sum TPM)
vs log10 mass. Also r of how many map genes are detected — if that tracks
the set score, the decline is coverage, not level.

Universe: honest nearest-CDS matrix, all genes with a number in at least
one species. Genes on a KEGG map but never in the matrix are 0 everywhere
and do not move r.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import DATA, RES, SPECIES, SPECIES_ORDER, load_geneset, normalize_symbol

OUT = RES / "honest_set_scan"
M_PATH = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
KEGG = DATA / "genesets" / "kegg_hsa_all_pathways.tsv"
LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in SPECIES_ORDER], float))
YC = LOGM - LOGM.mean()
YSS = float(YC @ YC)
N_NULL = 2000
RNG = np.random.default_rng(20260921)
MIN_N = 8
SIXTEEN = [
    "ATP5F1E", "ATP6V0E1", "ATP6V0E2", "COX6A1", "COX6C", "NDUFA5",
    "NDUFA7", "NDUFAB1", "NDUFB2", "NDUFB5", "NDUFB6", "NDUFB8",
    "NDUFC2", "SDHD", "UQCR10", "UQCRC1",
]


def r_of(y: np.ndarray) -> float:
    z = np.asarray(y, float) - np.mean(y)
    den = np.sqrt(float(z @ z) * YSS)
    return float((z @ YC) / den) if den else np.nan


def fit(y) -> tuple[float, float]:
    a = stats.linregress(LOGM, np.asarray(y, float))
    return float(a.rvalue), float(a.pvalue)


def kegg_bin(pid: str) -> str:
    try:
        n = int(str(pid)[3:])
    except ValueError:
        return "named"
    if n < 2000:
        return "метаболизм"
    if n < 3500:
        return "генетика"
    if n < 4200:
        return "клетка"
    if n < 5000:
        return "организм"
    return "болезнь"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(M_PATH, sep="\t", index_col=0)[SPECIES_ORDER]
    mat = raw.to_numpy(float)
    present = np.isfinite(mat) & (mat > 0)
    mat0 = np.nan_to_num(mat, nan=0.0)
    genes = list(raw.index)
    pos = {g: i for i, g in enumerate(genes)}
    print(f"матрица {mat.shape[0]} генов (хотя бы один вид), нули на дырах", flush=True)

    catalog: list[tuple[str, str, list[str]]] = [
        ("set16", "16 packed OXPHOS", list(SIXTEEN)),
        ("ruslan_133", "его OXPHOS 133", sorted({normalize_symbol(g) for g in load_geneset("ruslan_oxphos_133")} - {""})),
    ]
    kegg = pd.read_csv(KEGG, sep="\t")
    for rec in kegg.itertuples(index=False):
        gs = [normalize_symbol(x) for x in str(rec.genes).split(";") if x.strip()]
        catalog.append((rec.pathway, rec.name, [g for g in gs if g]))

    rows = []
    for sid, label, glist in catalog:
        idx = np.array([pos[g] for g in glist if g in pos], dtype=int)
        if len(idx) < MIN_N:
            continue
        sums = mat0[idx].sum(axis=0)
        if np.any(sums <= 0):
            continue
        r, p = fit(np.log10(sums))
        nd = present[idx].sum(axis=0).astype(float)
        rn, pn = fit(nd)
        rows.append(dict(
            id=sid, набор=label, bin=kegg_bin(sid),
            n_в_матрице=int(len(idx)), n_в_карте=len(glist),
            r=r, p=p, r_обнаружения=rn, p_обнаружения=pn,
            **{f"n_{sp}": int(nd[i]) for i, sp in enumerate(SPECIES_ORDER)},
            **{f"sum_{sp}": float(sums[i]) for i, sp in enumerate(SPECIES_ORDER)},
        ))

    t = pd.DataFrame(rows)

    # null: random n genes from the matrix, zeros already in mat0
    sizes = sorted({int(n) for n in t["n_в_матрице"].unique()})
    null_r: dict[int, np.ndarray] = {}
    G = mat0.shape[0]
    for n in sizes:
        rs = np.empty(N_NULL)
        for i in range(N_NULL):
            pick = RNG.choice(G, size=n, replace=False)
            rs[i] = r_of(np.log10(mat0[pick].sum(axis=0).clip(min=1e-12)))
        null_r[n] = rs
        print(f"  нуль n={n} готов", flush=True)

    t["p_нуля"] = [float((np.abs(null_r[int(n)]) >= abs(r)).mean())
                   for n, r in zip(t["n_в_матрице"], t["r"])]
    kegg_t = t[t.id.str.startswith("hsa", na=False)].copy()
    kegg_t["bh"] = stats.false_discovery_control(kegg_t["p_нуля"].fillna(1.0).to_numpy(), method="bh")
    t = t.merge(kegg_t[["id", "bh"]], on="id", how="left")
    t = t.sort_values("r")
    t.to_csv(OUT / "kegg_zeros_vs_mass.tsv", sep="\t", index=False)

    print(f"\nкарт {len(kegg_t)}  BH<0.05 {(kegg_t.bh < 0.05).sum()}")
    print("\n16 / 133 / 00190 / 04714 / фолат / терпены:")
    for sid in ("set16", "ruslan_133", "hsa00190", "hsa04714", "hsa00670", "hsa00900", "hsa00020", "hsa00010"):
        sub = t[t.id == sid]
        if sub.empty:
            continue
        r = sub.iloc[0]
        print(f"  {sid:12s} n={int(r.n_в_матрице):3d} r={r.r:+.3f} нуль {r.p_нуля:.4f}  "
              f"r_обнар {r.r_обнаружения:+.3f}  n={int(r.n_Zebra)}…{int(r.n_Emu)}")

    print("\nметаболизм, 10 вниз / 10 вверх:")
    met = t[t.bin == "метаболизм"]
    for rec in met.nsmallest(10, "r").itertuples():
        print(f"  ↓ {rec.id:10s} n={rec.n_в_матрице:3d} r={rec.r:+.3f} обнар {rec.r_обнаружения:+.3f}  {str(rec.набор)[:40]}")
    for rec in met.nlargest(10, "r").itertuples():
        print(f"  ↑ {rec.id:10s} n={rec.n_в_матрице:3d} r={rec.r:+.3f} обнар {rec.r_обнаружения:+.3f}  {str(rec.набор)[:40]}")

    print("\nBH<0.05 (нули):")
    hit = kegg_t[kegg_t.bh < 0.05].sort_values("r")
    if hit.empty:
        print("  нет")
    else:
        for rec in hit.itertuples():
            print(f"  {rec.id:10s} r={rec.r:+.3f} BH={rec.bh:.4f} обнар {rec.r_обнаружения:+.3f}  {str(rec.набор)[:40]}")
    print(f"wrote {OUT / 'kegg_zeros_vs_mass.tsv'}")


if __name__ == "__main__":
    main()
