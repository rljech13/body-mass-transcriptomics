#!/usr/bin/env python3
"""Reproduce the colleagues' per-gene trait-correlation heatmap and check its statistics.

The colleagues' `sizes.Rmd` correlates every gene against four species traits
(body mass, longevity, incubation, genome size) with Spearman, then stars the
result. Two things in that block change what the picture means:

1. Multiple testing. The line that applies BH across genes is commented out in
   the Rmd; the active line stars raw p-values.
2. Pseudoreplication. The main block aggregates to `bird_cell` groups, so the
   correlation runs over 15 points (7 species x cell types plus chicken PGC)
   while there are only 7 independent trait values. Every species contributes
   two or three rows carrying the same body mass.

This script reproduces the computation exactly as written, then reruns it with
the correction restored and with one point per species, and reports how many
stars survive each step.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES, normalize_symbol  # noqa: E402

CROSS = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/cross-birds")
OUT = RES / "colleagues_matrix_check"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")

TRAITS = pd.DataFrame({
    "species": ["chicken", "duck", "emu", "peacock", "quail", "turkey", "zebra"],
    "body_mass": [912.5, 1121.0, 36200.0, 4093.75, 96.6, 5811.0, 11.93],
    "incubation": [20.5, 26.5, 56.0, 28.0, 17.0, 28.0, 14.15],
    "longevity": [30.0, 29.1, 16.6, 25.0, 14.58, 13.0, 14.5],
    "genome_size": [1.223, 1.428, 1.594, 1.130, 1.379, 1.281, 1.223],
}).set_index("species")

# file.csv column order, as declared in sizes.Rmd
GROUPS = ([("chicken", "Fibroblast")] * 6 + [("chicken", "PGC")] * 2 + [("chicken", "Epiblast")] * 3
          + [("quail", "Fibroblast")] * 6 + [("quail", "Epiblast")] * 3
          + [("turkey", "Fibroblast")] * 6 + [("turkey", "Epiblast")] * 3
          + [("peacock", "Fibroblast")] * 6 + [("peacock", "Epiblast")] * 4
          + [("duck", "Fibroblast")] * 6 + [("duck", "Epiblast")] * 3
          + [("zebra", "Fibroblast")] * 4 + [("zebra", "Epiblast")] * 3
          + [("emu", "Fibroblast")] * 6 + [("emu", "Epiblast")] * 4)

MT_GENES = """ABCB6 ABCB7 ABCB8 ABCF2 ACADM SDHAF3 AIFM2 AIFM3 ALAS2 ALDH4A1 ATP5F1A ATP5F1B
ATP5F1C ATP5F1D ATP5F1E ATP5PB ATP5MC1 ATP5MC2 ATP5MC3 ATP5PF ATP5PO BAX BCKDHA BCKDHB BCKDK
BCL2 BCS1L BNIP3 TAMM41 CASP7 CASQ1 ADAP2 COX15 COX18 COX6B2 CS CYCS DBT DNAJA3 ETFA ETFB FIS1
GATM GRPEL1 HADHB HCCS HSD3B1 HSD3B2 HTRA2 IMMT CCAR2 MAOB MARS2 MCL1 MFN2 MPV17 MRPL10 MRPL12
MRPL23 MRPL32 MRPL40 MRPL41 MRPL42 MRPL51 MRPL52 MRPL55 MRPS10 MRPS11 MRPS12 MRPS15 MRPS16
MRPS18A MRPS18C MRPS21 MRPS22 MRPS24 MRPS28 MRPS35 KGD4 MSTO1 MTX2 MT-ND1 NDUFA1 NDUFA13 NDUFA2
NDUFA6 NDUFA9 NDUFAB1 NDUFS1 NDUFS2 NDUFS3 NDUFS4 NDUFS7 NDUFS8 NDUFV1 NFS1 NNT NR3C1 OGDH OPA1
OXA1L PHB1 PHB2 PIN4 PITRM1 PMPCA POLG2 PPOX RAB11FIP5 RAF1 RHOT1 RHOT2 SDHA SDHD SLC25A1
SLC25A11 SLC25A15 SLC25A22 SLC25A3 SLC25A6 SUPV3L1 SURF1 TFB2M TIMM10 TIMM13 TIMM17A TIMM17B
TIMM23 TIMM44 TIMM50 TIMM8B TIMM9 TOMM22 TOMM34 UCP3 UQCRB UQCRC1 UQCRH VDAC1 VDAC2 VDAC3""".split()

TRAIT_COLS = ["body_mass", "incubation", "longevity", "genome_size"]


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def spearman_block(mat: pd.DataFrame, traits: pd.DataFrame) -> dict[str, pd.DataFrame]:
    res = {}
    for t in TRAIT_COLS:
        x = traits[t].values
        rows = []
        for vals in mat.values:
            s = stats.spearmanr(vals, x)
            rows.append((float(s.statistic), float(s.pvalue)))
        res[t] = pd.DataFrame(rows, index=mat.index, columns=["rho", "p"])
    return res


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(CROSS / "counts.csv")["0"].map(normalize_symbol).values
    raw = pd.read_csv(CROSS / "file.csv", index_col=0)
    raw.index = genes
    meta = pd.DataFrame(GROUPS, columns=["species", "cell"])
    meta["bird_cell"] = meta["species"] + "_" + meta["cell"]
    raw.columns = meta["bird_cell"].values

    # exactly as the Rmd: mean per bird_cell group, correlate over those groups
    grouped = raw.T.groupby(level=0, sort=False).mean().T
    g_meta = pd.DataFrame({"bird_cell": grouped.columns})
    g_meta["species"] = g_meta["bird_cell"].str.split("_").str[0]
    g_traits = TRAITS.reindex(g_meta["species"].values)
    g_traits.index = grouped.columns
    print(f"as written: {grouped.shape[1]} groups over {g_meta['species'].nunique()} species "
          f"({len(grouped)} genes)", flush=True)

    mt = [g for g in grouped.index if g in set(MT_GENES)]
    print(f"mt_genes present in the matrix: {len(mt)} of {len(set(MT_GENES))}", flush=True)

    report: dict = {"n_groups_as_written": int(grouped.shape[1]),
                    "n_independent_species": 7, "n_mt_genes": len(mt)}

    # ---- 1. as written: 15 groups, raw p ----
    as_written = spearman_block(grouped.loc[mt], g_traits)
    # ---- 2. as written but BH across all genes, as the commented-out line intended ----
    all_genes_p = spearman_block(grouped, g_traits)

    # ---- 3. one point per species (fibroblasts only), raw p and BH ----
    fib_cols = [c for c in grouped.columns if c.endswith("Fibroblast")]
    fib = grouped[fib_cols]
    fib.columns = [c.split("_")[0] for c in fib_cols]
    f_traits = TRAITS.reindex(fib.columns)
    f_traits.index = fib.columns
    per_species = spearman_block(fib.loc[mt], f_traits)
    per_species_all = spearman_block(fib, f_traits)

    rows = []
    for trait in TRAIT_COLS:
        aw = as_written[trait]
        n_raw = int((aw["p"] < 0.05).sum())
        bh_full = bh(all_genes_p[trait]["p"].values)
        bh_mt = pd.Series(bh_full, index=grouped.index).reindex(mt)
        n_bh = int((bh_mt < 0.05).sum())
        ps = per_species[trait]
        n_raw_sp = int((ps["p"] < 0.05).sum())
        bh_sp = pd.Series(bh(per_species_all[trait]["p"].values), index=fib.index).reindex(mt)
        n_bh_sp = int((bh_sp < 0.05).sum())
        rows.append(dict(trait=trait, n_mt_genes=len(mt),
                         stars_as_written=n_raw,
                         pct_as_written=round(100 * n_raw / len(mt), 1),
                         after_BH_as_written=n_bh,
                         stars_one_point_per_species=n_raw_sp,
                         after_BH_one_point_per_species=n_bh_sp,
                         min_p_as_written=float(aw["p"].min()),
                         min_p_per_species=float(ps["p"].min())))
        print(f"{trait:12s}: stars as written {n_raw:3d}/{len(mt)} ({100*n_raw/len(mt):4.1f}%), "
              f"after BH {n_bh:3d}; one point per species {n_raw_sp:3d}, after BH {n_bh_sp:3d}")
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "colleagues_heatmap_reproduction.tsv", sep="\t", index=False)
    report["per_trait"] = rows

    # how much of the apparent signal is pseudoreplication alone
    n_eff = stats.spearmanr(np.arange(15), np.arange(15)).pvalue
    report["min_attainable_p"] = dict(
        n15=float(stats.spearmanr(np.arange(15), np.arange(15)).pvalue),
        n7=float(stats.spearmanr(np.arange(7), np.arange(7)).pvalue))
    print(f"\nsmallest attainable Spearman p: n=15 -> {report['min_attainable_p']['n15']:.2e}, "
          f"n=7 -> {report['min_attainable_p']['n7']:.2e}")

    # a permutation check: shuffle species labels, count stars the same way
    rng = np.random.default_rng(11)
    perm_counts = []
    sp_order = list(TRAITS.index)
    for _ in range(200):
        perm = {s: t for s, t in zip(sp_order, rng.permutation(sp_order))}
        pt = TRAITS.reindex([perm[s] for s in g_meta["species"].values])
        pt.index = grouped.columns
        x = pt["body_mass"].values
        cnt = 0
        for vals in grouped.loc[mt].values:
            if stats.spearmanr(vals, x).pvalue < 0.05:
                cnt += 1
        perm_counts.append(cnt)
    report["permutation_stars_body_mass"] = dict(
        observed=rows[0]["stars_as_written"], null_median=float(np.median(perm_counts)),
        null_q95=float(np.quantile(perm_counts, 0.95)), n_perm=200)
    print(f"body mass stars: observed {rows[0]['stars_as_written']}, "
          f"label-shuffled null median {np.median(perm_counts):.0f} "
          f"(95th pct {np.quantile(perm_counts, 0.95):.0f})")

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), constrained_layout=True)
    ax = axes[0]
    x = np.arange(len(TRAIT_COLS))
    w = 0.38
    ax.bar(x - w / 2, tab.stars_as_written, w, color="#c1121f", label="as written (15 groups, raw p)")
    ax.bar(x + w / 2, tab.after_BH_as_written, w, color="#457b9d", label="same, BH across genes")
    ax.set_xticks(x)
    ax.set_xticklabels(TRAIT_COLS, rotation=20, ha="right")
    ax.set(ylabel=f"mitochondrial genes with a star (of {len(mt)})",
           title="stars in the colleagues' heatmap\nbefore and after the commented-out correction")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    ax.hist(perm_counts, bins=25, color="#adb5bd", edgecolor="none")
    ax.axvline(rows[0]["stars_as_written"], color="#c1121f", lw=2)
    ax.set(xlabel="mitochondrial genes with raw p<0.05 vs body mass",
           ylabel="species-label shuffles",
           title="observed star count against shuffled species labels")

    fig.suptitle("Reproducing the colleagues' trait-correlation heatmap", fontsize=12)
    for p in (FIG / "fig_colleagues_heatmap_reproduction.png",
              OUT / "fig_colleagues_heatmap_reproduction.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "heatmap_reproduction.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
