#!/usr/bin/env python3
"""Is the CAGE fibroblast decline a property of OXPHOS or of highly expressed genes?

Script 26 found that OXPHOS, complex IV, ATP synthase and the cytosolic ribosome
all fall with body mass, while TCA, mitoribosome, proteasome, glycolysis and heat
shock all rise. The falling sets are the most abundant genes in the library, and
the script 25 null matched abundance only at decile resolution. Two tests:

1. Abundance gradient. Split all genes into fine abundance bins and fit the same
   model per bin. If r vs mass is a monotone function of expression level, the
   OXPHOS result is a property of where the set sits in the distribution.
2. Rank test. Replace the level by the within-species percentile rank of the set.
   Ranks are invariant to any monotone rescaling of a species library, so they
   remove every normalisation, denominator and library composition effect.
   If OXPHOS keeps its relative position, the decline is about the shape of the
   distribution, not about OXPHOS losing ground.
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES, SPECIES, SPECIES_ORDER, normalize_symbol, nuclear_oxphos  # noqa: E402

OUT = RES / "explain_cage_decline"
CACHE = RES / "reverify_confounds" / "gene_matrices.pkl"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
N_BINS = 20

LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))
AUGUSTUS = {"Peacock", "Turkey", "Emu"}
IS_AUG = np.array([sp in AUGUSTUS for sp in SPECIES_ORDER])


def fit(y, x=LOGM):
    y, x = np.asarray(y, float), np.asarray(x, float)
    ok = np.isfinite(y) & np.isfinite(x)
    if ok.sum() < 3 or np.ptp(y[ok]) == 0:
        return np.nan, np.nan
    res = stats.linregress(x[ok], y[ok])
    return float(res.rvalue), float(res.pvalue)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with CACHE.open("rb") as fh:
        layers = pickle.load(fh)
    fib, epi = layers["CAGE_fib"], layers["CAGE_epi"]
    oxphos = {normalize_symbol(g) for g in nuclear_oxphos()} - {""}

    wide = pd.DataFrame({sp: fib[sp][0] for sp in SPECIES_ORDER})
    wide_epi = pd.DataFrame({sp: epi[sp][0] for sp in SPECIES_ORDER})
    report: dict = {}

    # ---------- 1. abundance gradient ----------
    core = wide[wide.notna().sum(axis=1) == 7].copy()
    ref = core["Chicken"]
    core["bin"] = pd.qcut(ref.rank(method="first"), N_BINS, labels=False)
    rows = []
    for b, grp in core.groupby("bin"):
        means = grp[SPECIES_ORDER].mean(axis=0).values
        r, p = fit(means)
        rows.append(dict(bin=int(b), n_genes=int(len(grp)),
                         median_chicken_level=float(grp["Chicken"].median()),
                         n_oxphos=int(grp.index.isin(oxphos).sum()), r=r, p=p))
    grad = pd.DataFrame(rows)
    grad.to_csv(OUT / "abundance_gradient.tsv", sep="\t", index=False)
    rho = stats.spearmanr(grad["bin"], grad["r"])
    ox_bins = core.loc[core.index.isin(oxphos), "bin"]
    report["abundance_gradient"] = dict(
        n_core_genes=int(len(core)), n_bins=N_BINS,
        spearman_bin_vs_r=float(rho.statistic), spearman_p=float(rho.pvalue),
        r_in_top_bin=float(grad.iloc[-1]["r"]), r_in_bottom_bin=float(grad.iloc[0]["r"]),
        n_top_bins_with_r_below_minus_0p8=int((grad.tail(5)["r"] < -0.8).sum()),
        oxphos_median_bin=float(ox_bins.median()),
        oxphos_frac_in_top_quartile_bins=float((ox_bins >= N_BINS * 0.75).mean()))
    print("=== 1. abundance gradient (genes detected in all 7 species) ===")
    print(grad.to_string(index=False))
    print(f"spearman(bin, r) = {rho.statistic:+.3f} p={rho.pvalue:.4g}; "
          f"OXPHOS median bin {ox_bins.median():.0f}/{N_BINS-1}, "
          f"{(ox_bins >= N_BINS*0.75).mean():.0%} in top quartile of bins")

    # ---------- 2. within-species rank test ----------
    rank_rows = []
    for label, mat in (("fib", wide), ("epi", wide_epi)):
        pct = mat.rank(axis=0, pct=True)
        present = [g for g in mat.index if g in oxphos]
        ox_pct = pct.loc[present].mean(axis=0).values
        r, p = fit(ox_pct)
        rank_rows.append(dict(layer=label, metric="mean percentile rank of OXPHOS",
                              r=r, p=p, values=np.round(ox_pct, 4).tolist()))
        lvl = mat.loc[present].mean(axis=0).values
        rl, pl = fit(lvl)
        rank_rows.append(dict(layer=label, metric="mean level of OXPHOS",
                              r=rl, p=pl, values=np.round(lvl, 2).tolist()))
    ranks = pd.DataFrame(rank_rows)
    ranks.to_csv(OUT / "rank_test.tsv", sep="\t", index=False)
    report["rank_test"] = rank_rows
    print("\n=== 2. within-species rank test ===")
    for row in rank_rows:
        print(f"  {row['layer']:4s} {row['metric']:32s} r={row['r']:+.3f} p={row['p']:.4g}")
        print(f"        {row['values']}")

    # ---------- 3. shape of the library, top-share ----------
    shape_rows = []
    for q, lab in ((0.99, "top 1%"), (0.95, "top 5%"), (0.90, "top 10%"), (0.50, "top 50%")):
        share = []
        for sp in SPECIES_ORDER:
            v = wide[sp].dropna().values
            thr = np.quantile(v, q)
            share.append(float(v[v >= thr].sum() / v.sum()))
        r, p = fit(share)
        shape_rows.append(dict(quantile=lab, r=r, p=p, values=np.round(share, 4).tolist()))
        print(f"\n{lab:8s} of genes carry share r={r:+.3f} p={p:.4g}  {np.round(share,3).tolist()}")
    report["library_shape"] = shape_rows
    pd.DataFrame(shape_rows).to_csv(OUT / "library_shape.tsv", sep="\t", index=False)

    # ---------- 4. platform decomposition with honest df ----------
    present = np.array([g for g in wide.index if g in oxphos])
    y = wide.loc[present].mean(axis=0).values
    X_full = np.column_stack([np.ones(7), IS_AUG.astype(float), LOGM])
    X_grp = np.column_stack([np.ones(7), IS_AUG.astype(float)])
    def rss(X):
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        return float(((y - X @ beta) ** 2).sum())
    rss_full, rss_grp, rss_null = rss(X_full), rss(X_grp), rss(np.ones((7, 1)))
    f_stat = ((rss_grp - rss_full) / 1) / (rss_full / (7 - 3))
    p_mass_given_group = float(1 - stats.f.cdf(f_stat, 1, 7 - 3))
    report["platform_decomposition"] = dict(
        var_explained_by_platform_alone=float(1 - rss_grp / rss_null),
        var_explained_by_platform_plus_mass=float(1 - rss_full / rss_null),
        f_mass_given_platform=float(f_stat), p_mass_given_platform=p_mass_given_group,
        df_resid=4)
    print(f"\n=== 4. platform vs mass ===")
    print(f"  R2 platform only          {1 - rss_grp/rss_null:.3f}")
    print(f"  R2 platform + mass        {1 - rss_full/rss_null:.3f}")
    print(f"  mass given platform: F={f_stat:.3f} p={p_mass_given_group:.4f} (df 1,4)")

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    ax = axes[0]
    ax.axhline(0, color="#333", lw=1)
    ax.plot(grad["bin"], grad["r"], "o-", color="#457b9d")
    ax.axvspan(ox_bins.quantile(0.25), ox_bins.quantile(0.75), color="#c1121f", alpha=0.12)
    ax.set(xlabel=f"expression bin (0 = lowest, {N_BINS-1} = highest)",
           ylabel="r of bin mean vs log10 mass",
           title=f"1. the fall is a function of expression level\nspearman={rho.statistic:+.2f} "
                 f"p={rho.pvalue:.3g}; shaded = OXPHOS IQR")

    ax = axes[1]
    pct = wide.rank(axis=0, pct=True)
    ox_pct = pct.loc[present].mean(axis=0).values
    lvl = wide.loc[present].mean(axis=0).values
    ax2 = ax.twinx()
    ax.plot(LOGM, lvl, "o-", color="#c1121f", label="mean level")
    ax2.plot(LOGM, ox_pct, "s--", color="#2a9d8f", label="mean percentile rank")
    rl, pl = fit(lvl)
    rr, pr = fit(ox_pct)
    ax.set(xlabel="log10 body mass (g)", ylabel="OXPHOS mean level")
    ax2.set_ylabel("OXPHOS mean percentile rank")
    ax.set_title(f"2. level r={rl:+.2f} p={pl:.3g}\nrank r={rr:+.2f} p={pr:.3g}")
    for sp, x, yv in zip(SPECIES_ORDER, LOGM, lvl):
        ax.annotate(sp, (x, yv), fontsize=7, xytext=(3, 4), textcoords="offset points")

    ax = axes[2]
    for row in shape_rows:
        ax.plot(LOGM, row["values"], "o-", label=f"{row['quantile']} (r={row['r']:+.2f})")
    ax.set(xlabel="log10 body mass (g)", ylabel="share of library",
           title="3. shape of the promoter library")
    ax.legend(fontsize=7, frameon=False)

    fig.suptitle("Is the OXPHOS decline about OXPHOS or about the top of the distribution?", fontsize=11)
    for p in (FIG / "fig_abundance_gradient_rank_test.png", OUT / "fig_abundance_gradient_rank_test.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "abundance_and_rank_summary.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
