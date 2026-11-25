#!/usr/bin/env python3
"""Reproduce the colleagues' 'corr epi fib' figures and find what drives them.

The Dropbox deck shows, for each gene set, log10 of a per-species quantity against
log10 body mass, separately for epiblast and fibroblasts:

  mitochondrial genes  epi r=-0.820 p=0.024   fib r=-0.903 p=0.005
  cell cycle genes     epi r=-0.117 p=0.803   fib r=-0.219 p=0.637

The reported claim is that the mitochondrial correlation holds and looks especially
convincing against the flat cell-cycle control. This script identifies the metric,
reproduces both numbers, and decomposes them into the number of genes of the set
that the species detects and the mean level over those genes.
"""
from __future__ import annotations

import json
import re
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

RMD = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/sizes.Rmd")
CACHE = RES / "reverify_confounds" / "gene_matrices.pkl"
OUT = RES / "colleagues_matrix_check"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")

LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))
TARGET = {"mitochondrial": {"epi": (-0.820, 0.024), "fib": (-0.903, 0.005)},
          "cell cycle": {"epi": (-0.117, 0.803), "fib": (-0.219, 0.637)}}


def fit(y: np.ndarray) -> tuple[float, float]:
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan
    res = stats.linregress(LOGM[ok], y[ok])
    return float(res.rvalue), float(res.pvalue)


def gene_lists() -> dict[str, list[str]]:
    txt = RMD.read_text(errors="ignore")
    out = {}
    for var, key in (("mt_genes", "mitochondrial"), ("cell_cycle_genes", "cell cycle"),
                     ("oxphos_genes", "oxphos (their list)")):
        m = re.search(rf"^{var}\s*<-\s*c\((.*?)\)\s*$", txt, re.S | re.M)
        if not m:
            continue
        # the Rmd leaves many cell-cycle entries unquoted; take every bare token
        toks = re.findall(r"[A-Za-z][A-Za-z0-9\-]*", m.group(1))
        out[key] = sorted({normalize_symbol(t) for t in toks} - {""})
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    import pickle
    with CACHE.open("rb") as fh:
        layers = pickle.load(fh)
    fib, epi = layers["CAGE_fib"], layers["CAGE_epi"]
    wide = {"fib": pd.DataFrame({sp: fib[sp][0] for sp in SPECIES_ORDER}),
            "epi": pd.DataFrame({sp: epi[sp][0] for sp in SPECIES_ORDER})}

    sets = gene_lists()
    sets["OXPHOS (our canon)"] = sorted({normalize_symbol(g) for g in nuclear_oxphos()} - {""})
    for k, v in sets.items():
        print(f"{k}: {len(v)} symbols in the list", flush=True)

    report: dict = {"metrics": {}, "decomposition": {}}
    rows = []
    for name, genes in sets.items():
        for tis in ("fib", "epi"):
            m = wide[tis]
            hit = m.reindex(genes)
            n_det = hit.notna().sum(axis=0).values.astype(float)
            mean_det = hit.mean(axis=0).values
            total = hit.sum(axis=0).values
            for metric, vals in (("log10 sum over detected genes", np.log10(total)),
                                 ("log10 mean over detected genes", np.log10(mean_det)),
                                 ("n detected genes", n_det)):
                r, p = fit(vals)
                rows.append(dict(gene_set=name, tissue=tis, metric=metric, r=r, p=p,
                                 values=np.round(vals, 4).tolist()))
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "setsum_metric_scan.tsv", sep="\t", index=False)

    print("\n=== which metric reproduces the deck ===")
    for name in ("mitochondrial", "cell cycle"):
        for tis in ("epi", "fib"):
            tr, tp = TARGET[name][tis]
            sub = res[(res.gene_set == name) & (res.tissue == tis)]
            best = sub.iloc[(sub["r"] - tr).abs().argsort()].iloc[0]
            print(f"  {name:13s} {tis}: deck r={tr:+.3f} p={tp:.3f}  ->  "
                  f"closest is '{best.metric}' r={best.r:+.3f} p={best.p:.4f}")
            report["metrics"][f"{name}|{tis}"] = dict(
                deck_r=tr, deck_p=tp, best_metric=best.metric,
                reproduced_r=float(best.r), reproduced_p=float(best.p))

    print("\n=== decomposition of log10(sum) = log10(n detected) + log10(mean) ===")
    for name in ("mitochondrial", "cell cycle"):
        for tis in ("fib", "epi"):
            sub = res[(res.gene_set == name) & (res.tissue == tis)].set_index("metric")
            n = np.array(sub.loc["n detected genes", "values"])
            mean = np.array(sub.loc["log10 mean over detected genes", "values"])
            r_s, p_s = fit(np.log10(np.array(sub.loc["log10 sum over detected genes", "values"]) * 0 + 1)
                           + np.log10(n) + mean)
            r_n, p_n = fit(np.log10(n))
            r_m, p_m = fit(mean)
            print(f"  {name:13s} {tis}: sum r={r_s:+.3f} = n_detected r={r_n:+.3f} (p={p_n:.4f}) "
                  f"+ mean r={r_m:+.3f} (p={p_m:.4f});  n = {n.astype(int).tolist()}")
            report["decomposition"][f"{name}|{tis}"] = dict(
                sum_r=r_s, n_detected_r=r_n, n_detected_p=p_n, mean_r=r_m, mean_p=p_m,
                n_detected=n.astype(int).tolist())

    # how much of the set is lost, per set, in the three Augustus species
    print("\n=== detection loss by set (share of the list found) ===")
    loss = []
    for name, genes in sets.items():
        m = wide["fib"].reindex(genes)
        frac = m.notna().sum(axis=0).values / len(genes)
        r, p = fit(frac)
        loss.append(dict(gene_set=name, n_list=len(genes), r=r, p=p,
                         frac=np.round(frac, 3).tolist()))
        print(f"  {name:20s} n={len(genes):4d} found-fraction vs mass r={r:+.3f} p={p:.4f}  "
              f"{np.round(frac,2).tolist()}")
    pd.DataFrame(loss).to_csv(OUT / "detection_loss_by_set.tsv", sep="\t", index=False)
    report["detection_loss"] = loss

    # ---------- figure ----------
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), constrained_layout=True)
    for j, name in enumerate(("mitochondrial", "cell cycle")):
        for i, tis in enumerate(("epi", "fib")):
            ax = axes[j, i]
            sub = res[(res.gene_set == name) & (res.tissue == tis)].set_index("metric")
            y = np.array(sub.loc["log10 sum over detected genes", "values"])
            r, p = fit(y)
            ax.scatter(LOGM, y, s=60, color="#1d3557", zorder=3)
            s = stats.linregress(LOGM, y)
            xs = np.linspace(LOGM.min(), LOGM.max(), 10)
            ax.plot(xs, s.intercept + s.slope * xs, color="#457b9d", lw=2)
            for sp, x, yy in zip(SPECIES_ORDER, LOGM, y):
                ax.annotate(sp, (x, yy), fontsize=7, xytext=(3, 4), textcoords="offset points")
            tr, tp = TARGET[name][tis]
            ax.set(xlabel="log10 body mass (g)", ylabel=f"log10 sum, {tis}",
                   title=f"{name} / {tis}: reproduced r={r:+.3f} p={p:.3f}\n(deck: r={tr:+.3f} p={tp:.3f})")
    fig.suptitle("Reproducing the deck: the plotted quantity is the summed signal of the set", fontsize=12)
    for p in (FIG / "fig_colleagues_setsum_reproduction.png",
              OUT / "fig_colleagues_setsum_reproduction.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
    ax = axes[0]
    for name, col in (("mitochondrial", "#c1121f"), ("cell cycle", "#2a9d8f")):
        sub = res[(res.gene_set == name) & (res.tissue == "fib")].set_index("metric")
        n = np.array(sub.loc["n detected genes", "values"])
        ax.plot(LOGM, n / n.max(), "o-", color=col, label=f"{name} (n list = {len(sets[name])})")
    ax.set(xlabel="log10 body mass (g)", ylabel="genes of the set detected, share of maximum",
           title="how much of each list the species detects")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    for name, col in (("mitochondrial", "#c1121f"), ("cell cycle", "#2a9d8f")):
        sub = res[(res.gene_set == name) & (res.tissue == "fib")].set_index("metric")
        mean = np.array(sub.loc["log10 mean over detected genes", "values"])
        ax.plot(LOGM, mean - mean[0], "o-", color=col, label=name)
    ax.axhline(0, color="#adb5bd", lw=1)
    ax.set(xlabel="log10 body mass (g)", ylabel="log10 mean over detected genes, relative to Zebra",
           title="level per detected gene")
    ax.legend(fontsize=8, frameon=False)
    fig.suptitle("The summed metric splits into gene count and level per gene", fontsize=12)
    for p in (FIG / "fig_colleagues_setsum_decomposition.png",
              OUT / "fig_colleagues_setsum_decomposition.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "setsum_reproduction.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
