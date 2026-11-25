#!/usr/bin/env python3
"""Ruslan-style log–log scatters: theme_light, lm, species labels, 'cor. ; p ='.

Y is log10 of the OXPHOS *sum in TPM*, not log10 of the library fraction.
birds - ESC.R: fibroblast = mean(colSums(norm.C of oxphos_genes)); y = log10(fibroblast).
RNA fraction f is the same number as sum/1e6, so log10(f) = log10(sum) - 6.
Plotting log10(f) made every RNA panel sit on a negative axis; that was scale,
not negative expression, and not his plotted values.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from scipy import stats

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
FIG = ROOT / "figures"
OUT = ROOT / "results" / "oxphos_species_thread" / "ruslan_metric_fib_panels"

# Birds CAGE: Ruslan list n=133, mean(colSums) on exp_tab (2026-09-04).
# Matches his slide: cor = -0.892, p = 0.007. Units = TPM sum, not mean.
BIRDS_CAGE = [
    ("zebra", 11.93, 29098.4),
    ("quail", 96.6, 18002.8),
    ("chicken", 912.5, 17047.4),
    ("duck", 1121.0, 17494.5),
    ("peacock", 4093.75, 7897.5),
    ("turkey", 5811.0, 13194.1),
    ("emu", 36200.0, 8055.4),
]
# Our RNA recount: nuclear set library fraction × 1e6 = TPM sum.
# Same ranking as log10(frac); axis now matches his log10_fib.
BIRDS_RNA = [
    ("zebra", 11.93, 16649.0),
    ("quail", 96.6, 12127.0),
    ("chicken", 912.5, 10624.0),
    ("duck", 1121.0, 15297.0),
    ("peacock", 4093.75, 8194.0),
    ("turkey", 5811.0, 9684.0),
    ("emu", 36200.0, 6636.0),
]
WAGNER = [
    ("opossum", 104.0, 7373.0),
    ("rat", 300.0, 12512.0),
    ("guinea pig", 728.0, 7501.0),
    ("rabbit", 1800.0, 7043.0),
    ("cat", 4000.0, 8070.0),
    ("rhesus", 8000.0, 13097.0),
    ("dog", 40000.0, 10164.0),
    ("human", 62035.0, 13781.0),
    ("sheep", 80000.0, 12811.0),
    ("pig", 130000.0, 5853.0),
    ("horse", 400000.0, 8464.0),
    ("cow", 750000.0, 13170.0),
]
GSE167580 = [
    ("mouse", 20.5, 12949.0),
    ("chicken", 912.5, 9894.0),
    ("rabbit", 1800.0, 13059.0),
    ("cat", 4000.0, 8383.0),
    ("rhesus", 8000.0, 16221.0),
    ("dog", 40000.0, 8065.0),
    ("human", 62035.0, 12182.0),
    ("sheep", 80000.0, 18170.0),
    ("pig", 130000.0, 8582.0),
    ("cow", 750000.0, 10559.0),
]
GSE167580_MAMMALS = [r for r in GSE167580 if r[0] != "chicken"]

LINE = "#3366FF"
POINT = "#1a1a1a"


def theme_light(ax) -> None:
    ax.set_facecolor("white")
    ax.grid(True, which="major", color="#D0D0D0", lw=0.6)
    ax.grid(True, which="minor", color="#ECECEC", lw=0.4)
    for sp in ax.spines.values():
        sp.set_color("#B0B0B0")
        sp.set_linewidth(0.8)
    ax.tick_params(colors="#4A4A4A", labelsize=8)
    ax.xaxis.label.set_color("#4A4A4A")
    ax.yaxis.label.set_color("#4A4A4A")
    ax.set_axisbelow(True)


def panel(ax, rows, ylabel: str) -> dict:
    names = [r[0] for r in rows]
    x = np.log10(np.array([r[1] for r in rows], float))
    y = np.log10(np.array([r[2] for r in rows], float))
    fit = stats.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 80)
    ax.plot(xs, fit.intercept + fit.slope * xs, color=LINE, lw=1.6, zorder=2)
    ax.scatter(x, y, s=22, c=POINT, zorder=3, linewidths=0)
    texts = [
        ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4)
        for name, xi, yi in zip(names, x, y)
    ]
    adjust_text(
        texts,
        ax=ax,
        x=x,
        y=y,
        expand=(1.15, 1.25),
        force_text=(0.25, 0.4),
        force_static=(0.15, 0.2),
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
    )
    theme_light(ax)
    ax.set_xlabel("log10_body_mass")
    ax.set_ylabel(ylabel or "log10_fib")
    ptxt = f"{fit.pvalue:.3f}" if fit.pvalue >= 0.001 else f"{fit.pvalue:.3g}"
    ax.set_title(f"cor. {fit.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
    return {"r": float(fit.rvalue), "p": float(fit.pvalue), "n": int(len(rows))}


def save(fig, stem: str) -> list[Path]:
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    paths = []
    for folder in (FIG, OUT):
        png = folder / f"{stem}.png"
        pdf = folder / f"{stem}.pdf"
        fig.savefig(png, dpi=200, bbox_inches="tight", facecolor="white")
        fig.savefig(pdf, bbox_inches="tight", facecolor="white")
        paths.extend([png, pdf])
    return paths


def single(rows, ylabel: str, stem: str) -> dict:
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    fit = panel(ax, rows, ylabel)
    fig.subplots_adjust(left=0.16, right=0.96, top=0.88, bottom=0.16)
    save(fig, stem)
    plt.close(fig)
    return fit


def species_means(blocks: list[list[tuple]]) -> list[tuple]:
    acc: dict[str, list] = {}
    mass: dict[str, float] = {}
    for rows in blocks:
        for name, m, frac in rows:
            acc.setdefault(name, []).append(frac)
            mass[name] = m
    return [(n, mass[n], float(np.mean(acc[n]))) for n in mass]


def glued() -> dict:
    all_pts = [("wagner", *r) for r in WAGNER] + [("gse167580", *r) for r in GSE167580_MAMMALS]
    x = np.log10(np.array([r[2] for r in all_pts], float))
    y = np.log10(np.array([r[3] for r in all_pts], float))
    fit_all = stats.linregress(x, y)
    means = species_means([WAGNER, GSE167580_MAMMALS])
    xm = np.log10(np.array([r[1] for r in means], float))
    ym = np.log10(np.array([r[2] for r in means], float))
    fit_mean = stats.linregress(xm, ym)

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    xs = np.linspace(x.min(), x.max(), 80)
    ax.plot(xs, fit_all.intercept + fit_all.slope * xs, color=LINE, lw=1.6, zorder=2)
    styles = {
        "wagner": dict(c=POINT, marker="o", s=28),
        "gse167580": dict(c="white", marker="s", s=26, edgecolors=POINT, linewidths=1.1),
    }
    texts = []
    for study, name, mass, frac in all_pts:
        xi, yi = np.log10(mass), np.log10(frac)
        ax.scatter([xi], [yi], zorder=3, **styles[study])
        texts.append(ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4))
    ax.scatter([], [], c=POINT, marker="o", s=28, label="Wagner SF")
    ax.scatter([], [], c="white", marker="s", s=26, edgecolors=POINT, linewidths=1.1, label="GSE167580")
    adjust_text(
        texts,
        ax=ax,
        x=x,
        y=y,
        expand=(1.12, 1.2),
        force_text=(0.22, 0.35),
        force_static=(0.12, 0.18),
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
    )
    theme_light(ax)
    ax.set_xlabel("log10_body_mass")
    ax.set_ylabel("log10_fib")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ptxt = f"{fit_all.pvalue:.3f}" if fit_all.pvalue >= 0.001 else f"{fit_all.pvalue:.3g}"
    ax.set_title(f"cor. {fit_all.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
    fig.subplots_adjust(left=0.15, right=0.97, top=0.88, bottom=0.15)
    save(fig, "fig_ruslan_style_wagner_gse")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    panel(ax, means, "log10_fib")
    fig.subplots_adjust(left=0.16, right=0.96, top=0.88, bottom=0.16)
    save(fig, "fig_ruslan_style_wagner_gse_mean")
    plt.close(fig)

    return {
        "all_n": int(len(all_pts)),
        "all_r": float(fit_all.rvalue),
        "all_p": float(fit_all.pvalue),
        "mean_n": int(len(means)),
        "mean_r": float(fit_mean.rvalue),
        "mean_p": float(fit_mean.pvalue),
        "means": means,
    }


def main() -> None:
    # glued() used to merge Wagner SF and GSE167580 into one regression. That is
    # forbidden by the project rule (panels of different studies stay separate)
    # and the four files it wrote have been deleted.
    w = single(WAGNER, "log10_fib", "fig_ruslan_style_wagner")
    g = single(GSE167580_MAMMALS, "log10_fib", "fig_ruslan_style_gse167580")
    print("Wagner standalone", w)
    print("GSE167580 mammals-only standalone", g)

    # Two-up like his "corr epi fib - OXPHOS.pdf" (6 x 3 in).
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    wt = panel(axes[0], WAGNER, "log10_fib")
    gt = panel(axes[1], GSE167580_MAMMALS, "log10_fib")
    axes[0].set_title(f"Wagner SF   cor. {wt['r']:.3f}; p = {wt['p']:.3f}", fontsize=10, pad=8)
    axes[1].set_title(f"GSE167580 mammals   cor. {gt['r']:.3f}; p = {gt['p']:.3f}", fontsize=10, pad=8)
    fig.subplots_adjust(wspace=0.28, left=0.08, right=0.99, top=0.86, bottom=0.16)
    save(fig, "fig_ruslan_style_mammals")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.4))
    c = panel(axes[0], BIRDS_CAGE, "log10_fib")
    r = panel(axes[1], BIRDS_RNA, "log10_fib")
    axes[0].set_title(f"GSE213253 CAGE   cor. {c['r']:.3f}; p = {c['p']:.3f}", fontsize=10, pad=8)
    axes[1].set_title(f"GSE213253 RNA   cor. {r['r']:.3f}; p = {r['p']:.3f}", fontsize=10, pad=8)
    fig.subplots_adjust(wspace=0.28, left=0.08, right=0.99, top=0.86, bottom=0.16)
    save(fig, "fig_ruslan_style_birds")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.6))
    specs = [
        (axes[0, 0], BIRDS_CAGE, "log10_fib", "birds CAGE"),
        (axes[0, 1], BIRDS_RNA, "log10_fib", "birds RNA"),
        (axes[1, 0], WAGNER, "log10_fib", "Wagner SF"),
        (axes[1, 1], GSE167580_MAMMALS, "log10_fib", "GSE167580 mammals"),
    ]
    for ax, rows, ylab, name in specs:
        fit = panel(ax, rows, ylab)
        ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
        ax.set_title(f"{name}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
    fig.subplots_adjust(wspace=0.32, hspace=0.38, left=0.08, right=0.99, top=0.95, bottom=0.08)
    save(fig, "fig_ruslan_style_four")
    plt.close(fig)

    print("Wagner two-up", wt)
    print("GSE167580 mammals two-up", gt)
    print("birds CAGE log10(sum TPM)~log10(mass)", c)
    print("birds RNA log10(sum TPM)~log10(mass)", r)


if __name__ == "__main__":
    main()
