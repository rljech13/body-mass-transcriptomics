#!/usr/bin/env python3
"""Nine locked sets, fibroblast and epiblast, Ruslan-style log-log.

Membership is the frozen 22.09 files. Score is the one in
Рабочие_наборы_без_пересечений: promoter CAGE already in
M_fib_nearest2kb / M_epi_nearest2kb, missing gene dropped, y = log10(sum TPM),
x = log10(adult mass). One grid per tissue. Sets are not pooled.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES, SPECIES, SPECIES_ORDER, normalize_symbol  # noqa: E402

KB = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/corpus_kleiber_2026-08-06/02_база_знаний")
M_FIB = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
M_EPI = RES / "story_epi_traits" / "M_epi_nearest2kb.tsv"
OUT = RES / "nine_set_slopes"
PRES = RES / "presentation"

# r, p as written in the 22.09 lock. Used only as a check, not as plotted numbers.
EXPECT = {
    "oxphos_16": ((-0.860, 0.013), (-0.220, 0.64)),
    "folate_12": ((0.980, 0.0001), (0.903, 0.005)),
    "ribo_5": ((0.899, 0.006), (0.089, 0.85)),
    "green_hubs_7": ((0.845, 0.017), (0.804, 0.029)),
    "green_51_rest": ((0.926, 0.003), (0.727, 0.064)),
    "carbon_7": ((0.921, 0.003), (0.962, 0.0005)),
    "amino_5": ((0.921, 0.003), (0.964, 0.0005)),
    "nafld_4": ((-0.966, 0.0004), (-0.565, 0.19)),
    "epi_lipid_3": ((0.620, 0.14), (0.969, 0.0003)),
}
ORDER = list(EXPECT)


def load40():
    spec = importlib.util.spec_from_file_location(
        "s40", Path(__file__).resolve().parent / "40_ruslan_style_scatters.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_sets() -> dict[str, tuple[str, list[str]]]:
    frames = [
        pd.read_csv(KB / "working_sets_disjoint_2026-09-22.tsv", sep="\t"),
        pd.read_csv(KB / "working_sets_overlap_ok_2026-09-22.tsv", sep="\t"),
    ]
    out = {}
    for df in frames:
        for sid, g in df.groupby("set_id", sort=False):
            genes = [normalize_symbol(x) for x in g["gene"]]
            genes = [x for x in genes if x]
            label = str(g["set_label"].iloc[0])
            out[str(sid)] = (label, genes)
    missing = [s for s in ORDER if s not in out]
    if missing:
        raise SystemExit(f"sets missing from the lock files: {missing}")
    return out


def rows_of(wide: pd.DataFrame, genes: list[str]) -> tuple[list[tuple], list[str]]:
    hit = [g for g in genes if g in wide.index]
    dropped = [g for g in genes if g not in wide.index]
    sums = wide.loc[hit].sum(axis=0)
    rows = []
    for sp in SPECIES_ORDER:
        rows.append((sp.lower(), SPECIES[sp]["mass_g"], float(sums[sp])))
    return rows, dropped


def main() -> None:
    s40 = load40()
    sets = load_sets()
    fib = pd.read_csv(M_FIB, sep="\t", index_col=0)[SPECIES_ORDER]
    epi = pd.read_csv(M_EPI, sep="\t", index_col=0)[SPECIES_ORDER]
    OUT.mkdir(parents=True, exist_ok=True)
    recs = []
    panels = {"fib": [], "epi": []}
    for sid in ORDER:
        label, genes = sets[sid]
        for tissue, wide, slot in (("fib", fib, 0), ("epi", epi, 1)):
            rows, dropped = rows_of(wide, genes)
            if any(r[2] <= 0 for r in rows):
                raise SystemExit(f"{sid} {tissue}: non-positive sum, log10 refused")
            x = [r[1] for r in rows]
            y = [r[2] for r in rows]
            fit = stats.pearsonr(np.log10(x), np.log10(y))
            exp_r, exp_p = EXPECT[sid][slot]
            if abs(fit.statistic - exp_r) > 0.015:
                raise SystemExit(
                    f"{sid} {tissue}: r={fit.statistic:+.3f}, lock {exp_r:+.3f}; dropped {dropped}"
                )
            recs.append(dict(
                set=sid, label=label, tissue=tissue, n_listed=len(genes),
                n_used=len(genes) - len(dropped), dropped=",".join(dropped),
                r=fit.statistic, p=fit.pvalue, lock_r=exp_r, lock_p=exp_p,
            ))
            panels[tissue].append((label, len(genes) - len(dropped), rows, float(fit.pvalue)))
            print(f"{sid:16s} {tissue:3s} n={len(genes)-len(dropped):3d} r={fit.statistic:+.3f} p={fit.pvalue:.4g}")

    for tissue in panels:
        ps = np.array([p for _, _, _, p in panels[tissue]], float)
        qs = stats.false_discovery_control(ps, method="bh")
        panels[tissue] = [
            (label, n, rows, q) for (label, n, rows, _), q in zip(panels[tissue], qs)
        ]

    for tissue, ylabel, stem in (
        ("fib", "log10_fib", "fig_nine_fib"),
        ("epi", "log10_epi", "fig_nine_epi"),
    ):
        fig, axes = plt.subplots(3, 3, figsize=(12.6, 10.2))
        for ax, (label, n, rows, q) in zip(axes.ravel(), panels[tissue]):
            fit = s40.panel(ax, rows, ylabel)
            qtxt = f"{q:.3f}" if q >= 0.001 else f"{q:.3g}"
            ax.set_title(f"{label}  n={n}\ncor. {fit['r']:.3f}; BH = {qtxt}", fontsize=9, pad=6)
        fig.subplots_adjust(wspace=0.38, hspace=0.55, left=0.06, right=0.99, top=0.96, bottom=0.05)
        for folder in (OUT, PRES):
            folder.mkdir(parents=True, exist_ok=True)
            fig.savefig(folder / f"{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
            fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
        plt.close(fig)
    table = pd.DataFrame(recs)
    bh_map = {
        (tissue, label): q
        for tissue, items in panels.items()
        for label, _, _, q in items
    }
    table["bh"] = [bh_map[(r.tissue, r.label)] for r in table.itertuples()]
    table.to_csv(OUT / "nine_set_slopes.tsv", sep="\t", index=False)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
