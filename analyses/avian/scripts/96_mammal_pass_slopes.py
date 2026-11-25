#!/usr/bin/env python3
"""One Ruslan-style grid of mammal sets that passed abundance-null BH.

Three panels stay separate regressions. Wagner has no passing set; that cell
is a note, not a slope. Genes and BH come from sets_vs_mass.tsv.
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
from oxphos_thread_lib import RES  # noqa: E402

HERE = Path(__file__).resolve().parent
SCAN = RES / "mammal_size_sets" / "sets_vs_mass.tsv"
MATDIR = RES / "mammal_honest"
OUT = RES / "mammal_size_sets"
PRES = RES / "presentation"

# Home panel only. Wagner is drawn as an empty cell.
PASS = [
    ("ma2016", "GOBP_REGULATION_OF_NEURON_MIGRATION", "neuron migration", "Ma 2016"),
    ("ma2016", "REACTOME_BASIGIN_INTERACTIONS", "basigin", "Ma 2016"),
    ("ma2016", "REACTOME_INTERLEUKIN_RECEPTOR_SHC_SIGNALING", "IL receptor SHC", "Ma 2016"),
    ("gse105160", "GOBP_NEGATIVE_REGULATION_OF_GLYCOGEN_METABOLIC_PROCESS", "glycogen, down", "GSE105160"),
    ("gse105160", "GOBP_POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION", "fatty acid oxidation", "GSE105160"),
    ("gse105160", "GOBP_PROTEIN_LOCALIZATION_TO_CILIARY_TRANSITION_ZONE", "ciliary transition zone", "GSE105160"),
    ("gse105160", "GOBP_REGULATION_OF_RESPONSE_TO_OSMOTIC_STRESS", "osmotic stress", "GSE105160"),
]

SHORT = {
    "little brown bat": "brown bat",
    "big brown bat": "big bat",
    "white-footed mouse": "white-foot",
    "deer mouse": "deer mouse",
    "meadow vole": "vole",
    "red squirrel": "red sq.",
    "fox squirrel": "fox sq.",
    "guinea pig": "guinea pig",
    "grass rat": "grass rat",
    "cotton rat": "cotton rat",
    "squirrel monkey": "sq. monkey",
    "pigtail macaque": "pigtail",
    "olive baboon": "baboon",
}


def load40():
    spec = importlib.util.spec_from_file_location("s40", HERE / "40_ruslan_style_scatters.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    s40 = load40()
    scan = pd.read_csv(SCAN, sep="\t")
    passed = scan[scan.bh < 0.05]
    if len(passed) != 7:
        raise SystemExit(f"expected 7 BH passes, found {len(passed)}")
    if set(passed.panel) != {"ma2016", "gse105160"}:
        raise SystemExit(f"unexpected panels: {sorted(passed.panel.unique())}")

    mats = {}
    masses = {}
    for panel in ("ma2016", "gse105160"):
        mats[panel] = pd.read_csv(MATDIR / f"{panel}_M.tsv", sep="\t", index_col=0)
        masses[panel] = pd.read_csv(MATDIR / f"{panel}_mass.tsv", sep="\t")

    fig, axes = plt.subplots(2, 4, figsize=(16.8, 9.4))
    flat = list(axes.ravel())
    slots = [0, 1, 2, 4, 5, 6, 7]
    for ax_i, (panel, sid, label, paper) in zip(slots, PASS):
        hit = passed[(passed.panel == panel) & (passed.id == sid)]
        if len(hit) != 1:
            raise SystemExit(f"missing row {panel} {sid}")
        rec = hit.iloc[0]
        genes = [g for g in str(rec.genes).split(";") if g]
        wide = mats[panel]
        missing = [g for g in genes if g not in wide.index]
        if missing:
            raise SystemExit(f"{sid} missing from matrix: {missing}")
        mass = masses[panel].sort_values("mass_g")
        rows = []
        for sp, mg in zip(mass.species, mass.mass_g):
            total = float(wide.loc[genes, sp].sum())
            if total <= 0:
                raise SystemExit(f"{sid} {sp}: non-positive sum")
            rows.append((SHORT.get(sp, sp), float(mg), total))
        x = np.log10([r[1] for r in rows])
        y = np.log10([r[2] for r in rows])
        fit = stats.pearsonr(x, y)
        if abs(fit.statistic - float(rec.r)) > 0.015 or len(genes) != int(rec.n):
            raise SystemExit(
                f"{sid}: r={fit.statistic:+.3f} file {float(rec.r):+.3f} n={len(genes)} file {int(rec.n)}"
            )
        ax = flat[ax_i]
        drawn = s40.panel(ax, rows, "log10_sum")
        extra = "\ndated PGLS 0.083" if sid == "REACTOME_BASIGIN_INTERACTIONS" else ""
        ax.set_title(
            f"{label}  n={len(genes)}\n{paper}  cor. {drawn['r']:.3f}; BH < 0.0005{extra}",
            fontsize=9,
            pad=6,
        )
        print(f"{paper:12s} {label:28s} n={len(genes):2d} r={drawn['r']:+.3f}")

    note = flat[3]
    note.set_axis_off()
    note.text(
        0.5, 0.5,
        "Wagner\nno set with BH < 0.05\n3 086 sets tested",
        ha="center", va="center", fontsize=12, color="#4A4A4A",
    )
    fig.subplots_adjust(wspace=0.42, hspace=0.62, left=0.04, right=0.99, top=0.90, bottom=0.06)
    fig.suptitle("Mammal sums that pass BH. Each panel is its own regression.", fontsize=12, y=0.98)
    for folder in (OUT, PRES):
        folder.mkdir(parents=True, exist_ok=True)
        fig.savefig(folder / "fig_mammal_pass.png", dpi=200, bbox_inches="tight", facecolor="white")
        fig.savefig(folder / "fig_mammal_pass.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {PRES / 'fig_mammal_pass.png'}")


if __name__ == "__main__":
    main()
