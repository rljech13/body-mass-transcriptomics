#!/usr/bin/env python3
"""Ruslan-style log–log scatters for the five locked story sets.

One panel per set. y = log10(sum of gene TPM), x = log10(adult mass).
Honest CAGE, genes present in all seven. Do not glue sets.

Reuses theme / lm / labels from scripts/40_ruslan_style_scatters.py.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import DATA, FIG, RES, SPECIES, SPECIES_ORDER  # noqa: E402

OUT = RES / "story_ruslan_scatters"
M_PATH = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
SETS = [
    ("oxphos_16", "OXPHOS 16", "fig_story_oxphos_16"),
    ("folate_3", "folate 3", "fig_story_folate_3"),
    ("ribo_5", "ribo 5", "fig_story_ribo_5"),
    ("green_hubs_7", "green hubs 7", "fig_story_green_hubs_7"),
    ("green_sig_51", "green 51", "fig_story_green_51"),
]


def load40():
    spec = importlib.util.spec_from_file_location(
        "s40", Path(__file__).resolve().parent / "40_ruslan_style_scatters.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def genes_of(stem: str) -> list[str]:
    return [ln.strip() for ln in (DATA / "genesets" / f"story_{stem}.txt").read_text().splitlines() if ln.strip()]


def rows_of(wide: pd.DataFrame, genes: list[str]) -> list[tuple[str, float, float]]:
    hit = [g for g in genes if g in wide.index]
    if len(hit) != len(genes):
        missing = set(genes) - set(hit)
        raise SystemExit(f"missing from matrix: {missing}")
    sums = wide.loc[hit].sum(axis=0)
    out = []
    for sp in SPECIES_ORDER:
        out.append((sp.lower(), SPECIES[sp]["mass_g"], float(sums[sp])))
    return out


def main() -> None:
    s40 = load40()
    wide = pd.read_csv(M_PATH, sep="\t", index_col=0)[SPECIES_ORDER]
    OUT.mkdir(parents=True, exist_ok=True)
    recs = []
    panels = []
    for stem, title, figstem in SETS:
        genes = genes_of(stem)
        rows = rows_of(wide, genes)
        fit = s40.single(rows, "log10_fib", figstem)
        # also copy into story folder
        for src_dir in (FIG, RES / "oxphos_species_thread" / "ruslan_metric_fib_panels"):
            for ext in (".png", ".pdf"):
                src = src_dir / f"{figstem}{ext}"
                if src.exists():
                    dst = OUT / f"{figstem}{ext}"
                    dst.write_bytes(src.read_bytes())
        recs.append(dict(
            set=stem, title=title, n=len(genes), r=fit["r"], p=fit["p"],
            **{f"sum_{sp.lower()}": rows[i][2] for i, sp in enumerate(SPECIES_ORDER)},
        ))
        panels.append((title, rows, fit))
        print(f"{stem:16s} n={len(genes):3d} r={fit['r']:+.3f} p={fit['p']:.4f}  sums={[round(r[2],1) for r in rows]}")

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(12.4, 7.2))
    axes = axes.ravel()
    for i, (title, rows, fit) in enumerate(panels):
        s40.panel(axes[i], rows, "log10_fib")
        ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
        axes[i].set_title(f"{title}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
    axes[5].set_visible(False)
    fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.06, right=0.99, top=0.94, bottom=0.08)
    s40.save(fig, "fig_story_five_ruslan")
    for ext in (".png", ".pdf"):
        src = FIG / f"fig_story_five_ruslan{ext}"
        if src.exists():
            (OUT / src.name).write_bytes(src.read_bytes())
    plt.close(fig)

    pd.DataFrame(recs).to_csv(OUT / "story_set_sums.tsv", sep="\t", index=False)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
