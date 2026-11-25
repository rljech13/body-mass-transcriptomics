#!/usr/bin/env python3
"""Cores of the seven mammal sets, cut the way the bird slices were cut.

A gene stays if its own Pearson r with log10 mass has p < 0.05 and the same
sign as the parent sum, on the paper where the parent passed. A core is drawn
only when at least three genes remain. Panels are not pooled. The home-panel
p of a core is the sum after that selection, not a catalog BH.
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
ROOT = HERE.parent
SCAN = RES / "mammal_size_sets" / "sets_vs_mass.tsv"
MATDIR = RES / "mammal_honest"
OUT = RES / "mammal_size_sets"
PRES = RES / "presentation"
KB = ROOT / "corpus_kleiber_2026-08-06" / "02_база_знаний"
BIRD = RES / "nine_set_slopes" / "nine_set_slopes.tsv"
MIN_CORE = 3

PARENTS = [
    ("ma2016", "REACTOME_BASIGIN_INTERACTIONS", "basigin", "Ma 2016"),
    ("ma2016", "REACTOME_INTERLEUKIN_RECEPTOR_SHC_SIGNALING", "IL receptor SHC", "Ma 2016"),
    ("gse105160", "GOBP_NEGATIVE_REGULATION_OF_GLYCOGEN_METABOLIC_PROCESS", "glycogen, down", "GSE105160"),
    ("gse105160", "GOBP_POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION", "fatty acid oxidation", "GSE105160"),
    ("gse105160", "GOBP_PROTEIN_LOCALIZATION_TO_CILIARY_TRANSITION_ZONE", "ciliary transition zone", "GSE105160"),
    ("ma2016", "GOBP_REGULATION_OF_NEURON_MIGRATION", "neuron migration", "Ma 2016"),
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
    "squirrel monkey": "sq. monkey",
    "pigtail macaque": "pigtail",
    "olive baboon": "baboon",
}


def load40():
    spec = importlib.util.spec_from_file_location("s40", HERE / "40_ruslan_style_scatters.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_panel(name: str) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    wide = pd.read_csv(MATDIR / f"{name}_M.tsv", sep="\t", index_col=0)
    mass = pd.read_csv(MATDIR / f"{name}_mass.tsv", sep="\t").sort_values("mass_g")
    logm = np.log10(mass.mass_g.to_numpy(float))
    return wide[mass.species.tolist()], mass, logm


def gene_fit(wide: pd.DataFrame, logm: np.ndarray, gene: str):
    y = np.log10(wide.loc[gene].to_numpy(float))
    fit = stats.linregress(logm, y)
    return float(fit.rvalue), float(fit.pvalue)


def core_of(wide, logm, genes: list[str], parent_r: float) -> list[str]:
    want = np.sign(parent_r)
    keep = []
    for g in genes:
        r, p = gene_fit(wide, logm, g)
        if np.sign(r) == want and p < 0.05:
            keep.append((p, g))
    keep.sort()
    return [g for _, g in keep]


def sum_rows(wide, mass, genes: list[str]) -> list[tuple] | None:
    present = [g for g in genes if g in wide.index]
    if len(present) < 1:
        return None
    block = wide.loc[present, mass.species.tolist()]
    if (block <= 0).any().any() or block.isna().any().any():
        present = [g for g in present if (wide.loc[g] > 0).all()]
    if not present:
        return None
    totals = wide.loc[present, mass.species.tolist()].sum(axis=0)
    rows = []
    for sp, mg, total in zip(mass.species, mass.mass_g, totals):
        if float(total) <= 0:
            return None
        rows.append((SHORT.get(sp, sp), float(mg), float(total)))
    return rows


def main() -> None:
    s40 = load40()
    scan = pd.read_csv(SCAN, sep="\t")
    panels = {name: load_panel(name) for name in ("ma2016", "wagner", "gse105160")}
    cores = []
    for panel, sid, label, paper in PARENTS:
        rec = scan[(scan.panel == panel) & (scan.id == sid)].iloc[0]
        wide, mass, logm = panels[panel]
        genes = [g for g in str(rec.genes).split(";") if g]
        missing = [g for g in genes if g not in wide.index]
        if missing:
            raise SystemExit(f"{sid} missing {missing}")
        kept = core_of(wide, logm, genes, float(rec.r))
        cores.append(dict(
            panel=panel, id=sid, label=label, paper=paper,
            n_parent=len(genes), parent_r=float(rec.r), genes=kept,
        ))
        print(f"{label:28s} {len(genes):2d} -> {len(kept):2d}  {', '.join(kept)}")

    drawn = [c for c in cores if len(c["genes"]) >= MIN_CORE]
    short = [c for c in cores if len(c["genes"]) < MIN_CORE]
    if len(drawn) != 5 or len(short) != 2:
        raise SystemExit(f"expected 5 cores and 2 short, got {len(drawn)} and {len(short)}")

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.2))
    flat = list(axes.ravel())
    for ax, core in zip(flat, drawn):
        wide, mass, _ = panels[core["panel"]]
        rows = sum_rows(wide, mass, core["genes"])
        if rows is None:
            raise SystemExit(f"no sum for {core['label']}")
        fit = s40.panel(ax, rows, "log10_sum")
        ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
        ax.set_title(
            f"{core['label']}  n={len(core['genes'])}\n{core['paper']}  cor. {fit['r']:.3f}; p = {ptxt}",
            fontsize=9, pad=6,
        )
        core["home_r"] = fit["r"]
        core["home_p"] = fit["p"]
        print(f"  slope {core['paper']:12s} r={fit['r']:+.3f} p={fit['p']:.4g}")
    note = flat[5]
    note.set_axis_off()
    note.text(
        0.5, 0.5,
        "Below the minimum of 3\nneuron migration  n=2\nosmotic stress  n=1",
        ha="center", va="center", fontsize=12, color="#4A4A4A",
    )
    fig.subplots_adjust(wspace=0.38, hspace=0.48, left=0.05, right=0.99, top=0.90, bottom=0.06)
    fig.suptitle("Cores after the same gene cut. Each panel is its own regression.", fontsize=12, y=0.98)
    for folder in (OUT, PRES):
        folder.mkdir(parents=True, exist_ok=True)
        fig.savefig(folder / "fig_mammal_cores.png", dpi=200, bbox_inches="tight", facecolor="white")
        fig.savefig(folder / "fig_mammal_cores.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)

    cross = []
    for core in drawn:
        for panel in ("ma2016", "wagner", "gse105160"):
            wide, mass, _ = panels[panel]
            rows = sum_rows(wide, mass, core["genes"])
            used = [g for g in core["genes"] if g in wide.index and (wide.loc[g] > 0).all()]
            if rows is None:
                r = p = np.nan
                n = 0
            else:
                fit = stats.pearsonr(np.log10([a[1] for a in rows]), np.log10([a[2] for a in rows]))
                r, p, n = float(fit.statistic), float(fit.pvalue), len(used)
            cross.append(dict(label=core["label"], panel=panel, n=n, r=r, p=p, genes=";".join(core["genes"])))
            print(f"  {core['label']:28s} {panel:12s} n={n:2d} r={r:+.3f} p={p:.4g}" if n else
                  f"  {core['label']:28s} {panel:12s} n=0")
    pd.DataFrame(cross).to_csv(OUT / "mammal_cores_vs_mass.tsv", sep="\t", index=False)

    bird = pd.read_csv(BIRD, sep="\t")
    fib_pass = set(bird.loc[(bird.tissue == "fib") & (bird.bh < 0.05), "set"])
    lists = {}
    for path in (
        KB / "working_sets_disjoint_2026-09-22.tsv",
        KB / "working_sets_overlap_ok_2026-09-22.tsv",
    ):
        df = pd.read_csv(path, sep="\t")
        for sid, g in df.groupby("set_id"):
            if sid in fib_pass:
                lists[sid] = set(g["gene"].astype(str))
    print("bird fibroblast BH<0.05:", sorted(fib_pass))
    core_sets = {c["label"]: set(c["genes"]) for c in drawn}
    print("\ncore vs core")
    labels = list(core_sets)
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            hit = core_sets[a] & core_sets[b]
            if hit:
                print(f"  {a} & {b}: {sorted(hit)}")
    print("  none" if not any(core_sets[a] & core_sets[b] for i, a in enumerate(labels) for b in labels[i+1:]) else "")
    print("\ncore vs bird fibroblast")
    any_bird = False
    for label, genes in core_sets.items():
        for sid, bg in lists.items():
            hit = genes & bg
            if hit:
                any_bird = True
                print(f"  {label} & {sid}: {sorted(hit)}")
    if not any_bird:
        print("  none")
    print(f"wrote {PRES / 'fig_mammal_cores.png'}")


if __name__ == "__main__":
    main()
