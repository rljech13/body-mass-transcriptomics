#!/usr/bin/env python3
"""Wagner / Ba skin fibroblasts on official KEGG 00190 and 04714.

Same score as scripts/42: 1e6 × set / library, missing gene = 0.
Same swapped-axis figure. Not pooled with Ma or GSE105160.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


HERE = Path(__file__).resolve().parent
s42 = _load("s42", HERE / "42_file_panels_ruslan.py")
s37 = _load("s37", HERE / "37_wagner_sf_oxphos.py")


def score_wagner(sets: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    wanted = {k: s42.wanted_symbols(v) for k, v in sets.items()}
    lib_rows = []
    for sp, sample, _srr in s37.SAMPLES:
        qf = s37.QUANT / sample / "quant.sf"
        if not qf.exists():
            print("MISSING", sample)
            continue
        t2g = pd.read_csv(s37.build_tx2gene(sp), sep="\t")
        q = pd.read_csv(qf, sep="\t")
        q = q.merge(t2g, left_on="Name", right_on="tx", how="left")
        q["symbol"] = q["symbol"].fillna("").map(lambda s: s42.normalize_symbol(s) if s else "")
        if "raw_symbol" in q.columns:
            miss = q["symbol"].eq("")
            q.loc[miss, "symbol"] = q.loc[miss, "raw_symbol"].fillna("").map(s42.normalize_symbol)
        counts = q.groupby("symbol", observed=True)["TPM"].sum()
        counts = counts[counts.index != ""]
        id2sym = {i: i for i in counts.index}
        for set_name, genes in wanted.items():
            sc = s42.score_counts(counts, id2sym, genes)
            lib_rows.append(dict(
                panel="wagner_sf", set=set_name, sample=sample,
                species=sp.replace("_", " "), organism=s37.LATIN[sp],
                mass_g=s37.MASS_G[sp], mass_source="script 37 MASS_G (Wagner panel)",
                **sc,
            ))
    libs = pd.DataFrame(lib_rows)
    sp = libs.groupby(["panel", "set", "species", "organism", "mass_g", "mass_source"], as_index=False).agg(
        n_libs=("sample", "nunique"),
        n_mapped=("n_mapped", "median"),
        n_expr_med=("n_expr", "median"),
        sum_tpm=("sum_tpm", "mean"),
    )
    return libs, sp


def main() -> None:
    s42.OUT.mkdir(parents=True, exist_ok=True)
    sets = s42.load_sets()
    print("sets", {k: len(v) for k, v in sets.items()})
    libs, tab = score_wagner(sets)
    libs.to_csv(s42.OUT / "wagner_libraries.tsv", sep="\t", index=False)
    tab.to_csv(s42.OUT / "wagner_species_scores.tsv", sep="\t", index=False)

    titles = {"kegg_00190": "KEGG hsa00190 OXPHOS", "kegg_04714": "KEGG hsa04714 thermogenesis"}
    for set_name, g in tab.groupby("set"):
        g = g[np.isfinite(g["sum_tpm"]) & (g["sum_tpm"] > 0)]
        fit = s42.fit_xy(g["mass_g"], g["sum_tpm"])
        print(f"FIT wagner_sf {set_name:22s} n={fit['n']} r={fit['r']:+.3f} p={fit['p']:.3g}")

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0))
    for ax, set_name in zip(axes, titles):
        g = tab[tab.set == set_name].sort_values("mass_g")
        fit = s42.draw_panel(ax, g.species.tolist(), g.mass_g.tolist(), g.sum_tpm.tolist())
        ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
        ax.set_title(f"{titles[set_name]}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
    fig.suptitle("Wagner / Ba · skin fibroblasts", fontsize=12, y=0.98)
    fig.subplots_adjust(wspace=0.28, left=0.08, right=0.99, top=0.84, bottom=0.14)
    s42.save_fig(fig, "fig_file_panels_wagner")
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
