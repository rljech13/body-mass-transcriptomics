#!/usr/bin/env python3
"""A6: CAGE promoter usage for metabolic genes across species and chicken ESC/fib."""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
DATA = ROOT / "data"
RES = ROOT / "results"
FIG = ROOT / "figures"

SYMBOL_SPECIES = {
    "Chicken": "GSE213253_Chicken_CAGE.xlsx",
    "Duck": "GSE213253_Duck_CAGE.xlsx",
    "Quail": "GSE213253_Quail_CAGE.xlsx",
    "ZebraFinch": "GSE213253_ZebraFinch_CAGE.xlsx",
}

MASS = {"ZebraFinch": 0.0157, "Quail": 0.12, "Duck": 1.1, "Chicken": 2.0}


def load_modules() -> dict[str, list[str]]:
    out = {}
    for p in (DATA / "genesets").glob("*.txt"):
        out[p.stem] = [ln.strip().upper() for ln in p.read_text().splitlines() if ln.strip()]
    return out


def promoter_stats(path: Path, genes: set[str]) -> pd.DataFrame:
    df = pd.read_excel(path)
    norms = [c for c in df.columns if str(c).startswith("norm.")]
    df["geneId"] = df["geneId"].astype(str).str.upper()
    sub = df[df["geneId"].isin(genes)].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["_expr"] = sub[norms].apply(pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
    # number of promoters per gene, Shannon diversity of promoter usage, dominant promoter fraction
    rows = []
    for g, gdf in sub.groupby("geneId"):
        vals = gdf["_expr"].astype(float).values
        vals = vals[vals > 0]
        if len(vals) == 0:
            continue
        p = vals / vals.sum()
        shannon = float(-(p * np.log2(p + 1e-12)).sum())
        rows.append(
            {
                "gene": g,
                "n_promoters": int(len(gdf)),
                "n_active_promoters": int((gdf["_expr"] > 0).sum()),
                "total_tpm": float(gdf["_expr"].sum()),
                "dominant_frac": float(p.max()),
                "shannon": shannon,
            }
        )
    return pd.DataFrame(rows)


def main():
    modules = load_modules()
    ox_genes = set(modules.get("OXPHOS", [])) | set(modules.get("TCA", [])) | set(modules.get("GLYCOLYSIS", []))

    # per-species promoter stats for metabolic genes
    all_rows = []
    for sp, fname in SYMBOL_SPECIES.items():
        print("CAGE promoters", sp)
        st = promoter_stats(DATA / fname, ox_genes)
        if st.empty:
            continue
        st["species"] = sp
        st["mass_kg"] = MASS[sp]
        all_rows.append(st)
        st.to_csv(RES / f"cage_promoter_stats_{sp}.tsv", sep="\t", index=False)
    big = pd.concat(all_rows, ignore_index=True)
    big.to_csv(RES / "cage_promoter_stats_metabolic_genes.tsv", sep="\t", index=False)

    # species-level summary: mean n_active promoters / shannon for OXPHOS genes with signal
    summary = (
        big.groupby("species")
        .agg(
            mass_kg=("mass_kg", "first"),
            n_genes=("gene", "nunique"),
            mean_n_active=("n_active_promoters", "mean"),
            mean_shannon=("shannon", "mean"),
            mean_dominant_frac=("dominant_frac", "mean"),
            mean_total_tpm=("total_tpm", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(RES / "cage_promoter_species_summary.tsv", sep="\t", index=False)

    # correlate promoter complexity metrics with mass
    fit_rows = []
    for col in ["mean_n_active", "mean_shannon", "mean_dominant_frac", "mean_total_tpm"]:
        x = np.log10(summary.mass_kg.values)
        y = summary[col].values
        if len(x) >= 3 and np.all(np.isfinite(y)):
            slope, intercept, r, p, se = stats.linregress(x, y)
            fit_rows.append({"metric": col, "slope": slope, "r": r, "p": p, "se": se})
    pd.DataFrame(fit_rows).to_csv(RES / "cage_promoter_vs_mass_fits.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), constrained_layout=True)
    for ax, col, title in zip(
        axes,
        ["mean_n_active", "mean_shannon", "mean_total_tpm"],
        ["Mean active promoters / gene", "Mean promoter Shannon", "Mean total promoter TPM"],
    ):
        x = np.log10(summary.mass_kg)
        y = summary[col]
        ax.scatter(x, y, s=60)
        for sp, xi, yi in zip(summary.species, x, y):
            ax.text(xi, yi, sp, fontsize=8)
        if len(x) >= 2:
            slope, intercept, r, p, _ = stats.linregress(x, y)
            xs = np.linspace(x.min(), x.max(), 20)
            ax.plot(xs, intercept + slope * xs, color="black", lw=1)
            ax.set_title(f"{title}\nb={slope:.2f}, p={p:.3f}")
        ax.set_xlabel("log10 mass (kg)")
    fig.suptitle("CAGE promoter usage for metabolic genes vs body mass (4 spp)", fontsize=11)
    fig.savefig(FIG / "fig_A6_cage_promoter_vs_mass.png", dpi=160)
    plt.close(fig)

    # Chicken-only: alternative promoter switching for OXPHOS — compare dominant promoter identity
    # across fibroblast samples vs using atlas ESC if peak-level available.
    # Here: within GSE chicken CAGE, report genes with multi-promoter usage (shannon > 1).
    ch = big[big.species == "Chicken"]
    multi = ch[ch.shannon > 1.0].sort_values("shannon", ascending=False)
    multi.head(100).to_csv(RES / "cage_chicken_multipromoter_metabolic_top100.tsv", sep="\t", index=False)

    # Cross-species: genes where dominant_frac is low (promoter switching potential) and conserved
    wide = big.pivot_table(index="gene", columns="species", values="dominant_frac")
    wide["mean_dom"] = wide.mean(axis=1)
    wide["n_spp"] = wide.drop(columns=["mean_dom"]).notna().sum(axis=1)
    wide = wide[wide.n_spp >= 3].sort_values("mean_dom")
    wide.to_csv(RES / "cage_promoter_dominant_frac_wide.tsv", sep="\t")

    note = {
        "Q8": "CAGE promoter usage quantified for OXPHOS/TCA/glycolysis genes across 4 symbol-annotated species.",
        "metrics": "n_active_promoters, Shannon diversity, dominant promoter fraction, total promoter TPM",
        "esc_fib_peak_level": "Peak-level ESC vs fib promoter switching not fully resolved without shared peak IDs across libraries; gene-level atlas TPM used in A4/A5.",
    }
    with open(RES / "A6_cage_regulation_notes.json", "w") as f:
        json.dump(note, f, indent=2)
    print(summary.to_string(index=False))
    print("Done A6.")


if __name__ == "__main__":
    main()
