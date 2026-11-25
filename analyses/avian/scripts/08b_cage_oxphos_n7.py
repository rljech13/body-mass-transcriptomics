#!/usr/bin/env python3
"""CAGE OXPHOS/mito promoter metrics for all 7 fibroblast species.

Maps GEO CAGE peaks → gene symbols via Ruslan DE peak IDs
(`*_DE_epi(up)_fib(down)_promoters.xlsx` column `genes`), so Turkey/Peacock/Emu
are included (not a 4-species subset).
"""
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
RES = ROOT / "results" / "paperA_diff_pack"
FIG = ROOT / "figures"
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")

SPECIES = {
    "Chicken": {
        "cage": "GSE213253_Chicken_CAGE.xlsx",
        "de": "Chicken_DE_epi(up)_fib(down)_promoters.xlsx",
        "mass_kg": 0.9125,
        "sp_key": "chicken",
    },
    "Duck": {
        "cage": "GSE213253_Duck_CAGE.xlsx",
        "de": "Duck_DE_epi(up)_fib(down)_promoters.xlsx",
        "mass_kg": 1.121,
        "sp_key": "duck",
    },
    "Quail": {
        "cage": "GSE213253_Quail_CAGE.xlsx",
        "de": "Quail_DE_epi(up)_fib(down)_promoters.xlsx",
        "mass_kg": 0.0966,
        "sp_key": "quail",
    },
    "ZebraFinch": {
        "cage": "GSE213253_ZebraFinch_CAGE.xlsx",
        "de": "Zebrafinch_DE_epi(up)_fib(down)_promoters.xlsx",
        "mass_kg": 0.01193,
        "sp_key": "zebra",
    },
    "Turkey": {
        "cage": "GSE213253_Turkey_CAGE.xlsx",
        "de": "Turkey_DE_epi(up)_fib(down)_promoters.xlsx",
        "mass_kg": 5.811,
        "sp_key": "turkey",
    },
    "Peacock": {
        "cage": "GSE213253_Peacock_CAGE.xlsx",
        "de": "Peacock_DE_epi(up)_fib(down)_promoters.xlsx",
        "mass_kg": 4.09375,
        "sp_key": "peacock",
    },
    "Emu": {
        "cage": "GSE213253_Emu_CAGE.xlsx",
        "de": "Emu_DE_epi(up)_fib(down)_promoters.xlsx",
        "mass_kg": 36.2,
        "sp_key": "emu",
    },
}


def load_modules() -> dict[str, set[str]]:
    out = {}
    for p in (DATA / "genesets").glob("*.txt"):
        out[p.stem] = {ln.strip().upper() for ln in p.read_text().splitlines() if ln.strip()}
    return out


def peak_to_symbol(de_path: Path) -> pd.Series:
    de = pd.read_excel(de_path, usecols=["id", "genes"])
    de["genes"] = de["genes"].astype(str).str.upper()
    de = de[de["genes"].notna() & (de["genes"] != "NAN")]
    return de.drop_duplicates("id").set_index("id")["genes"]


def promoter_stats_for_genes(cage: pd.DataFrame, genes: set[str]) -> pd.DataFrame:
    norms = [c for c in cage.columns if str(c).startswith("norm.")]
    sub = cage[cage["symbol"].isin(genes)].copy()
    if sub.empty or not norms:
        return pd.DataFrame()
    sub["_expr"] = sub[norms].apply(pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
    rows = []
    for g, gdf in sub.groupby("symbol"):
        vals = gdf["_expr"].astype(float).values
        pos = vals[vals > 0]
        if len(pos) == 0:
            continue
        p = pos / pos.sum()
        rows.append(
            {
                "gene": g,
                "n_promoters": int(len(gdf)),
                "n_active_promoters": int((gdf["_expr"] > 0).sum()),
                "total_tpm": float(gdf["_expr"].sum()),
                "dominant_frac": float(p.max()),
                "shannon": float(-(p * np.log2(p + 1e-12)).sum()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    modules = load_modules()
    ox = modules["OXPHOS"]
    oxmito = ox | modules["MITO_BIOGENESIS"] | modules["TCA"]

    coverage_rows = []
    all_gene_rows = []

    for sp, meta in SPECIES.items():
        print(f"=== {sp} ===", flush=True)
        peak2sym = peak_to_symbol(RUSLAN / meta["de"])
        cage = pd.read_excel(DATA / meta["cage"])
        cage["symbol"] = cage["CAGE_cluster_id"].map(peak2sym)
        # fallback: native geneId if already a symbol
        if "geneId" in cage.columns:
            native = cage["geneId"].astype(str).str.upper()
            cage.loc[cage["symbol"].isna() & native.isin(oxmito), "symbol"] = native[
                cage["symbol"].isna() & native.isin(oxmito)
            ]

        frac = float(cage["symbol"].notna().mean())
        ox_n = int(cage.loc[cage["symbol"].isin(ox), "symbol"].nunique())
        oxmito_n = int(cage.loc[cage["symbol"].isin(oxmito), "symbol"].nunique())
        coverage_rows.append(
            {
                "species": sp,
                "sp_key": meta["sp_key"],
                "mass_kg": meta["mass_kg"],
                "n_peaks": len(cage),
                "frac_peaks_with_symbol": frac,
                "n_OXPHOS_genes": ox_n,
                "n_oxmito_genes": oxmito_n,
            }
        )
        print(f"  mapped {frac:.1%} peaks; OXPHOS genes={ox_n}; oxmito={oxmito_n}", flush=True)

        for label, gene_set in [("OXPHOS", ox), ("oxmito", oxmito)]:
            st = promoter_stats_for_genes(cage, gene_set)
            if st.empty:
                continue
            st["species"] = sp
            st["sp_key"] = meta["sp_key"]
            st["mass_kg"] = meta["mass_kg"]
            st["module_set"] = label
            all_gene_rows.append(st)

    cov = pd.DataFrame(coverage_rows)
    cov.to_csv(RES / "arm4_cage_n7_coverage.tsv", sep="\t", index=False)
    assert len(cov) == 7, cov
    assert (cov["n_OXPHOS_genes"] >= 10).all(), cov

    big = pd.concat(all_gene_rows, ignore_index=True)
    big.to_csv(RES / "arm4_cage_n7_promoter_genes.tsv", sep="\t", index=False)

    fit_all = []
    for label in ["OXPHOS", "oxmito"]:
        sub = big[big["module_set"] == label]
        summ = (
            sub.groupby("species")
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
        summ.to_csv(RES / f"arm4_cage_n7_{label}_species_summary.tsv", sep="\t", index=False)
        print(f"\n{label} summary n={len(summ)}\n", summ.to_string(index=False), flush=True)
        for col in ["mean_n_active", "mean_shannon", "mean_dominant_frac", "mean_total_tpm"]:
            sl, icept, r, p, se = stats.linregress(np.log10(summ["mass_kg"]), summ[col])
            fit_all.append(
                {
                    "module_set": label,
                    "metric": col,
                    "n_species": len(summ),
                    "slope": sl,
                    "r": r,
                    "p": p,
                    "se": se,
                }
            )

    fits = pd.DataFrame(fit_all)
    fits.to_csv(RES / "arm4_cage_n7_vs_mass_ols.tsv", sep="\t", index=False)
    print("\nOLS\n", fits.to_string(index=False), flush=True)

    # figure OXPHOS n=7
    summ = pd.read_csv(RES / "arm4_cage_n7_OXPHOS_species_summary.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), constrained_layout=True)
    for ax, col, title in zip(
        axes,
        ["mean_n_active", "mean_shannon", "mean_total_tpm"],
        ["Active promoters", "Shannon", "Total TPM"],
    ):
        x = np.log10(summ["mass_kg"])
        y = summ[col]
        ax.scatter(x, y, s=50)
        for sp, xi, yi in zip(summ["species"], x, y):
            ax.text(xi, yi, sp, fontsize=7)
        sl, icept, _, p, _ = stats.linregress(x, y)
        xs = np.linspace(x.min(), x.max(), 40)
        ax.plot(xs, icept + sl * xs, color="C1")
        ax.set_title(f"{title}\nslope={sl:.3f}, p={p:.3f}")
        ax.set_xlabel("log10 mass (kg)")
    fig.suptitle("CAGE OXPHOS promoters vs mass (n=7 fibroblasts)")
    fig.savefig(FIG / "fig_arm4_cage_oxphos_promoters_vs_mass_n7.png", dpi=160)
    # replace old n=4 figure name used in reports
    fig.savefig(FIG / "fig_arm4_cage_oxphos_promoters_vs_mass.png", dpi=160)
    plt.close(fig)

    summary = {
        "n_species": 7,
        "mapping": "CAGE_cluster_id → Ruslan DE id → genes symbol",
        "coverage": cov.to_dict(orient="records"),
        "ols": fits.to_dict(orient="records"),
        "note": "Supersedes prior n=4 symbol-only CAGE OXPHOS fits.",
    }
    (RES / "arm4_n7_summary.json").write_text(json.dumps(summary, indent=2))
    # overwrite arm4_summary key bits
    old = RES / "arm4_summary.json"
    prev = json.loads(old.read_text()) if old.exists() else {}
    prev.update(
        {
            "n_species_usable_fib_cage": 7,
            "n4_deprecated": True,
            "n7": summary,
        }
    )
    old.write_text(json.dumps(prev, indent=2, default=str))
    print("DONE n=7", RES)


if __name__ == "__main__":
    main()
