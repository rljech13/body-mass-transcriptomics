#!/usr/bin/env python3
"""Paper A — alternative scaling (NOT Kleiber).

Ask: given the n=7 CAGE OXPHOS signal, what phenomenological formula
describes promoter architecture better than forced M^{-1/4}?

Models compared (AIC):
  - free power: log Y ~ α log M
  - log-linear: Y ~ a + b log10 M
  - incubation: Y ~ a + b incubation
  - longevity / genome size
  - promoter simplicity index S = 1 - H/Hmax
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
RES = ROOT / "results" / "paperA_diff_pack"
OUT = RES / "arm5_alt_scaling"
FIG = ROOT / "figures"

TRAITS = pd.DataFrame(
    [
        ("chicken", "Chicken", 912.5, 20.5, 30.0, 1.223),
        ("duck", "Duck", 1121.0, 26.5, 29.1, 1.428),
        ("emu", "Emu", 36200.0, 56.0, 16.6, 1.594),
        ("peacock", "Peacock", 4093.75, 28.0, 25.0, 1.130),
        ("quail", "Quail", 96.6, 17.0, 14.58, 1.379),
        ("turkey", "Turkey", 5811.0, 28.0, 13.0, 1.281),
        ("zebra", "ZebraFinch", 11.93, 14.15, 14.5, 1.223),
    ],
    columns=["sp_key", "species", "mass_g", "incubation_d", "max_longevity_y", "genome_size_gb"],
)


def aic_ols(y: np.ndarray, yhat: np.ndarray, k: int) -> float:
    n = len(y)
    rss = float(np.sum((y - yhat) ** 2))
    if rss <= 0:
        return -np.inf
    return n * np.log(rss / n) + 2 * k


def fit_models(y: np.ndarray, df: pd.DataFrame) -> list[dict]:
    rows = []
    logM = np.log10(df["mass_g"].values)
    logM_nat = np.log(df["mass_g"].values)

    # 1) Y ~ a + b log10 M  (log-linear)
    sl, icept, r, p, se = stats.linregress(logM, y)
    yhat = icept + sl * logM
    rows.append(
        {
            "model": "Y = a + b·log10(M)",
            "formula_family": "log-linear",
            "param": "b",
            "param_value": sl,
            "intercept": icept,
            "r": r,
            "p": p,
            "AIC": aic_ols(y, yhat, 2),
            "note": "phenomenological; b free (not Kleiber)",
        }
    )

    # 2) free power-law: log10 Y ~ α log10 M  => Y ∝ M^α
    if np.all(y > 0):
        a, b, r2, p2, se2 = stats.linregress(logM, np.log10(y))
        yhat2 = 10 ** (b + a * logM)
        rows.append(
            {
                "model": "Y ∝ M^α  (free α)",
                "formula_family": "power-law",
                "param": "α",
                "param_value": a,
                "intercept": b,
                "r": r2,
                "p": p2,
                "AIC": aic_ols(y, yhat2, 2),
                "note": f"α free; Kleiber cellular would be −0.25 (not imposed)",
            }
        )

    # 3) incubation
    x = df["incubation_d"].values.astype(float)
    sl, icept, r, p, se = stats.linregress(x, y)
    yhat = icept + sl * x
    rows.append(
        {
            "model": "Y = a + b·incubation_d",
            "formula_family": "life-history",
            "param": "b",
            "param_value": sl,
            "intercept": icept,
            "r": r,
            "p": p,
            "AIC": aic_ols(y, yhat, 2),
            "note": "development time, not metabolic scaling",
        }
    )

    # 4) longevity
    x = df["max_longevity_y"].values.astype(float)
    sl, icept, r, p, se = stats.linregress(x, y)
    yhat = icept + sl * x
    rows.append(
        {
            "model": "Y = a + b·longevity_y",
            "formula_family": "life-history",
            "param": "b",
            "param_value": sl,
            "intercept": icept,
            "r": r,
            "p": p,
            "AIC": aic_ols(y, yhat, 2),
            "note": "",
        }
    )

    # 5) genome size
    x = df["genome_size_gb"].values.astype(float)
    sl, icept, r, p, se = stats.linregress(x, y)
    yhat = icept + sl * x
    rows.append(
        {
            "model": "Y = a + b·genome_Gb",
            "formula_family": "genome",
            "param": "b",
            "param_value": sl,
            "intercept": icept,
            "r": r,
            "p": p,
            "AIC": aic_ols(y, yhat, 2),
            "note": "",
        }
    )

    # 6) incubation + log mass (2 predictors) — partial
    X = np.column_stack([np.ones(len(y)), logM, df["incubation_d"].values])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    # approximate p via F vs null
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    rows.append(
        {
            "model": "Y = a + b·log10(M) + c·incubation",
            "formula_family": "combined",
            "param": "b,c",
            "param_value": float(beta[1]),
            "intercept": float(beta[0]),
            "r": float(np.sqrt(max(r2, 0)) * np.sign(beta[1])),
            "p": np.nan,
            "AIC": aic_ols(y, yhat, 3),
            "note": f"c={beta[2]:.5g}; exploratory 2-predictor",
        }
    )
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    cage = pd.read_csv(RES / "arm4_cage_n7_OXPHOS_species_summary.tsv", sep="\t")
    df = cage.merge(TRAITS, on="species")
    # simplicity: for mean shannon across genes with typically 1–few promoters,
    # Hmax ≈ log2(n_active); use S = 1 - shannon / log2(mean_n_active.clip(2))
    hmax = np.log2(np.maximum(df["mean_n_active"].values, 2.0))
    df["simplicity_S"] = 1.0 - df["mean_shannon"].values / hmax
    df["log10_mass_g"] = np.log10(df["mass_g"])
    df.to_csv(OUT / "endpoints_with_traits_n7.tsv", sep="\t", index=False)

    endpoints = {
        "OXPHOS_mean_shannon": df["mean_shannon"].values,
        "OXPHOS_mean_total_tpm": df["mean_total_tpm"].values,
        "OXPHOS_mean_dominant_frac": df["mean_dominant_frac"].values,
        "OXPHOS_simplicity_S": df["simplicity_S"].values,
        "OXPHOS_mean_n_active": df["mean_n_active"].values,
    }

    # also pull gemcat mito / oxphos module logFC
    gem = pd.read_csv(
        ROOT / "results/gemcat_multispp_n7/multispp_mito_median_vs_mass.tsv", sep="\t"
    )
    oxm = pd.read_csv(RES / "arm3_oxphos_module_epi_fib.tsv", sep="\t")
    gdf = gem.drop(columns=["mass_g"], errors="ignore").merge(
        TRAITS, left_on="species", right_on="sp_key"
    )
    odf = oxm.drop(columns=["mass_g"], errors="ignore").merge(
        TRAITS, left_on="species", right_on="sp_key"
    )
    endpoints_extra = {
        "GEMCAT_mito_median_log2": (gdf, "median_log2_mito_ies1300"),
        "OXPHOS_module_logFC_epi_fib": (odf, "logFC_epi_minus_fib"),
    }

    all_fits = []
    for name, y in endpoints.items():
        for row in fit_models(y, df):
            row["endpoint"] = name
            all_fits.append(row)

    for name, (tab, col) in endpoints_extra.items():
        y = tab[col].values.astype(float)
        # need trait columns aligned — rebuild mini df
        mini = tab.copy()
        for row in fit_models(y, mini):
            row["endpoint"] = name
            all_fits.append(row)

    fits = pd.DataFrame(all_fits)
    fits.to_csv(OUT / "alt_scaling_model_comparison.tsv", sep="\t", index=False)

    # best model per endpoint by AIC
    best = (
        fits.sort_values("AIC")
        .groupby("endpoint", as_index=False)
        .first()
        [["endpoint", "model", "formula_family", "param", "param_value", "r", "p", "AIC", "note"]]
    )
    best.to_csv(OUT / "alt_scaling_best_by_AIC.tsv", sep="\t", index=False)

    # free α for Shannon and TPM specifically
    alpha_rows = []
    for col, label in [
        ("mean_shannon", "Shannon"),
        ("mean_total_tpm", "total_TPM"),
        ("mean_dominant_frac", "dominant_frac"),
        ("simplicity_S", "simplicity_S"),
    ]:
        y = df[col].values.astype(float)
        logM = np.log10(df["mass_g"])
        if np.all(y > 0):
            alpha, icept, r, p, se = stats.linregress(logM, np.log10(y))
            alpha_rows.append(
                {
                    "endpoint": label,
                    "alpha_free": alpha,
                    "r": r,
                    "p": p,
                    "se": se,
                    "kleiber_cellular_ref": -0.25,
                    "delta_vs_m14": alpha - (-0.25),
                }
            )
    alphas = pd.DataFrame(alpha_rows)
    alphas.to_csv(OUT / "free_power_alpha.tsv", sep="\t", index=False)

    # gene-level: which OXPHOS genes have shannon/n_active correlating with mass
    genes = pd.read_csv(RES / "arm4_cage_n7_promoter_genes.tsv", sep="\t")
    genes = genes[genes["module_set"] == "OXPHOS"].copy()
    genes = genes.merge(TRAITS[["species", "mass_g"]], on="species")
    gene_fits = []
    for g, sub in genes.groupby("gene"):
        if len(sub) < 5:
            continue
        x = np.log10(sub["mass_g"])
        for metric in ["shannon", "n_active_promoters", "total_tpm", "dominant_frac"]:
            y = sub[metric].values.astype(float)
            if np.std(y) < 1e-12:
                continue
            sl, icept, r, p, se = stats.linregress(x, y)
            gene_fits.append(
                {"gene": g, "metric": metric, "n_spp": len(sub), "slope": sl, "r": r, "p": p}
            )
    gf = pd.DataFrame(gene_fits)
    gf.to_csv(OUT / "gene_level_oxphos_cage_vs_mass.tsv", sep="\t", index=False)
    # top drivers of Shannon decline
    shan = gf[gf.metric == "shannon"].sort_values("p")
    shan.head(20).to_csv(OUT / "top_genes_shannon_vs_mass.tsv", sep="\t", index=False)

    # figure: Shannon formulas
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.7), constrained_layout=True)
    y = df["mean_shannon"].values
    logM = np.log10(df["mass_g"].values)

    # panel 1: log-linear
    ax = axes[0]
    ax.scatter(logM, y, s=50)
    sl, icept, r, p, _ = stats.linregress(logM, y)
    xs = np.linspace(logM.min(), logM.max(), 50)
    ax.plot(xs, icept + sl * xs, color="C1")
    for _, row in df.iterrows():
        ax.text(np.log10(row.mass_g), row.mean_shannon, row.species, fontsize=7)
    ax.set_title(f"H = {icept:.3f} + ({sl:.4f})·log10(M)\nr={r:.2f}, p={p:.3g}")
    ax.set_xlabel("log10 mass (g)")
    ax.set_ylabel("mean OXPHOS Shannon H")

    # panel 2: free power
    ax = axes[1]
    ax.scatter(logM, np.log10(y), s=50)
    a, b, r, p, _ = stats.linregress(logM, np.log10(y))
    ax.plot(xs, b + a * xs, color="C1")
    ax.axhline(np.nan, color="none")
    # reference slope -1/4 on log-log of H (not meaningful physiologically but shown)
    ax.plot(xs, (b + a * logM.mean()) - 0.25 * (xs - logM.mean()), "k--", lw=0.8, label="α=−1/4 ref")
    ax.legend(fontsize=7)
    ax.set_title(f"H ∝ M^{a:.3f} (free α)\nr={r:.2f}, p={p:.3g}")
    ax.set_xlabel("log10 mass (g)")
    ax.set_ylabel("log10 Shannon")

    # panel 3: incubation
    ax = axes[2]
    x = df["incubation_d"].values
    ax.scatter(x, y, s=50)
    sl, icept, r, p, _ = stats.linregress(x, y)
    xs2 = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs2, icept + sl * xs2, color="C1")
    for _, row in df.iterrows():
        ax.text(row.incubation_d, row.mean_shannon, row.species, fontsize=7)
    ax.set_title(f"H vs incubation\nslope={sl:.4f}, p={p:.3g}")
    ax.set_xlabel("incubation (days)")
    ax.set_ylabel("mean OXPHOS Shannon H")

    fig.suptitle("Alternative (non-Kleiber) descriptions of OXPHOS promoter diversity")
    fig.savefig(FIG / "fig_arm5_alt_scaling_shannon.png", dpi=160)
    plt.close(fig)

    # pick proposed formula: best AIC among Shannon models that are simple
    shan_fits = fits[fits.endpoint == "OXPHOS_mean_shannon"].sort_values("AIC")
    proposed = shan_fits.iloc[0].to_dict()

    # explicit recommended formula from log-linear (most interpretable for H near 0)
    sl, icept, r, p, se = stats.linregress(logM, y)
    recommendation = {
        "proposed_formula": f"H_OXPHOS ≈ {icept:.4f} + ({sl:.4f}) · log10(M_g)",
        "or_power": f"H_OXPHOS ∝ M^{alphas.loc[alphas.endpoint=='Shannon','alpha_free'].values[0]:.3f}",
        "interpretation": (
            "Promoter diversity of OXPHOS genes declines approximately log-linearly "
            "with adult body mass across 7 cultured fibroblast species. "
            "This is a promoter-architecture scaling, not organismal BMR Kleiber."
        ),
        "best_AIC_model": proposed,
        "free_alpha_shannon": float(alphas.loc[alphas.endpoint == "Shannon", "alpha_free"].iloc[0]),
        "kleiber_not_used": True,
        "n": 7,
    }
    (OUT / "proposed_formula.json").write_text(json.dumps(recommendation, indent=2, default=str))

    print("BEST BY AIC:")
    print(best.to_string(index=False))
    print("\nFREE ALPHA:")
    print(alphas.to_string(index=False))
    print("\nPROPOSED:")
    print(json.dumps(recommendation, indent=2))
    print("\nTop Shannon-vs-mass genes:")
    print(shan.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
