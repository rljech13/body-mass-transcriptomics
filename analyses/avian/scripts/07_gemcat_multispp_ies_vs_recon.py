#!/usr/bin/env python3
"""Same design as Ruslan (per-species epi/fib FC → GEMCAT → mito summary ~ mass),
but primary network = iES1300; Recon3D recomputed as baseline.

Why we previously had n=1 for GEMCAT: A5 used chicken CAGE atlas ESC vs fib only.
Multi-spp epi/fib FC exists in Ruslan bird_sizes; this script uses that design.
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
DATA = ROOT / "data"
RES = ROOT / "results" / "gemcat_multispp_n7"
FIG = ROOT / "figures"
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")
IES = DATA / "ies1300/S1_File/S1_File in sbml format.xml"
RECON = DATA / "Recon3D.xml"
ENV = {**os.environ, "LD_PRELOAD": "/usr/lib/x86_64-linux-gnu/libsqlite3.so.0"}

SPECIES = [
    ("chicken", "Chicken"),
    ("duck", "Duck"),
    ("emu", "Emu"),
    ("peacock", "Peacock"),
    ("quail", "Quail"),
    ("turkey", "Turkey"),
    ("zebra", "Zebrafinch"),
]

# adult body mass in g — same table as Ruslan sizes.Rmd
MASS_G = {
    "chicken": 912.5,
    "duck": 1121.0,
    "emu": 36200.0,
    "peacock": 4093.75,
    "quail": 96.6,
    "turkey": 5811.0,
    "zebra": 11.93,
}


def run_gemcat(model: Path, de_csv: Path, out_csv: Path, effect_col: str = "FC") -> None:
    cmd = ["gemcat", str(model), str(de_csv), "-e", effect_col, "-o", str(out_csv)]
    print("RUN:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=ENV)


def sort_gemcat(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, header=None)
    df = df.iloc[:, :2]
    df.columns = ["metabolite", "score"]
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    return df.dropna(subset=["score"])


def build_ies_fc_from_xlsx(sp_file: str, sym2chk: pd.DataFrame, ies_genes: set[str]) -> pd.DataFrame:
    xlsx = RUSLAN / f"{sp_file}_DE_epi(up)_fib(down)_promoters.xlsx"
    df = pd.read_excel(xlsx, usecols=["genes", "logFC"])
    df = df.dropna(subset=["genes"]).copy()
    df["genes"] = df["genes"].astype(str).str.upper()
    df = df.sort_values("logFC", key=lambda s: s.abs(), ascending=False).drop_duplicates("genes")
    m = df.merge(sym2chk, left_on="genes", right_on="symbol", how="inner")
    m = m[m["chicken_entrez"].isin(ies_genes)].copy()
    m["FC"] = np.power(2.0, m["logFC"].astype(float))
    # clip extreme FC that can explode network weights
    m["FC"] = m["FC"].clip(lower=1e-6, upper=1e6)
    out = m[["chicken_entrez", "FC"]].rename(columns={"chicken_entrez": "gene"})
    return out.drop_duplicates("gene")


def expand_ruslan_fixed_to_path(sp_file: str) -> Path:
    """Reuse Ruslan's already-mapped Recon3D AT FC files (same design)."""
    return RUSLAN / f"{sp_file}_FC_linear_fixedID.csv"


def mito_median_log2(scores: pd.DataFrame, network: str) -> float:
    """Median log2(score) over mitochondrial metabolites with score>0.

    Matching Ruslan: Recon mitochondria = ids ending `_m` (compartment_id m).
    iES1300: cobra compartment `m` / ids ending `[m]`.
    """
    s = scores.copy()
    if "metabolite" not in s.columns and "metabolite_id" in s.columns:
        s = s.rename(columns={"metabolite_id": "metabolite"})
    s["score"] = pd.to_numeric(s["score"], errors="coerce")
    s = s[np.isfinite(s["score"]) & (s["score"] > 0)]
    mets = s["metabolite"].astype(str)
    if network == "ies":
        mask = mets.str.endswith("[m]")
    elif network == "recon":
        mask = mets.str.endswith("_m")
    else:
        raise ValueError(network)
    sub = s.loc[mask]
    if len(sub) < 10:
        raise RuntimeError(f"Too few mito mets for {network}: {len(sub)}")
    return float(np.median(np.log2(sub["score"].to_numpy())))


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    import cobra

    ies_model = cobra.io.read_sbml_model(str(IES))
    ies_genes = {g.id for g in ies_model.genes}

    sym2chk = pd.read_csv(DATA / "symbol_to_chicken_entrez.tsv", sep="\t")
    sym2chk["symbol"] = sym2chk["symbol"].astype(str).str.upper()
    sym2chk["chicken_entrez"] = sym2chk["chicken_entrez"].astype(int).astype(str)

    rows = []
    for sp, sp_file in SPECIES:
        print(f"\n=== {sp} ===", flush=True)
        # --- iES1300 primary ---
        ies_fc = build_ies_fc_from_xlsx(sp_file, sym2chk, ies_genes)
        ies_in = RES / f"{sp}_FC_ies1300_chicken_entrez.csv"
        ies_fc.to_csv(ies_in, index=False)
        ies_out = RES / f"{sp}_gemcat_ies1300.csv"
        if not ies_out.exists() or ies_out.stat().st_size < 10:
            run_gemcat(IES, ies_in, ies_out, effect_col="FC")
        ies_scores = sort_gemcat(ies_out)
        ies_scores.to_csv(RES / f"{sp}_gemcat_ies1300_sorted.tsv", sep="\t", index=False)
        ies_med = mito_median_log2(ies_scores, "ies")

        # --- Recon3D baseline: Ruslan linear FC files (same gene weights as him) ---
        recon_in = expand_ruslan_fixed_to_path(sp_file)
        recon_out = RES / f"{sp}_gemcat_recon3d.csv"
        if not recon_out.exists() or recon_out.stat().st_size < 10:
            run_gemcat(RECON, recon_in, recon_out, effect_col="FC")
        recon_scores = sort_gemcat(recon_out)
        recon_scores.to_csv(RES / f"{sp}_gemcat_recon3d_sorted.tsv", sep="\t", index=False)
        recon_med = mito_median_log2(recon_scores, "recon")

        # Ruslan's stored annotated result (exact mito compartment)
        ann = RUSLAN / f"result_annotated_{sp_file}.csv"
        if ann.exists():
            rus = pd.read_csv(ann)
            rus = rus[rus["compartment_name"] == "mitochondria"][
                ["metabolite_id", "score"]
            ].rename(columns={"metabolite_id": "metabolite"})
            rus_med = mito_median_log2(rus, "recon")
        else:
            rus = sort_gemcat(RUSLAN / f"result_{sp_file}.csv")
            rus_med = mito_median_log2(rus, "recon")

        rows.append(
            {
                "species": sp,
                "mass_g": MASS_G[sp],
                "log10_mass_g": np.log10(MASS_G[sp]),
                "n_genes_ies": len(ies_fc),
                "median_log2_mito_ies1300": ies_med,
                "median_log2_mito_recon3d_recompute": recon_med,
                "median_log2_mito_ruslan_result": rus_med,
            }
        )
        print(
            f"  iES n_genes={len(ies_fc)} mito_med={ies_med:.4f} | "
            f"Recon recompute={recon_med:.4f} | Ruslan stored={rus_med:.4f}",
            flush=True,
        )

    tab = pd.DataFrame(rows)
    tab.to_csv(RES / "multispp_mito_median_vs_mass.tsv", sep="\t", index=False)

    # OLS fits
    fit_rows = []
    for col, label in [
        ("median_log2_mito_ies1300", "iES1300_primary"),
        ("median_log2_mito_recon3d_recompute", "Recon3D_recompute"),
        ("median_log2_mito_ruslan_result", "Ruslan_stored_Recon3D"),
    ]:
        x = tab["log10_mass_g"].values
        y = tab[col].values
        slope, intercept, r, p, se = stats.linregress(x, y)
        fit_rows.append(
            {
                "model": label,
                "n": len(tab),
                "slope": slope,
                "intercept": intercept,
                "r": r,
                "p": p,
                "se": se,
            }
        )
        print(f"FIT {label}: slope={slope:.4f} p={p:.4g} r={r:.3f}")
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(RES / "multispp_mito_median_ols_fits.tsv", sep="\t", index=False)

    # figure
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
    for ax, col, title in [
        (axes[0], "median_log2_mito_ies1300", "iES1300 PRIMARY"),
        (axes[1], "median_log2_mito_recon3d_recompute", "Recon3D recompute"),
        (axes[2], "median_log2_mito_ruslan_result", "Ruslan stored Recon3D"),
    ]:
        ax.scatter(tab["log10_mass_g"], tab[col], s=50)
        for _, r in tab.iterrows():
            ax.text(r["log10_mass_g"], r[col], r["species"], fontsize=7, ha="left")
        slope, intercept, _, p, _ = stats.linregress(tab["log10_mass_g"], tab[col])
        xs = np.linspace(tab["log10_mass_g"].min(), tab["log10_mass_g"].max(), 50)
        ax.plot(xs, intercept + slope * xs, color="C1", lw=1.5)
        ax.set_title(f"{title}\nslope={slope:.3f}, p={p:.3f}")
        ax.set_xlabel("log10 adult body mass (g)")
        ax.set_ylabel("median log2 GEMCAT score (mito)")
    fig.suptitle("Same design (epi/fib FC → GEMCAT) — network comparison, n=7")
    fig.savefig(FIG / "fig_gemcat_multispp_ies_vs_recon_vs_mass.png", dpi=160)
    plt.close(fig)

    # short note
    note = RES / "DESIGN_NOTE.md"
    note.write_text(
        """# Multi-spp GEMCAT design note

## Why GEMCAT was chicken-only before
Paper-A step A5 used the **chicken CAGE atlas** ESC vs fibroblast contrast only.
Matched multi-spp ESC/epiblast was not in that atlas pipeline.

## This run
Same design as Ruslan `bird_sizes`:
per-species promoter DE **epi vs fib** → fold-change weights → GEMCAT →
median log2 score over mitochondrial metabolites → OLS vs log10(mass).

| Arm | Network | Gene weights |
|-----|---------|--------------|
| Primary | **iES1300** (chicken) | symbols from Ruslan xlsx → chicken Entrez ∩ iES genes |
| Baseline | Recon3D (human) | Ruslan `*_FC_linear_fixedID.csv` (recomputed gemcat) |
| Reference | Ruslan stored `result_*.csv` | his prior Recon3D outputs |

Primary question: does switching network from Recon3D to iES1300 change the Kleiber (mass) conclusion under the **same** multi-spp design?
"""
    )
    print("Wrote", RES)
    print(tab.to_string(index=False))
    print(fits.to_string(index=False))


if __name__ == "__main__":
    main()
