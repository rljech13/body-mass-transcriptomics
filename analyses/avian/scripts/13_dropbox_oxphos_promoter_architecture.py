#!/usr/bin/env python3
"""OXPHOS promoter architecture & regulation from Dropbox CAGE fib+epi tables (n=7).

Maps Augustus / mixed geneIds → chicken symbols via:
  1) chicken_genes column (Dropbox)
  2) native geneId when already a symbol
  3) Ruslan DE peak id → genes

Endpoints (fib and epi separately; plus log2 epi/fib where both exist):
  n_active, Shannon H, dominant_frac, simplicity S, total TPM

Controls (honest):
  - same-strand filter (major strand by TPM)
  - TPM threshold sensitivity
  - do NOT claim architecture simplification if TPM~mass dominates

Language: "на этих данных сдвига не видно" — not "сдвига нет".
"""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
DATA = ROOT / "data"
EXP = Path("/path/to/storage/kleiber_birds_gemcat/dropbox_20260808/unpacked/exp_tab")
# local mirror if present
EXP_LOCAL = DATA / "dropbox_20260808"
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")
RES = ROOT / "results" / "oxphos_promoter_architecture_n7"
FIG = ROOT / "figures"

SPECIES = {
    "Zebra": {"file": "Zebra_fib_epi_exp.txt", "de": "Zebrafinch_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 11.93},
    "Quail": {"file": "Quail_fib_epi_exp.txt", "de": "Quail_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 96.6},
    "Chicken": {"file": "Chicken_fib_epi_exp.txt", "de": "Chicken_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 912.5},
    "Duck": {"file": "Duck_fib_epi_exp.txt", "de": "Duck_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 1121.0},
    "Peacock": {"file": "Peacock_fib_epi_exp.txt", "de": "Peacock_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 4093.75},
    "Turkey": {"file": "Turkey_fib_epi_exp.txt", "de": "Turkey_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 5811.0},
    "Emu": {"file": "Emu_fib_epi_exp.txt", "de": "Emu_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 36200.0},
}

# HGNC old → preferred modern (when modern is in OXPHOS list)
PREFERRED_ALIAS = {
    "ATP5A1": "ATP5F1A",
    "ATP5B": "ATP5F1B",
    "ATP5C1": "ATP5F1C",
    "ATP5D": "ATP5F1D",
    "ATP5E": "ATP5F1E",
    "ATP5F1": "ATP5PB",
    "ATP5G1": "ATP5MC1",
    "ATP5G3": "ATP5MC3",
    "ATP5H": "ATP5PD",
    "ATP5J": "ATP5PF",
    "ATP5J2": "ATP5MF",
    "ATP5L": "ATP5MG",
    "ATP5O": "ATP5PO",
}


def load_oxphos() -> set[str]:
    return {ln.strip().upper() for ln in (DATA / "genesets/OXPHOS.txt").read_text().splitlines() if ln.strip()}


def is_augustus_id(s: str) -> bool:
    u = str(s).upper().strip()
    return bool(re.match(r"^G\d+(\.T\d+)?$", u)) or u in {"NAN", "NA", "", "NONE", "."}


def normalize_symbol(sym: str) -> str:
    """Uppercase + HGNC preferred alias; keep non-OXPHOS symbols too."""
    s = str(sym).upper().strip()
    if not s or s in {"NAN", "NA", "NONE", "."}:
        return ""
    if is_augustus_id(s):
        return ""
    return PREFERRED_ALIAS.get(s, s)


def exp_path(fname: str) -> Path:
    p = EXP / fname
    if p.exists():
        return p
    p2 = EXP_LOCAL / fname
    if p2.exists():
        return p2
    raise FileNotFoundError(fname)


def peak_to_symbol_de(de_path: Path) -> pd.Series:
    de = pd.read_excel(de_path, usecols=["id", "genes", "geneId"])
    for col in ("genes", "geneId"):
        de[col] = de[col].astype(str).str.upper()
    de["sym"] = de["genes"].where(~de["genes"].isin(["NAN", "NA", ""]), de["geneId"])
    de["sym"] = de["sym"].map(normalize_symbol)
    de = de[de["sym"] != ""].drop_duplicates("id")
    return de.set_index("id")["sym"]


def annotate_table(df: pd.DataFrame, id2sym: pd.Series) -> pd.DataFrame:
    out = df.copy()
    gid = out["geneId"].astype(str)
    cg = out["chicken_genes"].astype(str) if "chicken_genes" in out.columns else pd.Series([""] * len(out))
    sym = []
    for g, c in zip(gid, cg):
        cu = str(c).upper().strip()
        gu = str(g).upper().strip()
        pick = ""
        if cu and not is_augustus_id(cu):
            pick = cu
        elif not is_augustus_id(gu):
            pick = gu
        sym.append(normalize_symbol(pick) if pick else "")
    out["symbol"] = sym
    from_de = out["id"].map(id2sym)
    fill = out["symbol"].eq("") | out["symbol"].isna()
    out.loc[fill, "symbol"] = from_de[fill].fillna("").map(normalize_symbol)
    out["symbol"] = out["symbol"].fillna("").astype(str)
    return out


def lib_cols(df: pd.DataFrame, tissue: str) -> list[str]:
    """tissue: 'fib' → norm.C*, 'epi' → norm.E*"""
    pref = "norm.C" if tissue == "fib" else "norm.E"
    return [c for c in df.columns if str(c).startswith(pref)]


def metrics_for_gene(gdf: pd.DataFrame, expr_col: str, tpm_thr: float, same_strand: bool) -> dict:
    g = gdf.copy()
    g["_e"] = pd.to_numeric(g[expr_col], errors="coerce").fillna(0.0)
    g = g[g["_e"] >= tpm_thr]
    if g.empty:
        return {}
    if same_strand:
        maj = g.groupby("strand")["_e"].sum().idxmax()
        g = g[g["strand"] == maj]
    if g.empty:
        return {}
    vals = g["_e"].astype(float).values
    p = vals / vals.sum()
    H = float(-(p * np.log2(p + 1e-12)).sum())
    n = int(len(g))
    return {
        "n_active": n,
        "shannon": H,
        "dominant_frac": float(p.max()),
        "total_tpm": float(vals.sum()),
        "simplicity_S": float(1.0 - H / np.log2(max(n, 2))),
    }


def fit_mass(df: pd.DataFrame, col: str) -> dict:
    x = np.log10(df["mass_g"].values.astype(float))
    y = df[col].values.astype(float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return {"n": int(m.sum()), "slope": np.nan, "r": np.nan, "p": np.nan, "se": np.nan}
    slope, icept, r, p, se = stats.linregress(x[m], y[m])
    return {"n": int(m.sum()), "slope": float(slope), "intercept": float(icept), "r": float(r), "p": float(p), "se": float(se)}


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    ox = load_oxphos()
    thresholds = [0.0, 1.0, 5.0, 10.0]
    tissues = ("fib", "epi")

    coverage_rows = []
    sens_rows = []
    per_gene_rows = []
    maps_rows = []

    for sp, meta in SPECIES.items():
        print(f"=== {sp} ===", flush=True)
        df = pd.read_csv(exp_path(meta["file"]), sep="\t")
        id2 = peak_to_symbol_de(RUSLAN / meta["de"])
        df = annotate_table(df, id2)
        # keep promoter-class peaks for architecture (class column)
        if "class" in df.columns:
            prom = df[df["class"].astype(str).str.lower().eq("promoter")].copy()
        else:
            prom = df.copy()
        ox_mask = prom["symbol"].isin(ox)
        n_genes = int(prom.loc[ox_mask, "symbol"].nunique())
        n_peaks = int(ox_mask.sum())
        coverage_rows.append(
            {
                "species": sp,
                "mass_g": meta["mass_g"],
                "n_oxphos_genes": n_genes,
                "n_oxphos_promoter_peaks": n_peaks,
                "frac_peaks_mapped": float((df["symbol"] != "").mean()),
                "n_fib_libs": len(lib_cols(df, "fib")),
                "n_epi_libs": len(lib_cols(df, "epi")),
                "example_genes": ",".join(sorted(prom.loc[ox_mask, "symbol"].unique())[:12]),
            }
        )
        # save per-species peak map for OXPHOS
        submap = prom.loc[ox_mask, ["id", "seqnames", "strand", "tss", "geneId", "symbol", "class"]].copy()
        submap["species"] = sp
        maps_rows.append(submap)
        print(f"  OXPHOS genes={n_genes} peaks={n_peaks}", flush=True)

        ox_tab = prom.loc[ox_mask].copy()
        if ox_tab.empty:
            continue
        ox_tab["strand"] = ox_tab["strand"].astype(str)

        for tissue in tissues:
            cols = lib_cols(ox_tab, tissue)
            if not cols:
                continue
            ox_tab[f"_expr_{tissue}"] = ox_tab[cols].apply(pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)

        for tissue in tissues:
            ecol = f"_expr_{tissue}"
            if ecol not in ox_tab.columns:
                continue
            for thr in thresholds:
                for same in (False, True):
                    gmets = []
                    for gene, gdf in ox_tab.groupby("symbol"):
                        m = metrics_for_gene(gdf, ecol, thr, same_strand=same)
                        if not m:
                            continue
                        gmets.append({"gene": gene, **m})
                        if thr == 1.0 and same:
                            per_gene_rows.append(
                                {
                                    "species": sp,
                                    "mass_g": meta["mass_g"],
                                    "tissue": tissue,
                                    "gene": gene,
                                    **m,
                                }
                            )
                    if not gmets:
                        continue
                    gm = pd.DataFrame(gmets)
                    means = gm[["n_active", "shannon", "dominant_frac", "total_tpm", "simplicity_S"]].mean()
                    sens_rows.append(
                        {
                            "species": sp,
                            "mass_g": meta["mass_g"],
                            "tissue": tissue,
                            "tpm_thr": thr,
                            "same_strand": same,
                            "n_genes": len(gm),
                            **means.to_dict(),
                        }
                    )

        # FC architecture: per-gene log2(epi/fib) of total TPM (thr=1, same-strand)
        if "_expr_fib" in ox_tab.columns and "_expr_epi" in ox_tab.columns:
            for gene, gdf in ox_tab.groupby("symbol"):
                mf = metrics_for_gene(gdf, "_expr_fib", 1.0, True)
                me = metrics_for_gene(gdf, "_expr_epi", 1.0, True)
                if not mf or not me:
                    continue
                per_gene_rows.append(
                    {
                        "species": sp,
                        "mass_g": meta["mass_g"],
                        "tissue": "log2_epi_over_fib",
                        "gene": gene,
                        "n_active": me["n_active"] - mf["n_active"],
                        "shannon": me["shannon"] - mf["shannon"],
                        "dominant_frac": me["dominant_frac"] - mf["dominant_frac"],
                        "total_tpm": float(np.log2((me["total_tpm"] + 1e-3) / (mf["total_tpm"] + 1e-3))),
                        "simplicity_S": me["simplicity_S"] - mf["simplicity_S"],
                    }
                )

    cov = pd.DataFrame(coverage_rows)
    cov.to_csv(RES / "oxphos_coverage_n7.tsv", sep="\t", index=False)
    if maps_rows:
        pd.concat(maps_rows, ignore_index=True).to_csv(RES / "oxphos_peak_symbol_map_n7.tsv", sep="\t", index=False)

    sens = pd.DataFrame(sens_rows)
    sens.to_csv(RES / "species_means_threshold_sensitivity.tsv", sep="\t", index=False)
    per_gene = pd.DataFrame(per_gene_rows)
    per_gene.to_csv(RES / "per_gene_metrics_thr1_samestrand.tsv", sep="\t", index=False)

    # species means for primary setting + FC tissue row means
    primary = sens[(sens.tpm_thr == 1.0) & (sens.same_strand)].copy()
    # add FC species means from per_gene
    fc = per_gene[per_gene.tissue == "log2_epi_over_fib"]
    if len(fc):
        fc_means = (
            fc.groupby(["species", "mass_g"], as_index=False)
            .agg(
                n_genes=("gene", "nunique"),
                n_active=("n_active", "mean"),
                shannon=("shannon", "mean"),
                dominant_frac=("dominant_frac", "mean"),
                total_tpm=("total_tpm", "mean"),
                simplicity_S=("simplicity_S", "mean"),
            )
        )
        fc_means["tissue"] = "log2_epi_over_fib"
        fc_means["tpm_thr"] = 1.0
        fc_means["same_strand"] = True
        primary = pd.concat([primary, fc_means], ignore_index=True)
    primary.to_csv(RES / "species_means_primary_thr1_samestrand.tsv", sep="\t", index=False)

    # OLS vs mass
    fit_rows = []
    for thr in thresholds:
        for same in (False, True):
            for tissue in tissues:
                sub = sens[(sens.tpm_thr == thr) & (sens.same_strand == same) & (sens.tissue == tissue)]
                if len(sub) < 5:
                    continue
                for col in ["shannon", "simplicity_S", "dominant_frac", "total_tpm", "n_active"]:
                    f = fit_mass(sub, col)
                    fit_rows.append(
                        {"tpm_thr": thr, "same_strand": same, "tissue": tissue, "metric": col, **f, "n_spp": len(sub)}
                    )
    # FC fits
    if len(fc):
        sub = primary[primary.tissue == "log2_epi_over_fib"]
        for col in ["shannon", "simplicity_S", "dominant_frac", "total_tpm", "n_active"]:
            f = fit_mass(sub, col)
            fit_rows.append(
                {
                    "tpm_thr": 1.0,
                    "same_strand": True,
                    "tissue": "log2_epi_over_fib",
                    "metric": col,
                    **f,
                    "n_spp": len(sub),
                }
            )
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(RES / "ols_vs_log10mass.tsv", sep="\t", index=False)

    # ---------- figures ----------
    # 1) coverage
    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    order = list(SPECIES.keys())
    covp = cov.set_index("species").loc[order]
    ax.bar(order, covp["n_oxphos_genes"], color="#1d3557")
    ax.set_ylabel("OXPHOS genes with promoter peaks")
    ax.set_title("Dropbox CAGE: OXPHOS coverage after Augustus→symbol map (n=7)")
    ax.tick_params(axis="x", rotation=30)
    for i, v in enumerate(covp["n_oxphos_genes"]):
        ax.text(i, v + 0.5, str(int(v)), ha="center", fontsize=8)
    fig.savefig(FIG / "fig_oxphos_promoter_coverage_n7.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_promoter_coverage_n7.png", dpi=160)
    plt.close(fig)

    # 2) primary architecture vs mass (fib vs epi)
    metrics_plot = [
        ("total_tpm", "mean total TPM"),
        ("shannon", "mean Shannon H"),
        ("dominant_frac", "mean dominant fraction"),
        ("simplicity_S", "mean simplicity S"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(12.5, 6.2), constrained_layout=True)
    colors = {"fib": "#e76f51", "epi": "#2a9d8f"}
    for col_i, (metric, ylab) in enumerate(metrics_plot):
        for row_i, tissue in enumerate(tissues):
            ax = axes[row_i, col_i]
            sub = primary[primary.tissue == tissue].sort_values("mass_g")
            if sub.empty:
                ax.set_axis_off()
                continue
            x = np.log10(sub.mass_g.values)
            y = sub[metric].values.astype(float)
            ax.scatter(x, y, c=colors[tissue], s=45, zorder=3)
            for _, r in sub.iterrows():
                ax.text(np.log10(r.mass_g), r[metric], r.species[:3], fontsize=7, ha="left", va="bottom")
            f = fit_mass(sub, metric)
            if np.isfinite(f["slope"]):
                xx = np.linspace(x.min(), x.max(), 50)
                ax.plot(xx, f["slope"] * xx + f["intercept"], color=colors[tissue], lw=1.2, alpha=0.8)
                ax.set_title(f"{tissue}: p={f['p']:.3f}, r={f['r']:.2f}", fontsize=9)
            ax.set_xlabel("log10 mass (g)")
            ax.set_ylabel(ylab if row_i == 0 else "")
            if col_i == 0:
                ax.set_ylabel(f"{tissue}\n{ylab}")
    fig.suptitle(
        "OXPHOS promoter architecture vs mass (Dropbox CAGE, thr=1 TPM, same-strand)\n"
        "Honest framing: report which endpoint moves; TPM may confound Shannon",
        fontsize=10,
    )
    fig.savefig(FIG / "fig_oxphos_promoter_arch_fib_epi_vs_mass.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_promoter_arch_fib_epi_vs_mass.png", dpi=160)
    plt.close(fig)

    # 3) threshold sensitivity p-values
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for ax, tissue, title in zip(axes, tissues, ["fibroblasts (C)", "epiblast (E)"]):
        for same, ls, label in [(False, "-", "all strands"), (True, "--", "same-strand")]:
            subH = fits[(fits.tissue == tissue) & (fits.same_strand == same) & (fits.metric == "shannon")].sort_values(
                "tpm_thr"
            )
            subT = fits[(fits.tissue == tissue) & (fits.same_strand == same) & (fits.metric == "total_tpm")].sort_values(
                "tpm_thr"
            )
            if len(subH):
                ax.plot(subH.tpm_thr, subH.p, "o" + ls, color="#1d3557", label=f"Shannon ({label})")
            if len(subT):
                ax.plot(subT.tpm_thr, subT.p, "s" + ls, color="#e76f51", label=f"total TPM ({label})")
        ax.axhline(0.05, color="#999", ls=":", lw=1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("TPM threshold")
        ax.set_ylabel("OLS p vs log10 mass")
        ax.set_title(title)
        ax.legend(fontsize=7, loc="upper right")
    fig.suptitle("Threshold / strand sensitivity — architecture vs expression confound", fontsize=10)
    fig.savefig(FIG / "fig_oxphos_promoter_threshold_sensitivity_fib_epi.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_promoter_threshold_sensitivity_fib_epi.png", dpi=160)
    plt.close(fig)

    # 4) epi vs fib: Shannon and TPM per species
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2), constrained_layout=True)
    pf = primary[primary.tissue == "fib"].set_index("species")
    pe = primary[primary.tissue == "epi"].set_index("species")
    common = [s for s in order if s in pf.index and s in pe.index]
    for ax, metric, title in zip(
        axes, ["shannon", "total_tpm"], ["mean Shannon H", "mean total TPM"]
    ):
        xf = pf.loc[common, metric].values
        xe = pe.loc[common, metric].values
        ax.scatter(xf, xe, c="#1d3557", s=50)
        for s, a, b in zip(common, xf, xe):
            ax.text(a, b, s[:3], fontsize=8)
        lims = [min(xf.min(), xe.min()), max(xf.max(), xe.max())]
        pad = 0.05 * (lims[1] - lims[0] + 1e-9)
        ax.plot([lims[0] - pad, lims[1] + pad], [lims[0] - pad, lims[1] + pad], ":", color="#999")
        ax.set_xlabel(f"fib {title}")
        ax.set_ylabel(f"epi {title}")
        ax.set_title(title)
    fig.suptitle("OXPHOS promoters: epiblast vs fibroblast (species means)", fontsize=10)
    fig.savefig(FIG / "fig_oxphos_promoter_epi_vs_fib.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_promoter_epi_vs_fib.png", dpi=160)
    plt.close(fig)

    # summary JSON
    def pick_fit(tissue: str, metric: str) -> dict:
        row = fits[(fits.tpm_thr == 1.0) & (fits.same_strand) & (fits.tissue == tissue) & (fits.metric == metric)]
        if row.empty:
            return {}
        r = row.iloc[0]
        return {"slope": r.slope, "r": r.r, "p": r.p, "n": int(r.n)}

    summary = {
        "n_species": 7,
        "primary_setting": "tpm_thr=1, same_strand=True, promoter-class peaks",
        "mapping": "chicken_genes + native geneId + DE genes; HGNC alias canon",
        "coverage": cov.to_dict(orient="records"),
        "ols_primary": {
            "fib_shannon": pick_fit("fib", "shannon"),
            "fib_total_tpm": pick_fit("fib", "total_tpm"),
            "fib_simplicity_S": pick_fit("fib", "simplicity_S"),
            "epi_shannon": pick_fit("epi", "shannon"),
            "epi_total_tpm": pick_fit("epi", "total_tpm"),
            "epi_simplicity_S": pick_fit("epi", "simplicity_S"),
            "fc_shannon_delta": pick_fit("log2_epi_over_fib", "shannon"),
            "fc_log2tpm": pick_fit("log2_epi_over_fib", "total_tpm"),
        },
        "interpretation_note": (
            "На этих данных: смотреть, двигается ли Shannon/S независимо от total TPM. "
            "Если TPM~mass значим, а Shannon нет (или Shannon исчезает при same-strand) — "
            "не заявлять упрощение архитектуры промоторов."
        ),
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary["ols_primary"], indent=2))
    print("Wrote", RES)


if __name__ == "__main__":
    main()
