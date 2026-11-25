#!/usr/bin/env python3
"""Which OXPHOS genes/complexes drive fib TPM↓ with mass; multi-promoter context.

Uses per-gene metrics from results/oxphos_promoter_architecture_n7
(primary: thr=1 TPM, same-strand).

Outputs:
  - per-gene OLS vs log10(mass) for fib / epi
  - complex-level summaries
  - contribution of top drivers to species mean TPM
  - multi-promoter status of drivers
  - figures + short why-hypotheses JSON
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
RES_IN = ROOT / "results" / "oxphos_promoter_architecture_n7"
RES = ROOT / "results" / "oxphos_gene_drivers_vs_mass"
FIG = ROOT / "figures"

MTDNA = {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8"}


def complex_of(gene: str) -> str:
    g = gene.upper()
    if g in MTDNA:
        if g.startswith("ND"):
            return "CI_mtDNA"
        if g.startswith("COX"):
            return "CIV_mtDNA"
        if g.startswith("ATP"):
            return "CV_mtDNA"
        if g == "CYTB":
            return "CIII_mtDNA"
        return "mtDNA_other"
    if g.startswith("NDUF") or g in {"NDUFA1", "NDUFA10"}:
        return "CI_nuclear"
    if g.startswith("SDH"):
        return "CII"
    if g.startswith("UQCR") or g in {"CYC1", "CYCS"}:
        return "CIII_nuclear"
    if g.startswith("COX"):
        return "CIV_nuclear"
    if g.startswith("ATP5") or g.startswith("ATP"):
        return "CV_nuclear"
    return "other_OXPHOS"


def fit_gene(df: pd.DataFrame, ycol: str) -> dict:
    x = np.log10(df["mass_g"].values.astype(float))
    y = df[ycol].values.astype(float)
    m = np.isfinite(x) & np.isfinite(y) & (y > 0)
    if m.sum() < 5:
        return {"n_spp": int(m.sum()), "slope": np.nan, "r": np.nan, "p": np.nan}
    # use log10 TPM for scale-fair slopes
    y2 = np.log10(y[m])
    slope, icept, r, p, se = stats.linregress(x[m], y2)
    return {
        "n_spp": int(m.sum()),
        "slope_log10tpm": float(slope),
        "intercept": float(icept),
        "r": float(r),
        "p": float(p),
        "se": float(se),
        "mean_tpm": float(y[m].mean()),
        "median_tpm": float(np.median(y[m])),
    }


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    pg = pd.read_csv(RES_IN / "per_gene_metrics_thr1_samestrand.tsv", sep="\t")
    pg = pg[pg.tissue.isin(["fib", "epi"])].copy()
    pg["complex"] = pg["gene"].map(complex_of)

    gene_rows = []
    for tissue in ("fib", "epi"):
        sub = pg[pg.tissue == tissue]
        for gene, gdf in sub.groupby("gene"):
            f = fit_gene(gdf, "total_tpm")
            # architecture means
            gene_rows.append(
                {
                    "tissue": tissue,
                    "gene": gene,
                    "complex": complex_of(gene),
                    "mean_n_active": float(gdf["n_active"].mean()),
                    "mean_shannon": float(gdf["shannon"].mean()),
                    "mean_dominant_frac": float(gdf["dominant_frac"].mean()),
                    "frac_multipromoter": float((gdf["n_active"] >= 2).mean()),
                    **f,
                }
            )
    genes = pd.DataFrame(gene_rows)
    genes["fdr_bh"] = np.nan
    for tissue in ("fib", "epi"):
        m = genes.tissue == tissue
        p = genes.loc[m, "p"].to_numpy(dtype=float)
        ok = np.isfinite(p)
        q = np.full(p.shape, np.nan, dtype=float)
        if ok.any():
            pv = p[ok]
            order = np.argsort(pv)
            n = len(pv)
            bh = pv[order] * n / np.arange(1, n + 1)
            bh = np.minimum.accumulate(bh[::-1])[::-1]
            bh = np.clip(bh, 0, 1)
            qv = np.empty(n, dtype=float)
            qv[order] = bh
            q[ok] = qv
        genes.loc[m, "fdr_bh"] = q

    genes.to_csv(RES / "per_gene_ols_log10tpm_vs_mass.tsv", sep="\t", index=False)

    fib = genes[(genes.tissue == "fib") & (genes.n_spp >= 5)].sort_values("p")
    epi = genes[(genes.tissue == "epi") & (genes.n_spp >= 5)].sort_values("p")
    fib.to_csv(RES / "fib_gene_ranked.tsv", sep="\t", index=False)
    epi.to_csv(RES / "epi_gene_ranked.tsv", sep="\t", index=False)

    # complex-level: mean of per-gene slopes (genes with n>=5)
    crows = []
    for tissue, gtab in [("fib", fib), ("epi", epi)]:
        for cx, cdf in gtab.groupby("complex"):
            crows.append(
                {
                    "tissue": tissue,
                    "complex": cx,
                    "n_genes": len(cdf),
                    "median_slope": float(cdf.slope_log10tpm.median()),
                    "mean_slope": float(cdf.slope_log10tpm.mean()),
                    "frac_neg_slope": float((cdf.slope_log10tpm < 0).mean()),
                    "n_p_lt_0.05": int((cdf.p < 0.05).sum()),
                    "n_fdr_lt_0.1": int((cdf.fdr_bh < 0.1).sum()),
                    "median_mean_tpm": float(cdf.mean_tpm.median()),
                    "mean_frac_multipromoter": float(cdf.frac_multipromoter.mean()),
                }
            )
    complexes = pd.DataFrame(crows).sort_values(["tissue", "median_slope"])
    complexes.to_csv(RES / "complex_summary.tsv", sep="\t", index=False)

    # contribution: species mean TPM = mean over genes present.
    # Leave-one-out / top-driver share of cross-species variance
    fib_pg = pg[pg.tissue == "fib"]
    # wide TPM
    wide = fib_pg.pivot_table(index="gene", columns="species", values="total_tpm", aggfunc="mean")
    # species order by mass
    mass = fib_pg.groupby("species")["mass_g"].first().sort_values()
    spp = list(mass.index)
    wide = wide.reindex(columns=spp)

    # correlation of each gene with log10 mass (same genes as OLS)
    contrib = []
    x = np.log10(mass.values)
    sp_mean = wide.mean(axis=0, skipna=True)
    # variance of species means explained by projecting top genes
    for gene in wide.index:
        y = wide.loc[gene, spp].values.astype(float)
        m = np.isfinite(y) & (y > 0)
        if m.sum() < 5:
            continue
        r, p = stats.pearsonr(x[m], np.log10(y[m]))
        # weight by mean TPM (how much it pulls the module mean)
        contrib.append(
            {
                "gene": gene,
                "complex": complex_of(gene),
                "r_log10tpm_mass": float(r),
                "p": float(p),
                "mean_tpm": float(np.nanmean(y)),
                "weight_x_neg_r": float(np.nanmean(y) * max(-r, 0)),
                "n_spp": int(m.sum()),
            }
        )
    contrib = pd.DataFrame(contrib).sort_values("weight_x_neg_r", ascending=False)
    contrib["cum_weight_frac"] = contrib["weight_x_neg_r"].cumsum() / contrib["weight_x_neg_r"].sum()
    contrib.to_csv(RES / "fib_driver_contribution.tsv", sep="\t", index=False)

    # top drivers table join architecture
    top = contrib.head(20).merge(
        fib[["gene", "slope_log10tpm", "p", "fdr_bh", "mean_n_active", "mean_shannon", "frac_multipromoter"]],
        on="gene",
        how="left",
        suffixes=("", "_ols"),
    )
    top.to_csv(RES / "fib_top20_drivers.tsv", sep="\t", index=False)

    # nuclear vs mtDNA
    nuc_mt = []
    for tissue, gtab in [("fib", fib), ("epi", epi)]:
        gtab = gtab.copy()
        gtab["genome"] = np.where(gtab.complex.str.contains("mtDNA"), "mtDNA", "nuclear")
        for genome, cdf in gtab.groupby("genome"):
            nuc_mt.append(
                {
                    "tissue": tissue,
                    "genome": genome,
                    "n_genes": len(cdf),
                    "median_slope": float(cdf.slope_log10tpm.median()),
                    "mean_slope": float(cdf.slope_log10tpm.mean()),
                    "frac_p_lt_0.05": float((cdf.p < 0.05).mean()),
                    "median_tpm": float(cdf.mean_tpm.median()),
                }
            )
    pd.DataFrame(nuc_mt).to_csv(RES / "nuclear_vs_mtdna.tsv", sep="\t", index=False)

    # ---------- figures ----------
    # 1) volcano fib
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), constrained_layout=True)
    for ax, tissue, tab in zip(axes, ["fib", "epi"], [fib, epi]):
        ax.scatter(tab.slope_log10tpm, -np.log10(tab.p.clip(1e-6)), c="#1d3557", s=28, alpha=0.75)
        sig = tab[tab.p < 0.05]
        if len(sig):
            ax.scatter(sig.slope_log10tpm, -np.log10(sig.p.clip(1e-6)), c="#e76f51", s=36)
            for _, r in sig.nsmallest(12, "p").iterrows():
                ax.text(r.slope_log10tpm, -np.log10(max(r.p, 1e-6)), r.gene, fontsize=7)
        ax.axvline(0, color="#999", ls=":", lw=1)
        ax.axhline(-np.log10(0.05), color="#999", ls=":", lw=1)
        ax.set_xlabel("slope log10(TPM) ~ log10(mass)")
        ax.set_ylabel("−log10 p")
        ax.set_title(f"{tissue}: per-gene OXPHOS (n_spp≥5)")
    fig.suptitle("Which OXPHOS genes decline with mass?", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_gene_volcano_fib_epi.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_gene_volcano_fib_epi.png", dpi=160)
    plt.close(fig)

    # 2) complex box of slopes
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    order_cx = [
        "CI_nuclear",
        "CI_mtDNA",
        "CII",
        "CIII_nuclear",
        "CIII_mtDNA",
        "CIV_nuclear",
        "CIV_mtDNA",
        "CV_nuclear",
        "CV_mtDNA",
        "other_OXPHOS",
    ]
    for ax, tissue, tab in zip(axes, ["fib", "epi"], [fib, epi]):
        data, labels = [], []
        for cx in order_cx:
            vals = tab.loc[tab.complex == cx, "slope_log10tpm"].dropna().values
            if len(vals):
                data.append(vals)
                labels.append(f"{cx}\n(n={len(vals)})")
        if data:
            bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=True)
            for patch in bp["boxes"]:
                patch.set_facecolor("#a8dadc")
            ax.axhline(0, color="#999", ls=":", lw=1)
            ax.tick_params(axis="x", labelsize=7, rotation=35)
            ax.set_ylabel("slope log10(TPM)~log10(mass)")
            ax.set_title(tissue)
    fig.suptitle("OXPHOS complex slopes vs mass", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_complex_slopes_fib_epi.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_complex_slopes_fib_epi.png", dpi=160)
    plt.close(fig)

    # 3) top drivers heatmap (genes x species), z-scored log10 TPM
    top_genes = list(contrib.head(15)["gene"])
    mat = wide.loc[[g for g in top_genes if g in wide.index], spp]
    logm = np.log10(mat.clip(lower=1e-3))
    z = logm.sub(logm.mean(axis=1), axis=0).div(logm.std(axis=1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(8.2, 5.5), constrained_layout=True)
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(spp)))
    ax.set_xticklabels([f"{s}\n{mass[s]:.0f}g" for s in spp], fontsize=8)
    ax.set_yticks(range(len(z)))
    ax.set_yticklabels(z.index, fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.8, label="z(log10 TPM)")
    ax.set_title("Top fib drivers (by meanTPM × negative r with mass)")
    fig.savefig(FIG / "fig_oxphos_top_drivers_heatmap_fib.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_top_drivers_heatmap_fib.png", dpi=160)
    plt.close(fig)

    # 4) multipromoter vs slope
    fig, ax = plt.subplots(figsize=(5.8, 4.2), constrained_layout=True)
    ax.scatter(fib.frac_multipromoter, fib.slope_log10tpm, c="#1d3557", s=32, alpha=0.75)
    # highlight top drivers
    td = fib[fib.gene.isin(top_genes)]
    ax.scatter(td.frac_multipromoter, td.slope_log10tpm, c="#e76f51", s=50, label="top15 drivers")
    for _, r in td.iterrows():
        ax.text(r.frac_multipromoter, r.slope_log10tpm, r.gene, fontsize=7)
    ax.axhline(0, color="#999", ls=":", lw=1)
    ax.set_xlabel("fraction species with ≥2 active promoters")
    ax.set_ylabel("slope log10(TPM)~log10(mass)")
    ax.set_title("fib: multipromoter status vs mass-slope")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "fig_oxphos_multipromoter_vs_slope_fib.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_multipromoter_vs_slope_fib.png", dpi=160)
    plt.close(fig)

    # 5) species mean TPM: all vs without top5 drivers
    top5 = list(contrib.head(5)["gene"])
    rest = [g for g in wide.index if g not in top5]
    means_all = wide.mean(axis=0, skipna=True)
    means_rest = wide.loc[rest].mean(axis=0, skipna=True) if rest else means_all
    means_top5 = wide.loc[[g for g in top5 if g in wide.index]].mean(axis=0, skipna=True)
    fig, ax = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
    xx = np.log10(mass.values)
    ax.plot(xx, means_all.values, "o-", label="all OXPHOS genes", color="#1d3557")
    ax.plot(xx, means_rest.values, "s--", label="without top5 drivers", color="#2a9d8f")
    ax.plot(xx, means_top5.values, "^:", label="top5 drivers only", color="#e76f51")
    for s, xi, yi in zip(spp, xx, means_all.values):
        ax.text(xi, yi, s[:3], fontsize=7)
    ax.set_xlabel("log10 mass (g)")
    ax.set_ylabel("mean gene TPM (fib)")
    ax.set_title("Is the mass trend carried by a few genes?")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "fig_oxphos_leave_top5_species_mean_fib.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_leave_top5_species_mean_fib.png", dpi=160)
    plt.close(fig)

    # stats for leave-top5
    def ols(y):
        y = np.asarray(y, float)
        m = np.isfinite(y)
        sl, ic, r, p, se = stats.linregress(xx[m], y[m])
        return {"slope": float(sl), "r": float(r), "p": float(p)}

    leave = {
        "all": ols(means_all.values),
        "without_top5": ols(means_rest.values),
        "top5_only": ols(means_top5.values),
        "top5_genes": top5,
    }

    # fib vs epi slopes for same genes
    merged = fib.merge(epi, on="gene", suffixes=("_fib", "_epi"))
    merged.to_csv(RES / "fib_vs_epi_gene_slopes.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(5.5, 5.2), constrained_layout=True)
    ax.scatter(merged.slope_log10tpm_fib, merged.slope_log10tpm_epi, c="#1d3557", s=30, alpha=0.7)
    hi = merged[merged.gene.isin(top_genes)]
    ax.scatter(hi.slope_log10tpm_fib, hi.slope_log10tpm_epi, c="#e76f51", s=50)
    for _, r in hi.iterrows():
        ax.text(r.slope_log10tpm_fib, r.slope_log10tpm_epi, r.gene, fontsize=7)
    lim = np.nanmax(np.abs(merged[["slope_log10tpm_fib", "slope_log10tpm_epi"]].values)) * 1.05
    ax.plot([-lim, lim], [-lim, lim], ":", color="#999")
    ax.axhline(0, color="#ccc", lw=0.8)
    ax.axvline(0, color="#ccc", lw=0.8)
    ax.set_xlabel("fib slope")
    ax.set_ylabel("epi slope")
    ax.set_title("Same gene: mass-slope in fib vs epi")
    fig.savefig(FIG / "fig_oxphos_gene_slope_fib_vs_epi.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_gene_slope_fib_vs_epi.png", dpi=160)
    plt.close(fig)

    # why hypotheses from patterns
    n_fib_sig = int((fib.p < 0.05).sum())
    n_fib_fdr = int((fib.fdr_bh < 0.1).sum())
    n_epi_sig = int((epi.p < 0.05).sum())
    med_mp_drivers = float(top["frac_multipromoter"].median()) if len(top) else np.nan
    med_mp_all = float(fib["frac_multipromoter"].median())
    frac_neg = float((fib.slope_log10tpm < 0).mean())

    summary = {
        "n_genes_fib_nspp_ge5": int(len(fib)),
        "n_fib_p_lt_0.05": n_fib_sig,
        "n_fib_fdr_lt_0.1": n_fib_fdr,
        "n_epi_p_lt_0.05": n_epi_sig,
        "frac_fib_genes_neg_slope": frac_neg,
        "top5_drivers": top5,
        "leave_top5_species_mean_ols": leave,
        "median_multipromoter_frac_top20": med_mp_drivers,
        "median_multipromoter_frac_all": med_mp_all,
        "nuclear_vs_mtdna": pd.DataFrame(nuc_mt).to_dict(orient="records"),
        "complex_fib": complexes[complexes.tissue == "fib"].to_dict(orient="records"),
        "hypotheses_ranked": [
            {
                "id": "H1_culture_fib_expression",
                "claim": "Сигнал в основном про уровень OXPHOS TPM в culture fibroblasts, не про архитектуру промоторов.",
                "support": "Shannon/S null; fib TPM~mass сильный; большинство генов single-promoter (dominant_frac~1).",
            },
            {
                "id": "H2_broad_nuclear_not_one_complex",
                "claim": "Если негативный slope распространён по ядерным комплексам (не только mtDNA) — это скорее клеточный/метаболический setpoint культуры, не специфичный ЭТЦ-комплекс.",
                "support": "см. complex_summary + frac_neg_slope",
            },
            {
                "id": "H3_epi_buffering",
                "claim": "Эпибласт буферизует: те же гены не дают TPM~mass — регуляция OXPHOS в pluripotent-like состоянии менее завязана на body mass proxy.",
                "support": "epi gene volcano слабее; fib vs epi slope scatter уходит от диагонали для drivers",
            },
            {
                "id": "H4_not_promoter_switching",
                "claim": "Drivers не объясняются переключением промоторов: frac multipromoter у top drivers ≈ как у всех.",
                "support": "fig multipromoter vs slope",
            },
            {
                "id": "H5_annotation_coverage_confound",
                "claim": "Частичный confound: у тяжёлых видов (peacock/emu) меньше OXPHOS-генов в аннотации → module mean может смещаться, если missing genes skewed.",
                "support": "coverage 31–36 у peacock/emu vs 60+ у мелких; проверить leave-species и common-gene core",
            },
        ],
    }

    # common-gene core check: genes present in all 7 spp
    core_genes = wide.dropna(axis=0).index.tolist()
    if len(core_genes) >= 10:
        core_mean = wide.loc[core_genes, spp].mean(axis=0)
        summary["common_core"] = {"n_genes": len(core_genes), "ols": ols(core_mean.values), "genes_head": core_genes[:20]}
    else:
        # genes in >=6
        present_n = wide.notna().sum(axis=1)
        core6 = wide.index[present_n >= 6].tolist()
        core_mean = wide.loc[core6, spp].mean(axis=0, skipna=True)
        summary["common_core"] = {"n_genes_ge6": len(core6), "ols": ols(core_mean.values)}

    (RES / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    print("=== TOP fib drivers ===")
    print(top[["gene", "complex", "mean_tpm", "r_log10tpm_mass", "p", "frac_multipromoter", "cum_weight_frac"]].head(15).to_string(index=False))
    print("\n=== complex fib ===")
    print(complexes[complexes.tissue == "fib"].to_string(index=False))
    print("\n=== leave top5 ===")
    print(json.dumps(leave, indent=2))
    print("\n=== common core ===")
    print(json.dumps(summary["common_core"], indent=2))
    print("Wrote", RES)


if __name__ == "__main__":
    main()
