#!/usr/bin/env python3
"""(a) Expand OXPHOS mapping / fair nuclear-only recompute; (b) confounder residuals.

a) mtDNA: verify presence of mito chroms in Dropbox CAGE tables; document gaps.
   Nuclear: remap with expanded aliases + all peak classes (not only promoter).
   Recompute fib/epi mean TPM ~ mass under several gene sets.

b) Confounders: CAGE depth, n_peaks, assembly fragmentation in CAGE tables,
   NCBI assembly stats, n_oxphos_genes. Correlate with TPM and with mass-residuals.
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
EXP = Path("/path/to/storage/kleiber_birds_gemcat/dropbox_20260808/unpacked/exp_tab")
EXP_LOCAL = ROOT / "data" / "dropbox_20260808"
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")
MITO = ROOT / "ref" / "mito"
RES = ROOT / "results" / "oxphos_ortholog_confounders"
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

MTDNA = {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8"}
KNOWN_MT_ACC = {
    "Zebra": {"NC_007897.1"},
    "Chicken": {"NC_001323.1"},  # may not appear if MT not in galGal CAGE table
    "Duck": set(),
    "Quail": set(),
    "Peacock": {"NC_024533.1", "KF444060.1"},
    "Turkey": {"NC_010195.2", "PP840326.1"},
    "Emu": {"NC_002784.1"},
}

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
    "COX4": "COX4I1",
}


def load_ox() -> set[str]:
    return {ln.strip().upper() for ln in (ROOT / "data/genesets/OXPHOS.txt").read_text().splitlines() if ln.strip()}


def is_augustus(s: str) -> bool:
    u = str(s).upper().strip()
    return bool(re.match(r"^G\d+(\.T\d+)?$", u)) or u in {"NAN", "NA", "", "NONE", "."}


def norm_sym(s: str, ox: set[str]) -> str:
    s = str(s).upper().strip()
    if not s or is_augustus(s) or s in {"NAN", "NA", "NONE", "."}:
        return ""
    s2 = PREFERRED_ALIAS.get(s, s)
    if s2 in ox:
        return s2
    if s in ox:
        return s
    return ""


def exp_path(fname: str) -> Path:
    for p in (EXP / fname, EXP_LOCAL / fname):
        if p.exists():
            return p
    raise FileNotFoundError(fname)


def peak_to_de(de_path: Path, ox: set[str]) -> pd.Series:
    de = pd.read_excel(de_path, usecols=["id", "genes", "geneId"])
    for c in ("genes", "geneId"):
        de[c] = de[c].astype(str)
    de["sym"] = de["genes"].where(~de["genes"].str.upper().isin(["NAN", "NA", ""]), de["geneId"])
    de["sym"] = de["sym"].map(lambda x: norm_sym(x, ox))
    de = de[de["sym"] != ""].drop_duplicates("id")
    return de.set_index("id")["sym"]


def annotate(df: pd.DataFrame, id2: pd.Series, ox: set[str]) -> pd.DataFrame:
    out = df.copy()
    gid = out["geneId"].astype(str)
    cg = out["chicken_genes"].astype(str) if "chicken_genes" in out.columns else pd.Series([""] * len(out))
    sym = []
    for g, c in zip(gid, cg):
        pick = ""
        if str(c) and not is_augustus(c) and str(c).upper() not in {"NAN", "NA"}:
            pick = c
        elif not is_augustus(g):
            pick = g
        sym.append(norm_sym(pick, ox))
    out["symbol"] = sym
    fill = out["symbol"].eq("")
    out.loc[fill, "symbol"] = out.loc[fill, "id"].map(id2).fillna("").map(lambda x: norm_sym(x, ox))
    return out


def lib_cols(df: pd.DataFrame, tissue: str) -> list[str]:
    pref = "norm.C" if tissue == "fib" else "norm.E"
    return [c for c in df.columns if str(c).startswith(pref)]


def raw_cols(df: pd.DataFrame, tissue: str) -> list[str]:
    pref = "raw.C" if tissue == "fib" else "raw.E"
    return [c for c in df.columns if str(c).startswith(pref)]


def scaffold_len(name: str) -> float:
    m = re.search(r"len(\d+)", str(name))
    return float(m.group(1)) if m else np.nan


def ols(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return {"n": int(m.sum()), "slope": np.nan, "r": np.nan, "p": np.nan, "intercept": np.nan}
    sl, ic, r, p, se = stats.linregress(x[m], y[m])
    return {"n": int(m.sum()), "slope": float(sl), "intercept": float(ic), "r": float(r), "p": float(p), "se": float(se)}


def gene_tpm_table(df: pd.DataFrame, tissue: str, genes: set[str], promoter_only: bool) -> pd.DataFrame:
    cols = lib_cols(df, tissue)
    if not cols:
        return pd.DataFrame(columns=["gene", "total_tpm", "n_peaks"])
    sub = df[df["symbol"].isin(genes)].copy()
    if promoter_only and "class" in sub.columns:
        sub = sub[sub["class"].astype(str).str.lower().eq("promoter")]
    if sub.empty:
        return pd.DataFrame(columns=["gene", "total_tpm", "n_peaks"])
    sub["_e"] = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
    g = sub.groupby("symbol", as_index=False).agg(total_tpm=("_e", "sum"), n_peaks=("_e", "size"))
    g = g.rename(columns={"symbol": "gene"})
    return g


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    ox = load_ox()
    nuclear = ox - MTDNA

    mito_rows = []
    cov_rows = []
    gene_tpm_rows = []
    conf_rows = []

    order = list(SPECIES.keys())

    for sp, meta in SPECIES.items():
        print(f"=== {sp} ===", flush=True)
        df = pd.read_csv(exp_path(meta["file"]), sep="\t")
        id2 = peak_to_de(RUSLAN / meta["de"], ox)
        df = annotate(df, id2, ox)

        chroms = df["seqnames"].astype(str)
        # mtDNA chrom detection
        known = KNOWN_MT_ACC.get(sp, set())
        hit_known = sorted(set(chroms) & known)
        mito_name = chroms.str.contains(r"mito|mitochond|chrM|\bMT\b", case=False, regex=True)
        # zebra-style NC_00 mito already known; also any chrom with MTDNA gene symbols
        mt_sym = df["symbol"].isin(MTDNA)
        mt_chroms = sorted(chroms[mt_sym].unique())
        mito_rows.append(
            {
                "species": sp,
                "n_chroms_cage": int(chroms.nunique()),
                "known_mt_accessions": ",".join(sorted(known)) if known else "",
                "known_mt_present_in_cage": ",".join(hit_known) if hit_known else "",
                "n_peaks_mito_name": int(mito_name.sum()),
                "n_peaks_mtdna_symbol": int(mt_sym.sum()),
                "mtdna_genes_found": ",".join(sorted(df.loc[mt_sym, "symbol"].unique())),
                "mtdna_chroms": ",".join(mt_chroms[:10]),
                "status": (
                    "mtdna_peaks_present"
                    if mt_sym.any() or hit_known
                    else "mtdna_absent_from_cage_assembly"
                ),
            }
        )

        # coverage under mapping modes
        for mode, genes, prom_only in [
            ("all_oxphos_promoter", ox, True),
            ("all_oxphos_anyclass", ox, False),
            ("nuclear_promoter", nuclear, True),
            ("nuclear_anyclass", nuclear, False),
        ]:
            gtab = gene_tpm_table(df, "fib", genes, prom_only)
            cov_rows.append(
                {
                    "species": sp,
                    "mass_g": meta["mass_g"],
                    "mode": mode,
                    "n_genes": int(gtab["gene"].nunique()) if len(gtab) else 0,
                    "n_peaks": int(gtab["n_peaks"].sum()) if len(gtab) else 0,
                    "mean_tpm": float(gtab["total_tpm"].mean()) if len(gtab) else np.nan,
                    "median_tpm": float(gtab["total_tpm"].median()) if len(gtab) else np.nan,
                    "sum_tpm": float(gtab["total_tpm"].sum()) if len(gtab) else np.nan,
                }
            )
            for tissue in ("fib", "epi"):
                gtab_t = gene_tpm_table(df, tissue, genes, prom_only)
                for _, r in gtab_t.iterrows():
                    gene_tpm_rows.append(
                        {
                            "species": sp,
                            "mass_g": meta["mass_g"],
                            "tissue": tissue,
                            "mode": mode,
                            "gene": r["gene"],
                            "total_tpm": r["total_tpm"],
                            "n_peaks": r["n_peaks"],
                            "is_mtdna": r["gene"] in MTDNA,
                        }
                    )

        # confounders from CAGE table (assembly actually used)
        rawC = raw_cols(df, "fib")
        rawE = raw_cols(df, "epi")
        slens = chroms.map(scaffold_len)
        # fallback chrom size = max end
        chrom_max = df.groupby("seqnames")["end"].max()
        med_chrom = float(np.nanmedian(slens.dropna())) if slens.notna().any() else float(chrom_max.median())
        conf_rows.append(
            {
                "species": sp,
                "mass_g": meta["mass_g"],
                "log10_mass": float(np.log10(meta["mass_g"])),
                "n_chroms_cage": int(chroms.nunique()),
                "n_peaks": int(len(df)),
                "median_chrom_len_proxy": med_chrom,
                "n_fib_libs": len(rawC),
                "n_epi_libs": len(rawE),
                "fib_total_raw": float(df[rawC].sum().sum()) if rawC else np.nan,
                "fib_mean_lib_raw": float(df[rawC].sum().mean()) if rawC else np.nan,
                "epi_total_raw": float(df[rawE].sum().sum()) if rawE else np.nan,
                "epi_mean_lib_raw": float(df[rawE].sum().mean()) if rawE else np.nan,
                "frac_peaks_with_symbol": float((df["symbol"] != "").mean()),
                "n_oxphos_nuclear_promoter": int(
                    df.loc[
                        df["symbol"].isin(nuclear)
                        & df["class"].astype(str).str.lower().eq("promoter"),
                        "symbol",
                    ].nunique()
                ),
                "n_oxphos_all_promoter": int(
                    df.loc[
                        df["symbol"].isin(ox) & df["class"].astype(str).str.lower().eq("promoter"),
                        "symbol",
                    ].nunique()
                ),
                "n_mtdna_genes": int(df.loc[df["symbol"].isin(MTDNA), "symbol"].nunique()),
            }
        )

    mito = pd.DataFrame(mito_rows)
    mito.to_csv(RES / "mtdna_presence_in_cage.tsv", sep="\t", index=False)
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(RES / "coverage_by_mapping_mode.tsv", sep="\t", index=False)
    gtp = pd.DataFrame(gene_tpm_rows)
    gtp.to_csv(RES / "gene_tpm_expanded_modes.tsv", sep="\t", index=False)
    conf = pd.DataFrame(conf_rows).set_index("species").loc[order].reset_index()

    # merge NCBI assembly stats if present
    ncbi_path = RES / "ncbi_assembly_stats.tsv"
    if ncbi_path.exists():
        ncbi = pd.read_csv(ncbi_path, sep="\t")
        conf = conf.merge(
            ncbi[
                [
                    c
                    for c in [
                        "species",
                        "assemblyaccession",
                        "total_length",
                        "scaffold_n50",
                        "contig_n50",
                        "contig_count",
                        "scaffold_count",
                    ]
                    if c in ncbi.columns
                ]
            ],
            on="species",
            how="left",
        )
        conf["note_ncbi_assembly"] = (
            "NCBI latest ≠ necessarily CAGE mapping assembly; prefer n_chroms_cage"
        )

    # species means for OLS
    fit_rows = []
    mean_rows = []
    for mode in sorted(cov["mode"].unique()):
        for tissue in ("fib", "epi"):
            sub = gtp[(gtp["mode"] == mode) & (gtp["tissue"] == tissue)]
            if sub.empty:
                continue
            means = sub.groupby(["species", "mass_g"], as_index=False)["total_tpm"].mean()
            means = means.set_index("species").reindex(order)
            means["mass_g"] = [SPECIES[s]["mass_g"] for s in order]
            means = means.reset_index()
            # drop species with no genes in this mode
            means = means.dropna(subset=["total_tpm"])
            if means.empty:
                continue
            means["tissue"] = tissue
            means["mode"] = mode
            mean_rows.append(means)
            f = ols(np.log10(means.mass_g), means.total_tpm)
            fit_rows.append({"tissue": tissue, "mode": mode, "endpoint": "mean_tpm", "n_genes": np.nan, **f})

    # matched gene sets
    for tissue in ("fib", "epi"):
        sub = gtp[(gtp["mode"] == "nuclear_anyclass") & (gtp["tissue"] == tissue)]
        wide = sub.pivot_table(index="gene", columns="species", values="total_tpm", aggfunc="mean")
        wide = wide.reindex(columns=order)
        for label, mat in [
            ("nuclear_anyclass_ge6", wide.loc[wide.notna().sum(axis=1) >= 6]),
            ("nuclear_anyclass_core7", wide.dropna(axis=0)),
            (
                "nuclear_anyclass_emu_anchored",
                wide.loc[wide["Emu"].notna()] if "Emu" in wide.columns else wide.iloc[0:0],
            ),
        ]:
            if len(mat) < 5:
                continue
            means = mat.mean(axis=0, skipna=True).reindex(order)
            masses = [SPECIES[s]["mass_g"] for s in order]
            f = ols(np.log10(masses), means.values)
            fit_rows.append(
                {"tissue": tissue, "mode": label, "endpoint": "mean_tpm", "n_genes": len(mat), **f}
            )
            mean_rows.append(
                pd.DataFrame(
                    {
                        "species": order,
                        "mass_g": masses,
                        "total_tpm": means.values,
                        "tissue": tissue,
                        "mode": label,
                    }
                )
            )

    fits = pd.DataFrame(fit_rows)
    fits.to_csv(RES / "ols_mean_tpm_vs_mass_by_mode.tsv", sep="\t", index=False)
    means_all = pd.concat(mean_rows, ignore_index=True)
    means_all.to_csv(RES / "species_mean_tpm_by_mode.tsv", sep="\t", index=False)

    # primary fib nuclear means into conf for residual analysis
    fib_nuc = means_all[(means_all["tissue"] == "fib") & (means_all["mode"] == "nuclear_anyclass")].set_index("species")
    conf = conf.set_index("species")
    conf["fib_mean_tpm_nuclear"] = fib_nuc.reindex(conf.index)["total_tpm"]
    fib_all = means_all[(means_all["tissue"] == "fib") & (means_all["mode"] == "all_oxphos_anyclass")].set_index("species")
    conf["fib_mean_tpm_all"] = fib_all.reindex(conf.index)["total_tpm"]
    epi_nuc = means_all[(means_all["tissue"] == "epi") & (means_all["mode"] == "nuclear_anyclass")].set_index("species")
    conf["epi_mean_tpm_nuclear"] = epi_nuc.reindex(conf.index)["total_tpm"]

    # residual after mass
    x = conf["log10_mass"].values
    y = conf["fib_mean_tpm_nuclear"].values
    fmass = ols(x, y)
    yhat = fmass["slope"] * x + fmass["intercept"]
    conf["fib_tpm_resid_after_mass"] = y - yhat
    conf = conf.reset_index()
    conf.to_csv(RES / "confounders_by_species.tsv", sep="\t", index=False)

    # correlate confounders with TPM and with residual
    candidates = [
        "log10_mass",
        "n_oxphos_nuclear_promoter",
        "n_oxphos_all_promoter",
        "n_mtdna_genes",
        "n_chroms_cage",
        "n_peaks",
        "median_chrom_len_proxy",
        "fib_mean_lib_raw",
        "fib_total_raw",
        "frac_peaks_with_symbol",
        "total_length",
        "scaffold_n50",
        "contig_count",
    ]
    corr_rows = []
    for col in candidates:
        if col not in conf.columns:
            continue
        for ycol, ylab in [
            ("fib_mean_tpm_nuclear", "fib_tpm"),
            ("fib_tpm_resid_after_mass", "fib_tpm_resid"),
            ("epi_mean_tpm_nuclear", "epi_tpm"),
        ]:
            f = ols(conf[col], conf[ycol])
            # also Spearman
            m = np.isfinite(conf[col]) & np.isfinite(conf[ycol])
            if m.sum() >= 5:
                rho, sp = stats.spearmanr(conf.loc[m, col], conf.loc[m, ycol])
            else:
                rho, sp = np.nan, np.nan
            corr_rows.append(
                {
                    "predictor": col,
                    "response": ylab,
                    "pearson_r": f["r"],
                    "pearson_p": f["p"],
                    "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                    "spearman_p": float(sp) if np.isfinite(sp) else np.nan,
                    "n": f["n"],
                }
            )
    corrs = pd.DataFrame(corr_rows)
    corrs.to_csv(RES / "confounder_correlations.tsv", sep="\t", index=False)

    # partial: TPM ~ mass + n_oxphos + depth (simple sequential)
    # compare nested models via F? keep simple: report pairwise + mass-controlled residual corrs
    # multivariate OLS with available predictors
    pred_mv = ["log10_mass", "n_oxphos_nuclear_promoter", "n_chroms_cage", "fib_mean_lib_raw"]
    X = conf[pred_mv].to_numpy(dtype=float)
    yv = conf["fib_mean_tpm_nuclear"].to_numpy(dtype=float)
    # add intercept
    ok = np.isfinite(X).all(axis=1) & np.isfinite(yv)
    Xo = np.column_stack([np.ones(ok.sum()), X[ok]])
    yo = yv[ok]
    beta, *_ = np.linalg.lstsq(Xo, yo, rcond=None)
    yhat_mv = Xo @ beta
    ss_res = float(((yo - yhat_mv) ** 2).sum())
    ss_tot = float(((yo - yo.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    # univariate drop-one r2
    drop = {}
    for i, name in enumerate(pred_mv):
        cols = [j for j in range(X.shape[1]) if j != i]
        Xd = np.column_stack([np.ones(ok.sum()), X[ok][:, cols]])
        b, *_ = np.linalg.lstsq(Xd, yo, rcond=None)
        yh = Xd @ b
        r2d = 1 - float(((yo - yh) ** 2).sum()) / ss_tot
        drop[name] = {"r2_without": r2d, "delta_r2": r2 - r2d}

    mv = {
        "predictors": pred_mv,
        "beta_intercept": float(beta[0]),
        "betas": {n: float(b) for n, b in zip(pred_mv, beta[1:])},
        "r2": float(r2),
        "drop_one": drop,
        "n": int(ok.sum()),
    }
    (RES / "multivariate_fib_tpm.json").write_text(json.dumps(mv, indent=2))

    # ---------- figures ----------
    # 1) mtDNA presence + coverage modes
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    ax = axes[0]
    colors = ["#2a9d8f" if s == "mtdna_peaks_present" else "#e76f51" for s in mito.set_index("species").loc[order, "status"]]
    ax.bar(order, mito.set_index("species").loc[order, "n_peaks_mtdna_symbol"], color=colors)
    ax.set_title("mtDNA-symbol peaks in CAGE tables")
    ax.set_ylabel("n peaks")
    ax.tick_params(axis="x", rotation=35)
    ax = axes[1]
    for mode, sty, col in [
        ("all_oxphos_anyclass", "o-", "#1d3557"),
        ("nuclear_anyclass", "s-", "#e76f51"),
        ("nuclear_anyclass_ge6", "^-", "#2a9d8f"),
    ]:
        sub = means_all[(means_all["tissue"] == "fib") & (means_all["mode"] == mode)].copy()
        if sub.empty:
            continue
        sub = sub.set_index("species").reindex(order).dropna(subset=["total_tpm"])
        ax.plot(np.log10(sub.mass_g), sub.total_tpm, sty, color=col, label=mode)
    ax.legend(fontsize=7)
    ax.set_xlabel("log10 mass (g)")
    ax.set_ylabel("mean gene TPM (fib)")
    ax.set_title("Recompute after nuclear-only / matched sets")
    fig.suptitle("(a) mtDNA gap + expanded nuclear recompute", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_mtdna_gap_and_recompute.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_mtdna_gap_and_recompute.png", dpi=160)
    plt.close(fig)

    # 2) confounder panel
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.8), constrained_layout=True)
    plots = [
        ("n_oxphos_nuclear_promoter", "n nuclear OXPHOS genes"),
        ("n_chroms_cage", "n chroms in CAGE assembly"),
        ("fib_mean_lib_raw", "mean fib CAGE depth (raw)"),
        ("median_chrom_len_proxy", "median chrom len proxy"),
        ("scaffold_n50", "NCBI scaffold N50"),
        ("total_length", "NCBI genome length"),
    ]
    for ax, (col, title) in zip(axes.ravel(), plots):
        if col not in conf.columns or conf[col].isna().all():
            ax.set_axis_off()
            continue
        ax.scatter(conf[col], conf["fib_mean_tpm_nuclear"], c="#1d3557", s=45)
        for _, r in conf.iterrows():
            if np.isfinite(r[col]) and np.isfinite(r["fib_mean_tpm_nuclear"]):
                ax.text(r[col], r["fib_mean_tpm_nuclear"], r.species[:3], fontsize=7)
        f = ols(conf[col], conf["fib_mean_tpm_nuclear"])
        ax.set_title(f"{title}\nr={f['r']:.2f} p={f['p']:.3f}", fontsize=9)
        ax.set_ylabel("fib mean nuclear TPM")
        ax.set_xlabel(col)
    fig.suptitle("(b) Does fib OXPHOS TPM track confounders?", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_tpm_vs_confounders.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_tpm_vs_confounders.png", dpi=160)
    plt.close(fig)

    # 3) residual after mass vs key confounders
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8), constrained_layout=True)
    for ax, col in zip(axes, ["n_oxphos_nuclear_promoter", "n_chroms_cage", "fib_mean_lib_raw"]):
        ax.scatter(conf[col], conf["fib_tpm_resid_after_mass"], c="#e76f51", s=50)
        for _, r in conf.iterrows():
            ax.text(r[col], r["fib_tpm_resid_after_mass"], r.species[:3], fontsize=7)
        f = ols(conf[col], conf["fib_tpm_resid_after_mass"])
        ax.axhline(0, color="#999", ls=":", lw=1)
        ax.set_xlabel(col)
        ax.set_ylabel("TPM residual after mass")
        ax.set_title(f"r={f['r']:.2f} p={f['p']:.3f}")
    fig.suptitle("Residuals after regressing out log10(mass)", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_tpm_residuals_vs_confounders.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_tpm_residuals_vs_confounders.png", dpi=160)
    plt.close(fig)

    # 4) heatmap of correlations
    pivot = corrs.pivot(index="predictor", columns="response", values="pearson_r")
    fig, ax = plt.subplots(figsize=(6.2, 5.5), constrained_layout=True)
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title("Pearson r: predictors vs TPM / residual")
    fig.savefig(FIG / "fig_oxphos_confounder_corr_heatmap.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_confounder_corr_heatmap.png", dpi=160)
    plt.close(fig)

    summary = {
        "mtdna": mito.to_dict(orient="records"),
        "ols_key": fits.to_dict(orient="records"),
        "mass_fit_fib_nuclear": fmass,
        "multivariate": mv,
        "top_residual_corrs": corrs[corrs.response == "fib_tpm_resid"]
        .sort_values("pearson_p")
        .head(8)
        .to_dict(orient="records"),
        "conclusions": [
            "mtDNA CAGE peaks отсутствуют у Peacock/Turkey/Emu (сборка без MT) — это не failure ортолог-маппинга ядерных генов.",
            "Fair recompute: nuclear-only / ge6 / emu-anchored — см. ols table.",
            "Confounders: смотреть residual-after-mass; если n_oxphos_genes остаётся главным коррелятом residual — coverage всё ещё двигает сигнал.",
        ],
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))

    print("\n=== mtDNA status ===")
    print(mito[["species", "status", "n_peaks_mtdna_symbol", "mtdna_genes_found"]].to_string(index=False))
    print("\n=== OLS fib ===")
    print(fits[fits.tissue == "fib"][["mode", "slope", "r", "p"]].to_string(index=False))
    print("\n=== OLS epi ===")
    print(fits[fits.tissue == "epi"][["mode", "slope", "r", "p"]].to_string(index=False))
    print("\n=== residual corrs ===")
    print(corrs[corrs.response == "fib_tpm_resid"].sort_values("pearson_p").head(10).to_string(index=False))
    print("\n=== multivariate ===")
    print(json.dumps(mv, indent=2))
    print("Wrote", RES)


if __name__ == "__main__":
    main()
