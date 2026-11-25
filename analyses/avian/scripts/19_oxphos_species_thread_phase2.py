#!/usr/bin/env python3
"""Phase 2: fib+epi × set/module/regulon — Y1/Y2/Y3, matched, residual, coherence, coupling."""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from oxphos_thread_lib import (
    DRIVERS,
    EXP,
    FIG,
    MATCHED_HIGH_CODES,
    RES,
    RUSLAN,
    SALMON,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    filter_matched_cols,
    gene_by_lib_matrix,
    lib_cols,
    load_geneset,
    median_pairwise_spearman,
    ols_logmass,
    peak_to_symbol_de,
    residual_after_covariate,
)

warnings.filterwarnings("ignore", category=UserWarning)

DEF = RES / "definitions"
OUT = RES / "phase2"


def load_def(name: str) -> set[str]:
    return set(pd.read_csv(DEF / f"oxphos_{name}.tsv", sep="\t")["gene"].astype(str).str.upper())


def load_species_tables() -> dict[str, pd.DataFrame]:
    out = {}
    for sp, meta in SPECIES.items():
        df = pd.read_csv(EXP / meta["file"], sep="\t")
        id2 = peak_to_symbol_de(RUSLAN / meta["de"])
        out[sp] = annotate_table(df, id2)
    return out


def global_promoter_mean(df: pd.DataFrame, tissue: str, cols: list[str] | None = None) -> float:
    cols = cols or lib_cols(df, tissue)
    if not cols:
        return float("nan")
    sub = df[df["class"].astype(str).str.lower().eq("promoter")]
    if sub.empty:
        return float("nan")
    m = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return float(m.to_numpy().mean())


def global_promoter_sum(df: pd.DataFrame, tissue: str, cols: list[str] | None = None) -> float:
    cols = cols or lib_cols(df, tissue)
    if not cols:
        return float("nan")
    sub = df[df["class"].astype(str).str.lower().eq("promoter")]
    m = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # mean across libs of sum of peaks
    return float(m.sum(axis=0).mean())


def y_metrics(df: pd.DataFrame, genes: set[str], tissue: str, cols: list[str] | None = None) -> dict:
    cols = cols or lib_cols(df, tissue)
    mat = gene_by_lib_matrix(df, genes, tissue, cols=cols)
    mat = mat.loc[(mat > 0).any(axis=1)]
    if mat.empty:
        return dict(n_genes=0, Y1=np.nan, Y2=np.nan, gene_sum=np.nan, global_sum=np.nan)
    # per-gene mean across libs, then mean across genes
    gene_means = mat.mean(axis=1)
    y1 = float(gene_means.mean())
    gene_sum = float(mat.sum(axis=0).mean())  # mean over libs of sum of gene totals
    gsum = global_promoter_sum(df, tissue, cols)
    y2 = gene_sum / gsum if gsum and np.isfinite(gsum) and gsum > 0 else np.nan
    return dict(n_genes=int(len(mat)), Y1=y1, Y2=y2, gene_sum=gene_sum, global_sum=gsum)


def salmon_y3(genes: set[str]) -> pd.Series:
    sp = pd.read_csv(SALMON, sep="\t", index_col=0)
    # columns are species names matching salmon keys
    rows = {}
    for cage_sp, meta in SPECIES.items():
        col = meta["salmon"]
        if col not in sp.columns:
            rows[cage_sp] = np.nan
            continue
        present = [g for g in genes if g in sp.index]
        if not present:
            rows[cage_sp] = np.nan
            continue
        vals = np.log10(sp.loc[present, col].astype(float).clip(lower=0) + 1)
        rows[cage_sp] = float(vals.mean())
    return pd.Series(rows)


def coherence_species(df: pd.DataFrame, genes: set[str], tissue: str) -> float:
    mat = gene_by_lib_matrix(df, genes, tissue)
    mat = mat.loc[(mat > 0).any(axis=1)]
    return median_pairwise_spearman(mat)


def coupling_species(df: pd.DataFrame, genes: set[str], tissue: str) -> float:
    mat = gene_by_lib_matrix(df, genes, tissue)
    mat = mat.loc[(mat > 0).any(axis=1)]
    if mat.empty:
        return float("nan")
    mmean = mat.mean(axis=0)
    dmat = gene_by_lib_matrix(df, set(DRIVERS), tissue)
    rhos = []
    for d in DRIVERS:
        if d not in dmat.index:
            continue
        r, _ = stats.spearmanr(dmat.loc[d], mmean)
        if np.isfinite(r):
            rhos.append(r)
    return float(np.median(rhos)) if rhos else float("nan")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    defs = {k: load_def(k) for k in ("set", "module", "regulon")}
    glyc = load_geneset("GLYCOLYSIS")
    tables = load_species_tables()
    masses = pd.Series({sp: SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER})

    phen_rows = []
    coh_rows = []
    coup_rows = []

    for dname, genes in list(defs.items()) + [("GLYCOLYSIS", glyc)]:
        for tissue in ("fib", "epi"):
            y1s, y2s = {}, {}
            for sp in SPECIES_ORDER:
                df = tables[sp]
                m = y_metrics(df, genes, tissue)
                phen_rows.append(dict(definition=dname, tissue=tissue, species=sp, mass_g=masses[sp], **m))
                y1s[sp] = m["Y1"]
                y2s[sp] = m["Y2"]
                coh_rows.append(
                    dict(
                        definition=dname,
                        tissue=tissue,
                        species=sp,
                        mass_g=masses[sp],
                        coherence=coherence_species(df, genes, tissue),
                        n_libs=len(lib_cols(df, tissue)),
                    )
                )
                if dname != "GLYCOLYSIS":
                    coup_rows.append(
                        dict(
                            definition=dname,
                            tissue=tissue,
                            species=sp,
                            mass_g=masses[sp],
                            coupling=coupling_species(df, genes, tissue),
                        )
                    )
            # OLS
            for metric, series in [("Y1", y1s), ("Y2", y2s)]:
                s = pd.Series(series)
                fit = ols_logmass(s, masses)
                phen_rows.append(
                    dict(
                        definition=dname,
                        tissue=tissue,
                        species="OLS",
                        mass_g=np.nan,
                        n_genes=np.nan,
                        Y1=fit["slope"] if metric == "Y1" else np.nan,
                        Y2=fit["slope"] if metric == "Y2" else np.nan,
                        gene_sum=np.nan,
                        global_sum=np.nan,
                        metric=metric,
                        **{f"ols_{k}": v for k, v in fit.items()},
                    )
                )

    phen = pd.DataFrame(phen_rows)
    phen.to_csv(OUT / "phenotype_by_species.tsv", sep="\t", index=False)

    # clean OLS table
    ols_rows = []
    for dname in list(defs) + ["GLYCOLYSIS"]:
        for tissue in ("fib", "epi"):
            for metric in ("Y1", "Y2"):
                sub = phen[(phen.definition == dname) & (phen.tissue == tissue) & (phen.species != "OLS")]
                s = sub.set_index("species")[metric]
                fit = ols_logmass(s.reindex(SPECIES_ORDER), masses)
                ols_rows.append(dict(definition=dname, tissue=tissue, metric=metric, **fit))
    # Y3 salmon fib only
    for dname, genes in defs.items():
        y3 = salmon_y3(genes)
        fit = ols_logmass(y3.reindex(SPECIES_ORDER), masses)
        ols_rows.append(dict(definition=dname, tissue="fib", metric="Y3_salmon_mean_logTPM", **fit))
        for sp in SPECIES_ORDER:
            phen_rows.append(
                dict(
                    definition=dname,
                    tissue="fib",
                    species=sp,
                    mass_g=masses[sp],
                    n_genes=np.nan,
                    Y1=np.nan,
                    Y2=np.nan,
                    gene_sum=np.nan,
                    global_sum=np.nan,
                    Y3=y3.get(sp, np.nan),
                )
            )
    ols = pd.DataFrame(ols_rows)
    ols.to_csv(OUT / "ols_vs_log10mass.tsv", sep="\t", index=False)

    # matched High 1F1M for fib module (primary) + set/regulon
    match_rows = []
    for dname, genes in defs.items():
        y1m = {}
        for sp in SPECIES_ORDER:
            df = tables[sp]
            cols = lib_cols(df, "fib")
            mcols = filter_matched_cols(cols, MATCHED_HIGH_CODES.get(sp, []))
            if not mcols:
                # fallback: all fib
                mcols = cols
                used = "all_fib_fallback"
            else:
                used = "matched_high"
            m = y_metrics(df, genes, "fib", cols=mcols)
            match_rows.append(dict(definition=dname, species=sp, mass_g=masses[sp], cols=",".join(mcols), mode=used, **m))
            y1m[sp] = m["Y1"]
        fit = ols_logmass(pd.Series(y1m).reindex(SPECIES_ORDER), masses)
        match_rows.append(dict(definition=dname, species="OLS", mass_g=np.nan, cols="", mode="ols", n_genes=np.nan, Y1=fit["slope"], Y2=np.nan, gene_sum=np.nan, global_sum=np.nan, **{f"ols_{k}": v for k, v in fit.items()}))
    match = pd.DataFrame(match_rows)
    match.to_csv(OUT / "matched_high_Y1.tsv", sep="\t", index=False)

    # residual Y1 ~ mass | global mean (fib, each def)
    resid_rows = []
    for dname, genes in defs.items():
        y1, gmean = {}, {}
        for sp in SPECIES_ORDER:
            df = tables[sp]
            m = y_metrics(df, genes, "fib")
            y1[sp] = m["Y1"]
            gmean[sp] = global_promoter_mean(df, "fib")
        y1s = pd.Series(y1)
        gms = pd.Series(gmean)
        res = residual_after_covariate(y1s, gms)
        fit_raw = ols_logmass(y1s.reindex(SPECIES_ORDER), masses)
        fit_res = ols_logmass(res.reindex(SPECIES_ORDER), masses)
        resid_rows.append(dict(definition=dname, metric="Y1_raw", **fit_raw))
        resid_rows.append(dict(definition=dname, metric="Y1_resid_global", **fit_res))
        for sp in SPECIES_ORDER:
            resid_rows.append(dict(definition=dname, metric="per_species", species=sp, Y1=y1s[sp], global_mean=gms[sp], Y1_resid=res[sp], mass_g=masses[sp]))
    pd.DataFrame(resid_rows).to_csv(OUT / "residual_Y1_after_global.tsv", sep="\t", index=False)

    coh = pd.DataFrame(coh_rows)
    coh.to_csv(OUT / "coherence_by_species.tsv", sep="\t", index=False)
    coup = pd.DataFrame(coup_rows)
    coup.to_csv(OUT / "coupling_by_species.tsv", sep="\t", index=False)

    # OLS for coherence/coupling
    net_ols = []
    for dname in defs:
        for tissue in ("fib", "epi"):
            for metric, tab, col in [
                ("coherence", coh, "coherence"),
                ("coupling", coup, "coupling"),
            ]:
                sub = tab[(tab.definition == dname) & (tab.tissue == tissue)]
                if sub.empty:
                    continue
                s = sub.set_index("species")[col].reindex(SPECIES_ORDER)
                fit = ols_logmass(s, masses)
                net_ols.append(dict(definition=dname, tissue=tissue, metric=metric, **fit))
    pd.DataFrame(net_ols).to_csv(OUT / "network_scalars_ols.tsv", sep="\t", index=False)

    # matched OLS summary lines
    match_ols = []
    for dname in defs:
        sub = match[(match.definition == dname) & (match.species != "OLS")]
        fit = ols_logmass(sub.set_index("species")["Y1"].reindex(SPECIES_ORDER), masses)
        y3 = salmon_y3(defs[dname])
        fit3 = ols_logmass(y3.reindex(SPECIES_ORDER), masses)
        match_ols.append(dict(definition=dname, metric="Y1_matched_or_fallback", **fit))
        match_ols.append(dict(definition=dname, metric="Y3_salmon", **fit3))
    pd.DataFrame(match_ols).to_csv(OUT / "matched_vs_salmon_ols.tsv", sep="\t", index=False)

    # Figure: Y1 fib/epi for three defs + Y3
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), constrained_layout=True)
    x = np.log10(masses.reindex(SPECIES_ORDER).values)
    for j, dname in enumerate(("set", "module", "regulon")):
        ax = axes[0, j]
        sub = phen[(phen.definition == dname) & (phen.tissue == "fib") & (phen.species != "OLS")]
        y = sub.set_index("species")["Y1"].reindex(SPECIES_ORDER)
        ax.scatter(x, y.values, c="#e76f51", s=40)
        fit = ols_logmass(y, masses.reindex(SPECIES_ORDER))
        if np.isfinite(fit["slope"]):
            xx = np.linspace(x.min(), x.max(), 50)
            ax.plot(xx, fit["intercept"] + fit["slope"] * xx, color="#e76f51", lw=1)
        ax.set_title(f"fib Y1 {dname}\np={fit['p']:.3g}")
        ax.set_xlabel("log10 mass")
        if j == 0:
            ax.set_ylabel("mean gene TPM")
        ax = axes[1, j]
        sube = phen[(phen.definition == dname) & (phen.tissue == "epi") & (phen.species != "OLS")]
        ye = sube.set_index("species")["Y1"].reindex(SPECIES_ORDER)
        ax.scatter(x, ye.values, c="#2a9d8f", s=40)
        fite = ols_logmass(ye, masses.reindex(SPECIES_ORDER))
        if np.isfinite(fite["slope"]):
            xx = np.linspace(x.min(), x.max(), 50)
            ax.plot(xx, fite["intercept"] + fite["slope"] * xx, color="#2a9d8f", lw=1)
        ax.set_title(f"epi Y1 {dname}\np={fite['p']:.3g}")
        ax.set_xlabel("log10 mass")
        if j == 0:
            ax.set_ylabel("mean gene TPM")
    fig.suptitle("OXPHOS set / module / regulon-lite — CAGE mean TPM vs body mass", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_thread_Y1_fib_epi_three_defs.png", dpi=160)
    fig.savefig(OUT / "fig_oxphos_thread_Y1_fib_epi_three_defs.png", dpi=160)
    plt.close()

    # salmon vs matched for module
    fig, ax = plt.subplots(figsize=(5.2, 4.2), constrained_layout=True)
    y1m = match[(match.definition == "module") & (match.species != "OLS")].set_index("species")["Y1"].reindex(SPECIES_ORDER)
    y3 = salmon_y3(defs["module"]).reindex(SPECIES_ORDER)
    ax.scatter(np.log10(masses.values), y1m.values, label="CAGE Y1 matched/fallback", c="#e76f51")
    ax2 = ax.twinx()
    ax2.scatter(np.log10(masses.values), y3.values, label="Salmon Y3", c="#457b9d", marker="s")
    for sp, xi in zip(SPECIES_ORDER, np.log10(masses.values)):
        ax.text(xi, y1m[sp], sp[:3], fontsize=7)
    f1 = ols_logmass(y1m, masses)
    f3 = ols_logmass(y3, masses)
    ax.set_title(f"module: CAGE p={f1['p']:.3g} / Salmon p={f3['p']:.3g}")
    ax.set_xlabel("log10 mass")
    ax.set_ylabel("CAGE Y1")
    ax2.set_ylabel("Salmon mean log10(TPM+1)")
    fig.savefig(FIG / "fig_oxphos_thread_module_cage_vs_salmon.png", dpi=160)
    fig.savefig(OUT / "fig_oxphos_thread_module_cage_vs_salmon.png", dpi=160)
    plt.close()

    summary = {
        "ols": ols.to_dict(orient="records"),
        "matched_vs_salmon": match_ols,
        "network_ols": net_ols,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(ols.to_string(index=False))
    print("\nMatched vs Salmon:")
    print(pd.DataFrame(match_ols).to_string(index=False))


if __name__ == "__main__":
    main()
