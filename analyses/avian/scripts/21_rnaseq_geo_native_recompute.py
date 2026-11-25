#!/usr/bin/env python3
"""Recompute fibroblast RNA-seq OXPHOS ~ mass on official GSE213253 tables.

Why this exists: phase-2 Y3 used salmon→chicken mean log10(TPM+1) and called RNA
null. That was the wrong primary layer. GEO native RNA-seq is the published
quantification; Peacock/Turkey/Emu Augustus IDs map through Ruslan DE
transcriptId. GEO 'norm' library totals are not a shared TPM scale, so the
primary metric is within-table transcriptome fraction (OXPHOS / all genes).
Salmon remains a sensitivity analysis, not the headline.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    FIG,
    RES,
    SALMON,
    SPECIES,
    SPECIES_ORDER,
    load_geneset,
    load_geo_rnaseq,
    nuclear_oxphos,
    ols_logmass,
    spearman_logmass,
)

OUT = RES / "rnaseq_geo_recompute"
SALIB = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/results/salmon_all7_chicken/gene_tpm_by_library.tsv")


def load_defs() -> dict[str, set[str]]:
    mod = set(pd.read_csv(RES / "definitions" / "oxphos_module.tsv", sep="\t")["gene"].str.upper())
    return {
        "set": nuclear_oxphos(),
        "module": {g for g in mod if g},
        "GLYCOLYSIS": load_geneset("GLYCOLYSIS"),
    }


def fit_pair(y: pd.Series, masses: pd.Series) -> dict:
    ols = ols_logmass(y, masses)
    spr = spearman_logmass(y, masses)
    return {**{f"ols_{k}": v for k, v in ols.items()}, "spearman_rho": spr["rho"], "spearman_p": spr["p"]}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    defs = load_defs()
    masses = pd.Series({sp: SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER})

    geo = {sp: load_geo_rnaseq(sp) for sp in SPECIES_ORDER}
    qc = pd.DataFrame(
        [
            dict(
                species=sp,
                mass_g=masses[sp],
                mapped_how=geo[sp]["mapped_how"],
                n_raw=geo[sp]["n_raw"],
                n_mapped_rows=geo[sp]["n_mapped_rows"],
                n_symbols=geo[sp]["n_symbols"],
                mean_lib_total=float(geo[sp]["lib_total"].mean()),
            )
            for sp in SPECIES_ORDER
        ]
    )
    qc.to_csv(OUT / "geo_mapping_qc.tsv", sep="\t", index=False)

    present = {sp: set(geo[sp]["gene_mat"].index) for sp in SPECIES_ORDER}
    cores = {d: set.intersection(*[present[sp] & genes for sp in SPECIES_ORDER]) for d, genes in defs.items()}

    cov_rows = []
    for dname, genes in defs.items():
        for sp in SPECIES_ORDER:
            cov_rows.append(
                dict(
                    definition=dname,
                    species=sp,
                    n_list=len(genes),
                    n_detected=len(present[sp] & genes),
                    n_core7=len(cores[dname]),
                )
            )
    pd.DataFrame(cov_rows).to_csv(OUT / "gene_coverage.tsv", sep="\t", index=False)

    phen_rows = []
    lib_rows = []
    for dname, genes in defs.items():
        core = cores[dname]
        for sp in SPECIES_ORDER:
            gmat = geo[sp]["gene_mat"]
            tot = geo[sp]["lib_total"]
            detected = sorted(present[sp] & genes)
            core_g = sorted(core)
            gmean = gmat.mean(axis=1)
            num_det = float(gmean.reindex(detected).fillna(0).sum())
            num_core = float(gmean.reindex(core_g).fillna(0).sum()) if core_g else np.nan
            den = float(tot.mean())
            mapped_mean = float(gmean.mean())
            phen_rows.append(
                dict(
                    definition=dname,
                    species=sp,
                    mass_g=masses[sp],
                    n_detected=len(detected),
                    n_core=len(core_g),
                    frac_detected=num_det / den if den else np.nan,
                    frac_core=num_core / den if den and core_g else np.nan,
                    mean_linear_detected=float(gmean.reindex(detected).mean()) if detected else np.nan,
                    mean_log10_detected=float(np.log10(gmean.reindex(detected).clip(lower=0) + 1).mean()) if detected else np.nan,
                    mean_linear_core=float(gmean.reindex(core_g).mean()) if core_g else np.nan,
                    rel_core_over_mapped=(float(gmean.reindex(core_g).mean()) / mapped_mean) if core_g and mapped_mean else np.nan,
                    lib_total=den,
                    mapped_mean=mapped_mean,
                )
            )
            for col in gmat.columns:
                ser = gmat[col]
                lib_den = float(tot[col])
                lib_rows.append(
                    dict(
                        definition=dname,
                        species=sp,
                        library=col,
                        mass_g=masses[sp],
                        frac_detected=float(ser.reindex(detected).fillna(0).sum()) / lib_den if lib_den else np.nan,
                        frac_core=float(ser.reindex(core_g).fillna(0).sum()) / lib_den if lib_den and core_g else np.nan,
                        mean_log10_detected=float(np.log10(ser.reindex(detected).clip(lower=0) + 1).mean()) if detected else np.nan,
                    )
                )

    phen = pd.DataFrame(phen_rows)
    phen.to_csv(OUT / "phenotype_by_species.tsv", sep="\t", index=False)
    lib = pd.DataFrame(lib_rows)
    lib.to_csv(OUT / "phenotype_by_library.tsv", sep="\t", index=False)

    metrics = [
        "frac_detected",
        "frac_core",
        "mean_linear_detected",
        "mean_log10_detected",
        "mean_linear_core",
        "rel_core_over_mapped",
    ]
    ols_rows = []
    for dname in defs:
        sub = phen[phen.definition == dname].set_index("species").reindex(SPECIES_ORDER)
        for metric in metrics:
            y = sub[metric]
            fit = fit_pair(y, masses)
            ols_rows.append(dict(source="GEO", definition=dname, unit="species_mean_n7", metric=metric, **fit))
        # n=4 symbol-native only
        four = ["Zebra", "Quail", "Chicken", "Duck"]
        sub4 = sub.loc[four]
        m4 = masses.loc[four]
        for metric in ("frac_detected", "mean_log10_detected"):
            fit = fit_pair(sub4[metric], m4)
            ols_rows.append(dict(source="GEO", definition=dname, unit="species_mean_n4_symbols", metric=metric, **fit))

    # per-library n=14
    for dname in defs:
        sub = lib[lib.definition == dname]
        xmass = pd.Series(sub["mass_g"].values)
        for metric in ("frac_detected", "frac_core", "mean_log10_detected"):
            y = pd.Series(sub[metric].values)
            # fake index; ols_logmass uses values only
            y.index = range(len(y))
            xmass.index = y.index
            fit = fit_pair(y, xmass)
            ols_rows.append(dict(source="GEO", definition=dname, unit="library_n14", metric=metric, **fit))

    # Salmon sensitivity
    sal = pd.read_csv(SALMON, sep="\t", index_col=0)
    sal.index = sal.index.astype(str).str.upper()
    salib = pd.read_csv(SALIB, sep="\t", index_col=0)
    salib.index = salib.index.astype(str).str.upper()
    for dname, genes in defs.items():
        present_g = [g for g in genes if g in sal.index]
        tot = sal.sum(axis=0)
        num = sal.loc[present_g].sum(axis=0)
        y_frac, y_log = {}, {}
        for sp in SPECIES_ORDER:
            col = SPECIES[sp]["salmon"]
            y_frac[sp] = float(num[col] / tot[col])
            vals = sal.loc[present_g, col].astype(float).clip(lower=0)
            y_log[sp] = float(np.log10(vals + 1).mean())
        for metric, series in [("frac_detected", y_frac), ("mean_log10_detected", y_log)]:
            y = pd.Series(series).reindex(SPECIES_ORDER)
            fit = fit_pair(y, masses)
            ols_rows.append(dict(source="Salmon_chicken_proxy", definition=dname, unit="species_mean_n7", metric=metric, **fit))
        # per-lib salmon
        xs, yf, yl = [], [], []
        for col in salib.columns:
            sp = next((s for s in SPECIES_ORDER if col.startswith(SPECIES[s]["salmon"] + "_")), None)
            if sp is None:
                continue
            present_g = [g for g in genes if g in salib.index]
            v = salib.loc[present_g, col].astype(float).clip(lower=0)
            xs.append(masses[sp])
            yf.append(float(v.sum() / salib[col].sum()))
            yl.append(float(np.log10(v + 1).mean()))
        for metric, ys in [("frac_detected", yf), ("mean_log10_detected", yl)]:
            y = pd.Series(ys)
            xm = pd.Series(xs)
            y.index = xm.index = range(len(ys))
            fit = fit_pair(y, xm)
            ols_rows.append(dict(source="Salmon_chicken_proxy", definition=dname, unit="library_n14", metric=metric, **fit))

    ols = pd.DataFrame(ols_rows)
    ols.to_csv(OUT / "ols_spearman_vs_log10mass.tsv", sep="\t", index=False)

    # Figures: GEO fraction (primary) vs Salmon log (old Y3)
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.8), constrained_layout=True)
    x = np.log10(masses.reindex(SPECIES_ORDER).values)
    for j, dname in enumerate(("set", "module", "GLYCOLYSIS")):
        sub = phen[phen.definition == dname].set_index("species").reindex(SPECIES_ORDER)
        ax = axes[0, j]
        y = sub["frac_detected"]
        ax.scatter(x, y.values, c="#1d3557", s=50)
        fit = ols_logmass(y, masses)
        if np.isfinite(fit["slope"]):
            xx = np.linspace(x.min(), x.max(), 40)
            ax.plot(xx, fit["intercept"] + fit["slope"] * xx, color="#1d3557", lw=1.4)
        for sp, xi, yi in zip(SPECIES_ORDER, x, y.values):
            ax.text(xi, yi, sp[:3], fontsize=7)
        ax.set_title(f"GEO frac detected {dname}\nOLS p={fit['p']:.3g}")
        ax.set_xlabel("log10 mass (g)")
        if j == 0:
            ax.set_ylabel("OXPHOS / library total")
        ax = axes[1, j]
        y2 = sub["frac_core"]
        ax.scatter(x, y2.values, c="#457b9d", s=50)
        fit2 = ols_logmass(y2, masses)
        if np.isfinite(fit2["slope"]):
            xx = np.linspace(x.min(), x.max(), 40)
            ax.plot(xx, fit2["intercept"] + fit2["slope"] * xx, color="#457b9d", lw=1.4)
        ax.set_title(f"GEO frac common-core {dname}\nOLS p={fit2['p']:.3g}")
        ax.set_xlabel("log10 mass (g)")
        if j == 0:
            ax.set_ylabel("core / library total")
    fig.suptitle("GEO native RNA-seq: OXPHOS transcriptome fraction vs body mass", fontsize=11)
    fig.savefig(FIG / "fig_rnaseq_geo_fraction_vs_mass.png", dpi=160)
    fig.savefig(OUT / "fig_rnaseq_geo_fraction_vs_mass.png", dpi=160)
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), constrained_layout=True)
    sub = phen[phen.definition == "module"].set_index("species").reindex(SPECIES_ORDER)
    ax = axes[0]
    y = sub["frac_detected"]
    ax.scatter(x, y.values, c="#1d3557", s=55)
    fit = ols_logmass(y, masses)
    xx = np.linspace(x.min(), x.max(), 40)
    ax.plot(xx, fit["intercept"] + fit["slope"] * xx, color="#1d3557", lw=1.5)
    for sp, xi, yi in zip(SPECIES_ORDER, x, y.values):
        ax.text(xi, yi, sp[:3], fontsize=8)
    ax.set_title(f"GEO module fraction\nOLS r={fit['r']:.2f} p={fit['p']:.3g}")
    ax.set_xlabel("log10 mass (g)")
    ax.set_ylabel("module / library total")
    ax = axes[1]
    y3 = pd.Series(
        {
            sp: float(np.log10(sal.loc[[g for g in defs["module"] if g in sal.index], SPECIES[sp]["salmon"]].clip(lower=0) + 1).mean())
            for sp in SPECIES_ORDER
        }
    ).reindex(SPECIES_ORDER)
    ax.scatter(x, y3.values, c="#e76f51", s=55, marker="s")
    fit3 = ols_logmass(y3, masses)
    ax.plot(xx, fit3["intercept"] + fit3["slope"] * xx, color="#e76f51", lw=1.5)
    for sp, xi, yi in zip(SPECIES_ORDER, x, y3.values):
        ax.text(xi, yi, sp[:3], fontsize=8)
    ax.set_title(f"old Y3 Salmon mean log10(TPM+1)\nOLS r={fit3['r']:.2f} p={fit3['p']:.3g}")
    ax.set_xlabel("log10 mass (g)")
    ax.set_ylabel("mean log10(TPM+1)")
    fig.suptitle("Why the old RNA call was wrong: GEO fraction vs Salmon-to-chicken log-mean", fontsize=11)
    fig.savefig(FIG / "fig_rnaseq_geo_frac_vs_old_salmon_Y3.png", dpi=160)
    fig.savefig(OUT / "fig_rnaseq_geo_frac_vs_old_salmon_Y3.png", dpi=160)
    plt.close()

    primary = ols[
        (ols.source == "GEO")
        & (ols.unit == "species_mean_n7")
        & (ols.metric.isin(["frac_detected", "frac_core", "mean_log10_detected"]))
    ]
    summary = {
        "date": "2026-08-14",
        "primary_layer": "GSE213253 official RNAseq.xlsx, Augustus→symbol via Ruslan DE transcriptId",
        "primary_metric": "frac_detected = sum(detected OXPHOS gene means) / mean(library total)",
        "why_not_absolute_TPM": "GEO norm library totals differ ~400k (zebra) vs ~1.0M (peacock); not a shared TPM scale",
        "old_Y3_superseded": "salmon→chicken mean log10(TPM+1) remains null; that is a proxy, not the official RNA layer",
        "cores": {k: sorted(v) for k, v in cores.items()},
        "mapping_qc": qc.to_dict(orient="records"),
        "primary_fits": primary.to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(qc.to_string(index=False))
    print()
    show = ols[(ols.unit.isin(["species_mean_n7", "library_n14"])) & (ols.metric.isin(["frac_detected", "frac_core", "mean_log10_detected"]))]
    print(show.to_string(index=False))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
