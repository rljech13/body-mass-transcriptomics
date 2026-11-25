#!/usr/bin/env python3
"""Re-verify the OXPHOS-vs-mass decline against explanations the corpus did not test.

Reads only phenotype tables already produced by scripts 21 (GEO RNA) and 22 (CAGE
GEO-style). Nothing is recomputed from raw data, so any disagreement with the
corpus is a disagreement about interpretation of the same numbers.

Explanations tested:
  E1 gene-count / annotation coverage  (n_detected, variable vs fixed 7/7 core)
  E2 denominator inflation             (lib_total, mapped_mean vs mass)
  E3 single-species leverage           (leave-one-out)
  E4 phylogeny                         (PGLS under BM + Pagel lambda; Galliformes-only)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
RES = ROOT / "results" / "oxphos_species_thread"
OUT = RES / "reverify_confounds"

SPECIES_ORDER = ["Zebra", "Quail", "Chicken", "Duck", "Peacock", "Turkey", "Emu"]
MASS = {"Zebra": 11.93, "Quail": 96.6, "Chicken": 912.5, "Duck": 1121.0,
        "Peacock": 4093.75, "Turkey": 5811.0, "Emu": 36200.0}
# annotation provenance of the expression tables (Augustus = no official gene names)
AUGUSTUS = {"Peacock", "Turkey", "Emu"}
GALLIFORMES = ["Quail", "Chicken", "Peacock", "Turkey"]

# approximate crown divergence times (Mya) for a BM/lambda sensitivity check
TREE = {
    "root": 110.0,       # Palaeognathae (Emu) | Neognathae
    "neognath": 100.0,   # Neoaves (Zebra) | Galloanserae
    "galloanser": 90.0,  # Anseriformes (Duck) | Galliformes
    "phasianid": 37.0,   # Quail | rest
    "pav_gal": 32.0,     # Peacock | (Chicken, Turkey)
    "gal_mel": 27.0,     # Chicken | Turkey
}


def shared_age(a: str, b: str) -> float:
    """Time from the root down to the MRCA of a and b (shared branch length)."""
    if a == b:
        return TREE["root"]
    s = {a, b}
    if "Emu" in s:
        return 0.0
    if "Zebra" in s:
        return TREE["root"] - TREE["neognath"]
    if "Duck" in s:
        return TREE["root"] - TREE["galloanser"]
    if "Quail" in s:
        return TREE["root"] - TREE["phasianid"]
    if "Peacock" in s:
        return TREE["root"] - TREE["pav_gal"]
    return TREE["root"] - TREE["gal_mel"]


def vcv() -> pd.DataFrame:
    m = pd.DataFrame(0.0, index=SPECIES_ORDER, columns=SPECIES_ORDER)
    for a in SPECIES_ORDER:
        for b in SPECIES_ORDER:
            m.loc[a, b] = shared_age(a, b)
    return m / TREE["root"]


def ols(y: pd.Series, x: pd.Series) -> dict:
    yy, xx = y.astype(float).values, x.astype(float).values
    ok = np.isfinite(yy) & np.isfinite(xx)
    if ok.sum() < 3:
        return dict(n=int(ok.sum()), slope=np.nan, r=np.nan, p=np.nan)
    sl, ic, r, p, se = stats.linregress(xx[ok], yy[ok])
    return dict(n=int(ok.sum()), slope=float(sl), intercept=float(ic), r=float(r), p=float(p), se=float(se))


def partial_corr(y: pd.Series, x: pd.Series, z: pd.Series) -> dict:
    """Correlation of y and x after removing linear effect of z, with exact df."""
    d = pd.concat([y, x, z], axis=1).dropna()
    d.columns = ["y", "x", "z"]
    n = len(d)
    if n < 4:
        return dict(n=n, r=np.nan, p=np.nan)
    ry = d["y"] - np.polyval(np.polyfit(d["z"], d["y"], 1), d["z"])
    rx = d["x"] - np.polyval(np.polyfit(d["z"], d["x"], 1), d["z"])
    r = float(np.corrcoef(ry, rx)[0, 1])
    df = n - 3
    if df <= 0 or abs(r) >= 1:
        return dict(n=n, r=r, p=np.nan, df=df)
    t = r * np.sqrt(df / (1 - r * r))
    return dict(n=n, r=r, p=float(2 * stats.t.sf(abs(t), df)), df=df)


def gls_fit(y: np.ndarray, x: np.ndarray, C: np.ndarray) -> dict:
    X = np.column_stack([np.ones_like(x), x])
    Ci = np.linalg.inv(C)
    XtCiX = X.T @ Ci @ X
    beta = np.linalg.solve(XtCiX, X.T @ Ci @ y)
    resid = y - X @ beta
    n, k = len(y), X.shape[1]
    s2 = float(resid.T @ Ci @ resid) / (n - k)
    cov = s2 * np.linalg.inv(XtCiX)
    se = float(np.sqrt(cov[1, 1]))
    t = float(beta[1] / se)
    ll = -0.5 * (n * np.log(2 * np.pi * s2) + np.log(np.linalg.det(C)) + float(resid.T @ Ci @ resid) / s2)
    return dict(slope=float(beta[1]), se=se, t=t, p=float(2 * stats.t.sf(abs(t), n - k)), loglik=ll)


def pgls(y: pd.Series, x: pd.Series, C0: pd.DataFrame) -> dict:
    d = pd.concat([y, x], axis=1).dropna()
    sp = [s for s in SPECIES_ORDER if s in d.index]
    yv = d.loc[sp].iloc[:, 0].values.astype(float)
    xv = d.loc[sp].iloc[:, 1].values.astype(float)
    C = C0.loc[sp, sp].values
    bm = gls_fit(yv, xv, C)

    def neg_ll(lam: float) -> float:
        Cl = C.copy()
        off = ~np.eye(len(sp), dtype=bool)
        Cl[off] *= lam
        try:
            return -gls_fit(yv, xv, Cl)["loglik"]
        except np.linalg.LinAlgError:
            return 1e9

    res = minimize_scalar(neg_ll, bounds=(0.0, 1.0), method="bounded")
    lam = float(res.x)
    Cl = C.copy()
    Cl[~np.eye(len(sp), dtype=bool)] *= lam
    lm = gls_fit(yv, xv, Cl)
    return {"n": len(sp), "bm_slope": bm["slope"], "bm_p": bm["p"],
            "lambda_ml": lam, "lambda_slope": lm["slope"], "lambda_p": lm["p"]}


def loo(y: pd.Series, x: pd.Series) -> pd.DataFrame:
    rows = [dict(dropped="(none)", **{k: v for k, v in ols(y, x).items() if k in ("n", "slope", "r", "p")})]
    for sp in y.dropna().index:
        f = ols(y.drop(sp), x.drop(sp))
        rows.append(dict(dropped=sp, **{k: v for k, v in f.items() if k in ("n", "slope", "r", "p")}))
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mass = pd.Series(MASS).reindex(SPECIES_ORDER)
    logm = np.log10(mass)
    C0 = vcv()

    cage = pd.read_csv(RES / "cage_geo_style_recompute" / "phenotype_by_species.tsv", sep="\t")
    geo = pd.read_csv(RES / "rnaseq_geo_recompute" / "phenotype_by_species.tsv", sep="\t")

    panels = {}
    for tiss in ("fib", "epi"):
        for d in ("set", "module", "GLYCOLYSIS"):
            sub = cage[(cage.tissue == tiss) & (cage.definition == d)].set_index("species").reindex(SPECIES_ORDER)
            panels[(f"CAGE_{tiss}", d)] = sub
    for d in ("set", "module", "GLYCOLYSIS"):
        sub = geo[geo.definition == d].set_index("species").reindex(SPECIES_ORDER)
        panels[("GEO_RNA", d)] = sub

    # ---------- E1/E2: what actually tracks mass besides the phenotype ----------
    cov_rows = []
    for (layer, d), sub in panels.items():
        for col in ("n_detected", "lib_total", "mapped_mean", "n_core"):
            if col not in sub.columns:
                continue
            f = ols(sub[col], logm)
            s = stats.spearmanr(logm, sub[col], nan_policy="omit")
            cov_rows.append(dict(layer=layer, definition=d, covariate=col,
                                 r=f["r"], p=f["p"], rho=float(s.statistic), rho_p=float(s.pvalue)))
    cov = pd.DataFrame(cov_rows).drop_duplicates(subset=["layer", "covariate", "definition"])
    cov.to_csv(OUT / "covariates_vs_mass.tsv", sep="\t", index=False)

    # ---------- main comparison: variable-composition vs fixed 7/7 core ----------
    main_rows = []
    for (layer, d), sub in panels.items():
        for metric in ("frac_detected", "frac_soft", "frac_core", "Y1_detected", "Y1_core",
                       "mean_linear_detected", "mean_linear_core", "mean_log10_detected"):
            if metric not in sub.columns or sub[metric].notna().sum() < 4:
                continue
            f = ols(sub[metric], logm)
            s = stats.spearmanr(logm, sub[metric], nan_policy="omit")
            pg = pgls(sub[metric], logm, C0)
            pc = partial_corr(sub[metric], logm, sub["n_detected"]) if "n_detected" in sub else dict(r=np.nan, p=np.nan)
            gal = ols(sub[metric].reindex(GALLIFORMES), logm.reindex(GALLIFORMES))
            fixed = metric.endswith("_core")
            main_rows.append(dict(
                layer=layer, definition=d, metric=metric,
                gene_list="fixed_7of7" if fixed else "variable_per_species",
                ols_r=f["r"], ols_p=f["p"], rho=float(s.statistic), rho_p=float(s.pvalue),
                partial_r_given_ndetected=pc["r"], partial_p_given_ndetected=pc["p"],
                pgls_bm_p=pg["bm_p"], pgls_lambda_ml=pg["lambda_ml"], pgls_lambda_p=pg["lambda_p"],
                galliformes_n4_r=gal["r"], galliformes_n4_p=gal["p"],
            ))
    main = pd.DataFrame(main_rows)
    main.to_csv(OUT / "fixed_vs_variable_genelist.tsv", sep="\t", index=False)

    # ---------- E3: leave-one-out on the corpus headline metrics ----------
    headline = [("CAGE_fib", "set", "frac_detected"), ("CAGE_fib", "set", "Y1_detected"),
                ("CAGE_fib", "module", "Y1_detected"), ("GEO_RNA", "set", "frac_detected"),
                ("GEO_RNA", "module", "frac_detected"),
                ("CAGE_fib", "set", "frac_core"), ("GEO_RNA", "set", "frac_core")]
    loo_all = []
    for layer, d, metric in headline:
        sub = panels[(layer, d)]
        if metric not in sub.columns:
            continue
        t = loo(sub[metric], logm)
        t.insert(0, "metric", metric)
        t.insert(0, "definition", d)
        t.insert(0, "layer", layer)
        loo_all.append(t)
    loo_df = pd.concat(loo_all, ignore_index=True)
    loo_df.to_csv(OUT / "leave_one_out.tsv", sep="\t", index=False)

    # ---------- E1 quantified: numerator is a gene count ----------
    dec_rows = []
    for (layer, d), sub in panels.items():
        if "frac_core" not in sub.columns or "frac_detected" not in sub.columns:
            continue
        extra_genes = sub["n_detected"] - sub["n_core"]
        extra_frac = sub["frac_detected"] - sub["frac_core"]
        dec_rows.append(dict(
            layer=layer, definition=d,
            r_frac_vs_ndetected=ols(sub["frac_detected"], sub["n_detected"])["r"],
            p_frac_vs_ndetected=ols(sub["frac_detected"], sub["n_detected"])["p"],
            r_extrafrac_vs_extragenes=ols(extra_frac, extra_genes)["r"],
            p_extrafrac_vs_extragenes=ols(extra_frac, extra_genes)["p"],
            r_ndetected_vs_mass=ols(sub["n_detected"], logm)["r"],
            p_ndetected_vs_mass=ols(sub["n_detected"], logm)["p"],
        ))
    dec = pd.DataFrame(dec_rows)
    dec.to_csv(OUT / "numerator_is_gene_count.tsv", sep="\t", index=False)

    # ---------- annotation provenance vs mass ----------
    aug = pd.Series({sp: int(sp in AUGUSTUS) for sp in SPECIES_ORDER})
    prov = dict(
        r_augustus_vs_logmass=float(np.corrcoef(aug.values, logm.values)[0, 1]),
        p_augustus_vs_logmass=float(stats.pearsonr(aug.values, logm.values)[1]),
        rho_augustus_vs_logmass=float(stats.spearmanr(aug.values, logm.values).statistic),
        note="Augustus-annotated = Peacock, Turkey, Emu = the 3 heaviest species",
    )
    for (layer, d), sub in panels.items():
        if d != "set":
            continue
        prov[f"{layer}_ndetected_augustus_vs_named_t_p"] = float(
            stats.ttest_ind(sub.loc[list(AUGUSTUS), "n_detected"],
                            sub.loc[[s for s in SPECIES_ORDER if s not in AUGUSTUS], "n_detected"]).pvalue)
    (OUT / "annotation_provenance.json").write_text(json.dumps(prov, indent=2))

    # ---------- report ----------
    pd.set_option("display.width", 200, "display.max_columns", 40)
    print("=== covariates vs log10 mass ===")
    print(cov.to_string(index=False))
    print("\n=== fixed 7/7 core vs variable per-species gene list ===")
    print(main[["layer", "definition", "metric", "gene_list", "ols_r", "ols_p",
                "partial_p_given_ndetected", "pgls_lambda_p", "galliformes_n4_p"]].to_string(index=False))
    print("\n=== numerator decomposition ===")
    print(dec.to_string(index=False))
    print("\n=== leave-one-out ===")
    print(loo_df.to_string(index=False))
    print("\n=== annotation provenance ===")
    print(json.dumps(prov, indent=2))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
