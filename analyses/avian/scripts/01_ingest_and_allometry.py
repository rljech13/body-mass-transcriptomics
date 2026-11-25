#!/usr/bin/env python3
"""A1–A3: ingest GSE213253, OXPHOS module allometry, exploratory gene-subset search."""
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
RES = ROOT / "results"
FIG = ROOT / "figures"
for d in (RES, FIG, DATA / "genesets"):
    d.mkdir(parents=True, exist_ok=True)

# Adult body mass (kg). Sources: Dunning CRC Handbook / standard avian references.
# Documented for transparency; not from GSE213253 metadata.
SPECIES_META = {
    "ZebraFinch": {
        "mass_kg": 0.0157,
        "latin": "Taeniopygia guttata",
        "rnaseq": "GSE213253_Zebrafinch_RNAseq.xlsx",
        "cage": "GSE213253_ZebraFinch_CAGE.xlsx",
        "symbol_ok": True,
    },
    "Quail": {
        "mass_kg": 0.12,
        "latin": "Coturnix japonica",
        "rnaseq": "GSE213253_Quail_RNAseq.xlsx",
        "cage": "GSE213253_Quail_CAGE.xlsx",
        "symbol_ok": True,
    },
    "Chicken": {
        "mass_kg": 2.0,
        "latin": "Gallus gallus",
        "rnaseq": "GSE213253_Chicken_RNAseq.xlsx",
        "cage": "GSE213253_Chicken_CAGE.xlsx",
        "symbol_ok": True,
    },
    "Duck": {
        "mass_kg": 1.1,
        "latin": "Anas platyrhynchos",
        "rnaseq": "GSE213253_Duck_RNAseq.xlsx",
        "cage": "GSE213253_Duck_CAGE.xlsx",
        "symbol_ok": True,
    },
    "Peacock": {
        "mass_kg": 4.5,
        "latin": "Pavo cristatus",
        "rnaseq": "GSE213253_Peacock_RNAseq.xlsx",
        "cage": "GSE213253_Peacock_CAGE.xlsx",
        "symbol_ok": False,
    },
    "Turkey": {
        "mass_kg": 8.0,
        "latin": "Meleagris gallopavo",
        "rnaseq": "GSE213253_Turkey_RNAseq.xlsx",
        "cage": "GSE213253_Turkey_CAGE.xlsx",
        "symbol_ok": False,
    },
    "Emu": {
        "mass_kg": 35.0,
        "latin": "Dromaius novaehollandiae",
        "rnaseq": "GSE213253_Emu_RNAseq.xlsx",
        "cage": "GSE213253_Emu_CAGE.xlsx",
        "symbol_ok": False,
    },
}

# Core OXPHOS / mito / glycolysis symbols (human/chicken-shared naming).
# Pre-registered; not optimized against mass.
OXPHOS_CORE = sorted(
    {
        "NDUFA1", "NDUFA2", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA8", "NDUFA9", "NDUFA10",
        "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10",
        "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8",
        "NDUFV1", "NDUFV2", "NDUFV3",
        "SDHA", "SDHB", "SDHC", "SDHD",
        "UQCRC1", "UQCRC2", "UQCRB", "UQCRQ", "UQCRFS1", "CYC1", "CYCS",
        "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7B", "COX7C", "COX8A",
        "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO",
        "ATP5MC1", "ATP5MC3", "ATP5ME", "ATP5MF", "ATP5MG",
        # older chicken aliases often present in GEO tables
        "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E", "ATP5G1", "ATP5G3", "ATP5H", "ATP5J", "ATP5J2", "ATP5L", "ATP5O",
        "COX4", "COX7A2L", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8",
    }
)

MITO_BIOGENESIS = sorted(
    {
        "TFAM", "TFB1M", "TFB2M", "POLG", "POLRMT", "SSBP1", "TWNK",
        "PPARGC1A", "PPARGC1B", "NRF1", "GABPA", "ESRRA", "ESRRB",
        "PPRC1", "SIRT1", "SIRT3", "LONP1", "CLPP",
    }
)

GLYCOLYSIS = sorted(
    {
        "HK1", "HK2", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA", "ALDOB", "ALDOC",
        "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "PKLR", "LDHA", "LDHB",
        "SLC2A1", "SLC2A3", "PFKFB3",
    }
)

TCA = sorted(
    {
        "CS", "ACO2", "IDH2", "IDH3A", "IDH3B", "IDH3G", "OGDH", "DLST", "DLD",
        "SUCLA2", "SUCLG1", "SDHA", "SDHB", "FH", "MDH2", "GOT2", "MPC1", "MPC2",
    }
)

MODULES = {
    "OXPHOS": OXPHOS_CORE,
    "MITO_BIOGENESIS": MITO_BIOGENESIS,
    "GLYCOLYSIS": GLYCOLYSIS,
    "TCA": TCA,
}


def load_rnaseq_gene_mean(path: Path) -> pd.Series:
    df = pd.read_excel(path)
    norms = [c for c in df.columns if str(c).startswith("norm.")]
    if not norms:
        raise ValueError(f"No norm columns in {path}")
    g = df["gene"].astype(str).str.upper()
    mat = df[norms].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # mean across available sexes / replicates
    s = mat.mean(axis=1)
    out = pd.Series(s.values, index=g.values, name=path.stem)
    # collapse duplicate gene symbols by mean
    out = out.groupby(level=0).mean()
    return out


def load_cage_gene_sum(path: Path) -> pd.Series:
    df = pd.read_excel(path)
    norms = [c for c in df.columns if str(c).startswith("norm.")]
    gcol = "geneId"
    g = df[gcol].astype(str).str.upper()
    # drop unannotated
    keep = ~g.isin(["NAN", "NONE", ""]) & g.notna()
    df = df.loc[keep]
    g = g.loc[keep]
    mat = df[norms].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    # sum promoters per gene, then mean across samples
    tmp = mat.copy()
    tmp["_gene"] = g.values
    gene_sum = tmp.groupby("_gene")[norms].sum()
    return gene_sum.mean(axis=1)


def module_score(expr: pd.Series, genes: list[str], min_genes: int = 5) -> tuple[float, int]:
    present = [g for g in genes if g in expr.index]
    if len(present) < min_genes:
        return np.nan, len(present)
    vals = np.log2(expr.loc[present].astype(float) + 1.0)
    return float(vals.mean()), len(present)


def fit_loglog(mass: np.ndarray, score: np.ndarray) -> dict:
    x = np.log10(mass)
    y = score
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3:
        return {"n": int(len(x)), "slope": np.nan, "intercept": np.nan, "r": np.nan, "p": np.nan, "ci95": (np.nan, np.nan)}
    slope, intercept, r, p, se = stats.linregress(x, y)
    # approximate 95% CI
    tcrit = stats.t.ppf(0.975, df=max(len(x) - 2, 1))
    ci = (slope - tcrit * se, slope + tcrit * se)
    return {
        "n": int(len(x)),
        "slope": float(slope),
        "intercept": float(intercept),
        "r": float(r),
        "p": float(p),
        "se": float(se),
        "ci95": (float(ci[0]), float(ci[1])),
    }


def pgls_like_ols(mass: np.ndarray, score: np.ndarray) -> dict:
    """Ordinary log-log fit; note true PGLS needs a tree — reported as OLS with caveat."""
    return fit_loglog(mass, score)


def exploratory_subset_search(
    expr_by_sp: dict[str, pd.Series],
    mass_by_sp: dict[str, float],
    gene_pool: list[str],
    target_slope: float = -1.0 / 3.0,
    n_iter: int = 2000,
    subset_size: int = 25,
    seed: int = 42,
) -> dict:
    """Post-hoc search for gene subsets whose module score ~ log10(mass) near target slope.
    EXPLORATORY — not confirmatory.
    """
    rng = np.random.default_rng(seed)
    species = sorted(expr_by_sp)
    # genes present in >=3 species
    counts = {}
    for sp in species:
        for g in gene_pool:
            if g in expr_by_sp[sp].index and expr_by_sp[sp].loc[g] > 0:
                counts[g] = counts.get(g, 0) + 1
    pool = [g for g, c in counts.items() if c >= 3]
    if len(pool) < subset_size:
        subset_size = max(5, len(pool) // 2)
    best = {"abs_err": np.inf, "genes": [], "fit": {}}
    rows = []
    for i in range(n_iter):
        genes = list(rng.choice(pool, size=min(subset_size, len(pool)), replace=False))
        scores, masses = [], []
        for sp in species:
            sc, n = module_score(expr_by_sp[sp], genes, min_genes=3)
            if np.isfinite(sc):
                scores.append(sc)
                masses.append(mass_by_sp[sp])
        fit = fit_loglog(np.array(masses), np.array(scores))
        err = abs(fit["slope"] - target_slope) if np.isfinite(fit["slope"]) else np.inf
        rows.append({"iter": i, "slope": fit["slope"], "p": fit["p"], "r": fit["r"], "abs_err": err})
        if err < best["abs_err"]:
            best = {"abs_err": err, "genes": genes, "fit": fit, "target": target_slope}
    hist = pd.DataFrame(rows)
    return {"best": best, "history": hist, "pool_size": len(pool)}


def main():
    meta_df = pd.DataFrame(
        [
            {
                "species": sp,
                "mass_kg": m["mass_kg"],
                "log10_mass": np.log10(m["mass_kg"]),
                "latin": m["latin"],
                "symbol_ok": m["symbol_ok"],
            }
            for sp, m in SPECIES_META.items()
        ]
    )
    meta_df.to_csv(RES / "species_mass_table.tsv", sep="\t", index=False)

    # --- RNAseq gene means for symbol-ok species ---
    expr = {}
    for sp, m in SPECIES_META.items():
        if not m["symbol_ok"]:
            continue
        path = DATA / m["rnaseq"]
        print(f"Loading RNAseq {sp}: {path.name}")
        expr[sp] = load_rnaseq_gene_mean(path)

    # common ortholog core by shared gene symbol
    common = set.intersection(*(set(s.index) for s in expr.values()))
    common = sorted(g for g in common if re.match(r"^[A-Z][A-Z0-9\-]+$", g))
    print(f"Common symbol core (4 spp): {len(common)}")
    with open(RES / "common_gene_symbols_4spp.txt", "w") as f:
        f.write("\n".join(common) + "\n")

    # save gene-level matrix (common)
    mat = pd.DataFrame({sp: expr[sp].reindex(common).fillna(0.0) for sp in expr})
    mat.to_csv(RES / "rnaseq_gene_tpm_mean_common4.tsv", sep="\t")

    # module coverage
    cov_rows = []
    score_rows = []
    for mod, genes in MODULES.items():
        for sp in expr:
            sc, n = module_score(expr[sp], genes)
            cov_rows.append({"module": mod, "species": sp, "n_genes": n, "score_log2": sc})
            score_rows.append(
                {
                    "module": mod,
                    "species": sp,
                    "mass_kg": SPECIES_META[sp]["mass_kg"],
                    "score_log2": sc,
                    "n_genes": n,
                }
            )
    cov = pd.DataFrame(cov_rows)
    scores = pd.DataFrame(score_rows)
    cov.to_csv(RES / "module_gene_coverage.tsv", sep="\t", index=False)
    scores.to_csv(RES / "module_scores_by_species.tsv", sep="\t", index=False)

    # save gene sets
    for name, genes in MODULES.items():
        (DATA / "genesets" / f"{name}.txt").write_text("\n".join(genes) + "\n")

    # allometry fits
    fit_rows = []
    for mod in MODULES:
        sub = scores[scores.module == mod].dropna(subset=["score_log2"])
        fit = pgls_like_ols(sub.mass_kg.values, sub.score_log2.values)
        fit_rows.append(
            {
                "module": mod,
                **{k: fit[k] for k in ("n", "slope", "intercept", "r", "p", "se")},
                "ci95_lo": fit["ci95"][0],
                "ci95_hi": fit["ci95"][1],
                "target_m13": -1 / 3,
                "target_m14": -1 / 4,
                "abs_err_m13": abs(fit["slope"] + 1 / 3) if np.isfinite(fit["slope"]) else np.nan,
                "abs_err_m14": abs(fit["slope"] + 1 / 4) if np.isfinite(fit["slope"]) else np.nan,
                "note": "OLS on species means; not phylogenetic GLS (n=4)",
            }
        )
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(RES / "allometry_ols_fits.tsv", sep="\t", index=False)
    print(fits.to_string(index=False))

    # figures: allometry
    fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    targets = [(-1 / 3, "M^{-1/3}"), (-1 / 4, "M^{-1/4}")]
    for ax, mod in zip(axes.ravel(), MODULES):
        sub = scores[scores.module == mod].dropna(subset=["score_log2"])
        x = np.log10(sub.mass_kg.values)
        y = sub.score_log2.values
        ax.scatter(x, y, s=60, zorder=3)
        for sp, xi, yi in zip(sub.species, x, y):
            ax.text(xi, yi, sp, fontsize=8, ha="left", va="bottom")
        if len(x) >= 2 and np.all(np.isfinite(y)):
            slope, intercept, r, p, _ = stats.linregress(x, y)
            xs = np.linspace(x.min() - 0.1, x.max() + 0.1, 50)
            ax.plot(xs, intercept + slope * xs, color="black", lw=1.5, label=f"fit b={slope:.2f}")
            y0 = np.nanmean(y)
            x0 = np.nanmean(x)
            for t, lab in targets:
                ax.plot(xs, y0 + t * (xs - x0), ls="--", lw=1, alpha=0.7, label=lab)
            ax.set_title(f"{mod}: b={slope:.2f}, r={r:.2f}, p={p:.3f}")
        else:
            ax.set_title(f"{mod}: insufficient data")
        ax.set_xlabel("log10 adult mass (kg)")
        ax.set_ylabel("module score (mean log2(TPM+1))")
        ax.legend(fontsize=7, loc="best")
    fig.suptitle("Fibroblast RNAseq metabolic modules vs body mass (GSE213253)\nOLS; expect null in culture (Jimenez 2014)", fontsize=11)
    fig.savefig(FIG / "fig_A2_fibroblast_module_allometry.png", dpi=160)
    fig.savefig(FIG / "fig_A2_fibroblast_module_allometry.pdf")
    plt.close(fig)

    # per-gene slopes for OXPHOS pool (common genes)
    gene_slope_rows = []
    mass = np.array([SPECIES_META[sp]["mass_kg"] for sp in expr])
    sps = list(expr)
    for g in sorted(set(OXPHOS_CORE + MITO_BIOGENESIS + GLYCOLYSIS + TCA) & set(common)):
        y = np.array([np.log2(float(expr[sp].get(g, 0.0)) + 1.0) for sp in sps])
        if np.sum(y > 0) < 3:
            continue
        fit = fit_loglog(mass, y)
        gene_slope_rows.append({"gene": g, **{k: fit[k] for k in ("n", "slope", "r", "p")}})
    gene_slopes = pd.DataFrame(gene_slope_rows).sort_values("slope")
    gene_slopes.to_csv(RES / "per_gene_logtpm_vs_logmass_slopes.tsv", sep="\t", index=False)

    # A3 exploratory subset search toward -1/3 and -1/4
    pool = sorted(set(OXPHOS_CORE + MITO_BIOGENESIS + TCA + GLYCOLYSIS))
    mass_map = {sp: SPECIES_META[sp]["mass_kg"] for sp in expr}
    for target, tag in [(-1 / 3, "m13"), (-1 / 4, "m14")]:
        out = exploratory_subset_search(expr, mass_map, pool, target_slope=target, n_iter=3000, subset_size=20)
        best = out["best"]
        (RES / f"exploratory_best_genes_{tag}.txt").write_text("\n".join(sorted(best["genes"])) + "\n")
        out["history"].to_csv(RES / f"exploratory_search_history_{tag}.tsv", sep="\t", index=False)
        summary = {
            "target_slope": target,
            "best_slope": best["fit"].get("slope"),
            "best_p": best["fit"].get("p"),
            "best_r": best["fit"].get("r"),
            "abs_err": best["abs_err"],
            "n_genes": len(best["genes"]),
            "pool_size": out["pool_size"],
            "genes": sorted(best["genes"]),
            "disclaimer": "EXPLORATORY post-hoc gene-subset search; not a confirmatory Kleiber test.",
        }
        with open(RES / f"exploratory_best_summary_{tag}.json", "w") as f:
            json.dump(summary, f, indent=2)
        print(f"Exploratory {tag}: slope={summary['best_slope']}, err={summary['abs_err']:.3f}, genes={summary['n_genes']}")

        # plot best exploratory
        genes = best["genes"]
        xs, ys, labs = [], [], []
        for sp in expr:
            sc, _ = module_score(expr[sp], genes, min_genes=3)
            xs.append(np.log10(SPECIES_META[sp]["mass_kg"]))
            ys.append(sc)
            labs.append(sp)
        fig, ax = plt.subplots(figsize=(5.5, 4.2), constrained_layout=True)
        ax.scatter(xs, ys, s=70)
        for xi, yi, lab in zip(xs, ys, labs):
            ax.text(xi, yi, lab, fontsize=8)
        if np.all(np.isfinite(ys)):
            slope, intercept, r, p, _ = stats.linregress(xs, ys)
            xg = np.linspace(min(xs) - 0.1, max(xs) + 0.1, 40)
            ax.plot(xg, intercept + slope * xg, color="black")
            y0, x0 = np.mean(ys), np.mean(xs)
            ax.plot(xg, y0 + target * (xg - x0), "--", label=f"target {target:.3f}")
            ax.set_title(f"Exploratory subset → target {target:.3f}\nfit b={slope:.3f}, r={r:.2f} (POST-HOC)")
        ax.set_xlabel("log10 mass (kg)")
        ax.set_ylabel("subset score")
        ax.legend(fontsize=8)
        fig.savefig(FIG / f"fig_A3_exploratory_subset_{tag}.png", dpi=160)
        plt.close(fig)

    # CAGE gene-level for symbol-ok (for A6 later reuse)
    cage_expr = {}
    for sp, m in SPECIES_META.items():
        if not m["symbol_ok"]:
            continue
        print(f"Loading CAGE {sp}")
        cage_expr[sp] = load_cage_gene_sum(DATA / m["cage"])
    cage_common = set.intersection(*(set(s.index) for s in cage_expr.values()))
    cage_mat = pd.DataFrame({sp: cage_expr[sp].reindex(sorted(cage_common)).fillna(0.0) for sp in cage_expr})
    cage_mat.to_csv(RES / "cage_gene_promoter_sum_mean_common4.tsv", sep="\t")

    # CAGE OXPHOS allometry (promoter-sum scores)
    cage_score_rows = []
    for mod, genes in MODULES.items():
        for sp in cage_expr:
            sc, n = module_score(cage_expr[sp], genes)
            cage_score_rows.append(
                {"module": mod, "species": sp, "mass_kg": SPECIES_META[sp]["mass_kg"], "score_log2": sc, "n_genes": n}
            )
    cage_scores = pd.DataFrame(cage_score_rows)
    cage_scores.to_csv(RES / "cage_module_scores_by_species.tsv", sep="\t", index=False)
    cage_fits = []
    for mod in MODULES:
        sub = cage_scores[cage_scores.module == mod].dropna(subset=["score_log2"])
        fit = pgls_like_ols(sub.mass_kg.values, sub.score_log2.values)
        cage_fits.append({"module": mod, "assay": "CAGE_promoter_sum", **{k: fit[k] for k in ("n", "slope", "r", "p")}})
    pd.DataFrame(cage_fits).to_csv(RES / "cage_allometry_ols_fits.tsv", sep="\t", index=False)

    # note on missing species
    notes = {
        "symbol_ok_species": [sp for sp, m in SPECIES_META.items() if m["symbol_ok"]],
        "no_portable_symbols": [sp for sp, m in SPECIES_META.items() if not m["symbol_ok"]],
        "reason": "GEO RNAseq/CAGE for Peacock/Turkey/Emu use assembly gene IDs (g####), not gene symbols; excluded from cross-species symbol intersection.",
        "mass_sources": "Adult body mass (kg) from standard avian references (Dunning-style handbook values); not measured in GSE213253.",
        "exponents": "Mass-specific cellular rate targets are M^{-1/3} and M^{-1/4}; organismal BMR M^{2/3|3/4} is a different level.",
    }
    with open(RES / "A1_ingest_notes.json", "w") as f:
        json.dump(notes, f, indent=2)
    print("Done A1–A3 core.")


if __name__ == "__main__":
    main()
