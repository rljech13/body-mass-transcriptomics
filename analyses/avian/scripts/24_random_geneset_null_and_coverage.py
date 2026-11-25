#!/usr/bin/env python3
"""Is the OXPHOS-vs-mass decline specific to OXPHOS, or does any gene set do it?

The corpus tests OXPHOS against exactly one control set (GLYCOLYSIS, 16 genes).
Here the same phenotype recipe is applied to thousands of random gene sets drawn
from the same annotation, which gives the null distribution the claim needs.

Also decomposes the primary metric into its arithmetic parts:
    fraction = n_detected * mean_over_detected / library_total
and regresses each part on log10 body mass separately.

Layers: CAGE fib, CAGE epi (native Dropbox tables), GEO RNA-seq (GSE213253).
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    EXP,
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    load_geneset,
    load_geo_rnaseq,
    normalize_symbol,
    nuclear_oxphos,
    peak_to_symbol_de,
)

OUT = RES / "reverify_confounds"
CACHE = OUT / "gene_matrices.pkl"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
N_RANDOM = 3000
RNG = np.random.default_rng(20260820)

MASS = pd.Series({sp: SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER})
LOGM = np.log10(MASS.reindex(SPECIES_ORDER).values.astype(float))


def build_layers() -> dict:
    """{layer: {species: (gene_mean_series, library_total)}}"""
    if CACHE.exists():
        with CACHE.open("rb") as fh:
            return pickle.load(fh)
    layers = {"CAGE_fib": {}, "CAGE_epi": {}, "GEO_RNA": {}}
    for sp, meta in SPECIES.items():
        df = pd.read_csv(EXP / meta["file"], sep="\t")
        df = annotate_table(df, peak_to_symbol_de(RUSLAN / meta["de"]))
        df["symbol"] = df["symbol"].map(normalize_symbol)
        prom = df[df["class"].astype(str).str.lower().eq("promoter")]
        for tiss, key in (("fib", "CAGE_fib"), ("epi", "CAGE_epi")):
            cols = lib_cols(prom, tiss)
            mat = prom[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
            lib_total = float(mat.sum(axis=0).mean())
            keep = prom["symbol"].astype(str).str.len() > 0
            g = mat.loc[keep].copy()
            g["_sym"] = prom.loc[keep, "symbol"].values
            gene_mean = g.groupby("_sym", observed=True)[cols].sum().mean(axis=1)
            layers[key][sp] = (gene_mean, lib_total)
            print(f"  {sp} {key}: {len(gene_mean)} genes, total={lib_total:,.0f}", flush=True)
        geo = load_geo_rnaseq(sp)
        layers["GEO_RNA"][sp] = (geo["gene_mat"].mean(axis=1), float(geo["lib_total"].mean()))
        print(f"  {sp} GEO_RNA: {geo['n_symbols']} genes ({geo['mapped_how']})", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    with CACHE.open("wb") as fh:
        pickle.dump(layers, fh)
    return layers


def frac_for_set(layer: dict, genes: np.ndarray) -> np.ndarray:
    """Fraction of the library taken by a gene set; undetected genes contribute 0."""
    out = []
    for sp in SPECIES_ORDER:
        gmean, tot = layer[sp]
        out.append(float(gmean.reindex(genes).fillna(0.0).sum()) / tot)
    return np.array(out)


def fit_r_p(y: np.ndarray) -> tuple[float, float]:
    ok = np.isfinite(y)
    if ok.sum() < 4 or np.ptp(y[ok]) == 0:
        return np.nan, np.nan
    sl, ic, r, p, se = stats.linregress(LOGM[ok], y[ok])
    return float(r), float(p)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("=== building gene matrices ===", flush=True)
    layers = build_layers()

    oxphos = np.array(sorted({normalize_symbol(g) for g in nuclear_oxphos()} - {""}))
    glyco = np.array(sorted({normalize_symbol(g) for g in load_geneset("GLYCOLYSIS")} - {""}))
    module = np.array(sorted({normalize_symbol(g) for g in
                              pd.read_csv(RES / "definitions" / "oxphos_module.tsv", sep="\t")["gene"]} - {""}))
    print(f"canonical: OXPHOS {len(oxphos)}, module {len(module)}, glycolysis {len(glyco)}", flush=True)

    # ---------- 1. arithmetic decomposition of the primary metric ----------
    dec_rows = []
    for lname, layer in layers.items():
        for setname, genes in (("set", oxphos), ("module", module), ("GLYCOLYSIS", glyco)):
            n_det, mean_det, tot, frac = [], [], [], []
            for sp in SPECIES_ORDER:
                gmean, t = layer[sp]
                hit = gmean.reindex(genes).dropna()
                n_det.append(len(hit))
                mean_det.append(float(hit.mean()) if len(hit) else np.nan)
                tot.append(t)
                frac.append(float(hit.sum()) / t)
            for part, vals in (("n_detected", n_det), ("mean_over_detected", mean_det),
                               ("library_total", tot), ("fraction", frac)):
                r, p = fit_r_p(np.log10(np.array(vals, dtype=float)))
                dec_rows.append(dict(layer=lname, definition=setname, component=part,
                                     log10_r=r, log10_p=p,
                                     values=json.dumps([round(float(v), 6) for v in vals])))
    dec = pd.DataFrame(dec_rows)
    dec.to_csv(OUT / "metric_decomposition.tsv", sep="\t", index=False)

    # ---------- 2. random gene-set null ----------
    null_rows, null_dists = [], {}
    for lname, layer in layers.items():
        universe = set()
        for sp in SPECIES_ORDER:
            universe |= set(layers[lname][sp][0].index)
        universe -= set(oxphos) | set(glyco) | set(module)
        universe = np.array(sorted(universe))
        for setname, genes in (("set", oxphos), ("module", module), ("GLYCOLYSIS", glyco)):
            obs_r, obs_p = fit_r_p(frac_for_set(layer, genes))
            rs = np.empty(N_RANDOM)
            for i in range(N_RANDOM):
                draw = RNG.choice(universe, size=len(genes), replace=False)
                rs[i] = fit_r_p(frac_for_set(layer, draw))[0]
            rs = rs[np.isfinite(rs)]
            null_dists[f"{lname}|{setname}"] = rs
            pct = float((rs <= obs_r).mean())
            null_rows.append(dict(
                layer=lname, definition=setname, n_genes=len(genes),
                observed_r=obs_r, observed_ols_p=obs_p,
                null_median_r=float(np.median(rs)), null_mean_r=float(np.mean(rs)),
                null_frac_r_below_minus0p8=float((rs < -0.8).mean()),
                null_frac_ols_significant=float(np.nan),
                empirical_p_one_sided=pct,
                n_null=len(rs)))
            print(f"  {lname:9s} {setname:11s} obs r={obs_r:+.3f} (p={obs_p:.4g}) | "
                  f"null median r={np.median(rs):+.3f}, empirical p={pct:.3f}", flush=True)
    null = pd.DataFrame(null_rows)
    null.to_csv(OUT / "random_geneset_null.tsv", sep="\t", index=False)
    np.savez_compressed(OUT / "random_geneset_null_dists.npz", **null_dists)

    # ---------- 3. fixed-list threshold gradient ----------
    grad_rows = []
    for lname, layer in layers.items():
        for setname, genes in (("set", oxphos), ("module", module)):
            counts = pd.Series(0, index=genes, dtype=int)
            for sp in SPECIES_ORDER:
                counts += pd.Series(1, index=genes).where(
                    pd.Series(genes, index=genes).isin(layer[sp][0].index), 0)
            for k in range(1, 8):
                keep = counts[counts >= k].index.to_numpy()
                if len(keep) < 4:
                    continue
                r, p = fit_r_p(frac_for_set(layer, keep))
                grad_rows.append(dict(layer=lname, definition=setname,
                                      min_species_with_gene=k, n_genes_in_list=len(keep),
                                      ols_r=r, ols_p=p))
    grad = pd.DataFrame(grad_rows)
    grad.to_csv(OUT / "detection_threshold_gradient.tsv", sep="\t", index=False)

    # ---------- figure ----------
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.2), constrained_layout=True)
    for j, lname in enumerate(("CAGE_fib", "CAGE_epi", "GEO_RNA")):
        ax = axes[0, j]
        rs = null_dists[f"{lname}|set"]
        obs = null[(null.layer == lname) & (null.definition == "set")].iloc[0]
        ax.hist(rs, bins=45, color="#adb5bd", edgecolor="none")
        ax.axvline(obs.observed_r, color="#c1121f", lw=2)
        ax.axvline(float(np.median(rs)), color="#264653", lw=1.2, ls="--")
        ax.set_title(f"{lname}: random 65-gene sets\nOXPHOS r={obs.observed_r:+.2f}, "
                     f"empirical p={obs.empirical_p_one_sided:.3f}", fontsize=9)
        ax.set_xlabel("OLS r of set fraction vs log10 mass")
        if j == 0:
            ax.set_ylabel("random gene sets")

        ax = axes[1, j]
        g = grad[(grad.layer == lname) & (grad.definition == "set")]
        ax.plot(g.min_species_with_gene, -np.log10(g.ols_p), "o-", color="#1d3557")
        ax.axhline(-np.log10(0.05), color="#c1121f", ls=":", lw=1.2)
        for _, row in g.iterrows():
            ax.annotate(f"{int(row.n_genes_in_list)}g", (row.min_species_with_gene, -np.log10(row.ols_p)),
                        textcoords="offset points", xytext=(4, 4), fontsize=7)
        ax.set_title(f"{lname}: signal vs how many species must have the gene", fontsize=9)
        ax.set_xlabel("gene kept if detected in ≥ k of 7 species")
        if j == 0:
            ax.set_ylabel("-log10 p (fraction vs mass)")
    fig.suptitle("OXPHOS-vs-mass decline against a random gene-set null and against detection completeness",
                 fontsize=11)
    fig.savefig(FIG / "fig_oxphos_random_null_and_coverage.png", dpi=160)
    fig.savefig(OUT / "fig_oxphos_random_null_and_coverage.png", dpi=160)
    plt.close()

    pd.set_option("display.width", 200, "display.max_columns", 30)
    print("\n=== decomposition (log10 of each component vs log10 mass) ===")
    print(dec[dec.component != "fraction"][["layer", "definition", "component", "log10_r", "log10_p"]]
          .to_string(index=False))
    print("\n=== random gene-set null ===")
    print(null[["layer", "definition", "n_genes", "observed_r", "observed_ols_p",
                "null_median_r", "empirical_p_one_sided"]].to_string(index=False))
    print("\n=== detection-threshold gradient ===")
    print(grad.to_string(index=False))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
