#!/usr/bin/env python3
"""Matched nulls for the OXPHOS-vs-mass decline.

Script 24 showed the plain random-set null is not neutral: in the GEO RNA tables
a random gene set already declines with mass (median r ~ -0.45), and the OXPHOS
signal in CAGE disappears once the gene list is restricted to genes detected in
all seven species. Both facts can be produced by detection completeness or by
expression level rather than by OXPHOS biology.

This script draws null gene sets matched to the tested set on
  (a) detection breadth  - in how many of the 7 species the gene is seen
  (b) abundance          - chicken expression decile
and reports the observed statistic against those matched nulls, for both the
fraction metric and the Y1 (mean over detected genes) metric.
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
from oxphos_thread_lib import RES, SPECIES, SPECIES_ORDER, load_geneset, normalize_symbol, nuclear_oxphos  # noqa: E402

OUT = RES / "reverify_confounds"
CACHE = OUT / "gene_matrices.pkl"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
N_NULL = 4000
RNG = np.random.default_rng(4242)

MASS = pd.Series({sp: SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER})
LOGM = np.log10(MASS.reindex(SPECIES_ORDER).values.astype(float))


def fit_r(y: np.ndarray) -> float:
    ok = np.isfinite(y)
    if ok.sum() < 4 or np.ptp(y[ok]) == 0:
        return np.nan
    return float(stats.linregress(LOGM[ok], y[ok]).rvalue)


def fit_r_p(y: np.ndarray) -> tuple[float, float]:
    ok = np.isfinite(y)
    if ok.sum() < 4 or np.ptp(y[ok]) == 0:
        return np.nan, np.nan
    res = stats.linregress(LOGM[ok], y[ok])
    return float(res.rvalue), float(res.pvalue)


def metrics(layer: dict, genes: np.ndarray) -> dict[str, np.ndarray]:
    frac, y1 = [], []
    for sp in SPECIES_ORDER:
        gmean, tot = layer[sp]
        hit = gmean.reindex(genes).dropna()
        frac.append(float(hit.sum()) / tot)
        y1.append(float(hit.mean()) if len(hit) else np.nan)
    return {"fraction": np.array(frac), "Y1": np.array(y1)}


def build_gene_profile(layer: dict) -> pd.DataFrame:
    """Per gene: detection breadth across 7 species + chicken abundance decile."""
    breadth = {}
    for sp in SPECIES_ORDER:
        for g in layer[sp][0].index:
            breadth[g] = breadth.get(g, 0) + 1
    prof = pd.DataFrame({"symbol": list(breadth.keys()), "breadth": list(breadth.values())})
    chick = layer["Chicken"][0]
    prof["abundance"] = prof["symbol"].map(chick).astype(float)
    seen = prof["abundance"].notna()
    prof["abund_decile"] = -1
    if seen.sum() > 20:
        prof.loc[seen, "abund_decile"] = pd.qcut(
            prof.loc[seen, "abundance"].rank(method="first"), 10, labels=False).astype(int)
    return prof.set_index("symbol")


def matched_draw(prof: pd.DataFrame, target: np.ndarray, keys: list[str]) -> np.ndarray | None:
    """One null set with the same joint distribution of `keys` as `target`."""
    tgt = prof.reindex(target).dropna(subset=keys)
    if tgt.empty:
        return None
    pool = prof.drop(index=[g for g in target if g in prof.index], errors="ignore")
    out = []
    for _, grp in tgt.groupby(keys, observed=True):
        mask = np.ones(len(pool), dtype=bool)
        for k, v in zip(keys, (grp.iloc[0][k] for k in keys)):
            mask &= (pool[k].values == v)
        cand = pool.index.values[mask]
        if len(cand) == 0:
            return None
        out.append(RNG.choice(cand, size=len(grp), replace=len(cand) < len(grp)))
    return np.concatenate(out)


def main() -> None:
    with CACHE.open("rb") as fh:
        layers = pickle.load(fh)

    oxphos = np.array(sorted({normalize_symbol(g) for g in nuclear_oxphos()} - {""}))
    module = np.array(sorted({normalize_symbol(g) for g in
                              pd.read_csv(RES / "definitions" / "oxphos_module.tsv", sep="\t")["gene"]} - {""}))
    glyco = np.array(sorted({normalize_symbol(g) for g in load_geneset("GLYCOLYSIS")} - {""}))

    rows, dists = [], {}
    for lname, layer in layers.items():
        prof = build_gene_profile(layer)
        prof.to_csv(OUT / f"gene_profile_{lname}.tsv", sep="\t")
        for setname, genes in (("set", oxphos), ("module", module), ("GLYCOLYSIS", glyco)):
            present = np.array([g for g in genes if g in prof.index])
            for null_kind, keys in (("unmatched", None),
                                    ("matched_breadth", ["breadth"]),
                                    ("matched_breadth_abundance", ["breadth", "abund_decile"])):
                for metric in ("fraction", "Y1"):
                    obs_r, obs_p = fit_r_p(metrics(layer, present)[metric])
                    null = []
                    fails = 0
                    for _ in range(N_NULL):
                        if keys is None:
                            draw = RNG.choice(prof.index.values, size=len(present), replace=False)
                        else:
                            draw = matched_draw(prof, present, keys)
                            if draw is None:
                                fails += 1
                                continue
                        r = fit_r(metrics(layer, draw)[metric])
                        if np.isfinite(r):
                            null.append(r)
                    null = np.array(null)
                    if len(null) < 100:
                        continue
                    emp = float((null <= obs_r).mean())
                    dists[f"{lname}|{setname}|{null_kind}|{metric}"] = null
                    rows.append(dict(layer=lname, definition=setname, metric=metric,
                                     null_kind=null_kind, n_genes=len(present),
                                     observed_r=obs_r, naive_ols_p=obs_p,
                                     null_median_r=float(np.median(null)),
                                     null_q05_r=float(np.quantile(null, 0.05)),
                                     empirical_p=emp, n_null=len(null), n_draw_failed=fails))
                    print(f"  {lname:9s} {setname:11s} {metric:8s} {null_kind:26s} "
                          f"obs r={obs_r:+.3f} null med={np.median(null):+.3f} emp p={emp:.4f}", flush=True)

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "matched_null_results.tsv", sep="\t", index=False)
    np.savez_compressed(OUT / "matched_null_dists.npz", **dists)

    # figure: headline claims against the strongest null
    show = [("CAGE_fib", "set", "fraction"), ("CAGE_fib", "set", "Y1"),
            ("CAGE_fib", "module", "Y1"), ("GEO_RNA", "set", "fraction"),
            ("GEO_RNA", "module", "fraction"), ("CAGE_epi", "set", "fraction")]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for ax, (lname, setname, metric) in zip(axes.ravel(), show):
        key = f"{lname}|{setname}|matched_breadth_abundance|{metric}"
        if key not in dists:
            key = f"{lname}|{setname}|matched_breadth|{metric}"
        null = dists[key]
        r = res[(res.layer == lname) & (res.definition == setname) & (res.metric == metric)
                & (res.null_kind == key.split("|")[2])].iloc[0]
        ax.hist(null, bins=45, color="#adb5bd", edgecolor="none")
        ax.axvline(r.observed_r, color="#c1121f", lw=2)
        ax.axvline(float(np.median(null)), color="#264653", lw=1.2, ls="--")
        ax.set_title(f"{lname} {setname} {metric}\nnaive p={r.naive_ols_p:.4g} → matched p={r.empirical_p:.3f}",
                     fontsize=9)
        ax.set_xlabel("OLS r vs log10 mass")
    fig.suptitle("OXPHOS decline against nulls matched on detection breadth and abundance", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_matched_null.png", dpi=160)
    fig.savefig(OUT / "fig_oxphos_matched_null.png", dpi=160)
    plt.close()

    pd.set_option("display.width", 220, "display.max_columns", 30)
    print("\n=== matched nulls ===")
    print(res[["layer", "definition", "metric", "null_kind", "n_genes", "observed_r",
               "naive_ols_p", "null_median_r", "empirical_p"]].to_string(index=False))
    (OUT / "matched_null_summary.json").write_text(json.dumps(
        {"n_null_requested": N_NULL, "seed": 4242,
         "headline": res[res.null_kind == "matched_breadth_abundance"].to_dict(orient="records")},
        indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
