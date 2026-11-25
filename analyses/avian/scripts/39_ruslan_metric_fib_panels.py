#!/usr/bin/env python3
"""Ruslan metric on every local cultured-fibroblast panel.

Metric: nuclear OXPHOS after alias collapse; per library, mean TPM of genes
with TPM>0; species = mean of libraries; Pearson of log10(mean TPM) vs
log10(adult mass). Panels are never pooled.

Local only: GSE213253 CAGE, GSE213253 GEO RNA, FANTOM5 human/mouse fib,
Wagner salmon SF, GSE167580 primary-like. Ma / whale / elephant are not here.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    EXP,
    PREFERRED_ALIAS,
    RES,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    load_geo_rnaseq,
    nuclear_oxphos,
    normalize_symbol,
    ols_logmass,
    peak_to_symbol_de,
    spearman_logmass,
)
from oxphos_thread_lib import RUSLAN  # noqa: E402

OUT = RES / "ruslan_metric_fib_panels"


def oxphos_targets() -> set[str]:
    genes = nuclear_oxphos()
    extra = {PREFERRED_ALIAS.get(g, g) for g in genes}
    return {normalize_symbol(g) for g in genes | extra if normalize_symbol(g)}


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def mean_present(vals: pd.Series) -> tuple[float, int]:
    v = pd.to_numeric(vals, errors="coerce").fillna(0.0)
    present = v[v > 0]
    if present.empty:
        return float("nan"), 0
    return float(present.mean()), int(len(present))


def fit_species(df: pd.DataFrame, panel: str) -> dict:
    y = np.log10(df["mean_tpm"].astype(float))
    fit = ols_logmass(pd.Series(y.values), df["mass_g"])
    sp = spearman_logmass(pd.Series(y.values), df["mass_g"])
    return {
        "panel": panel,
        "n": fit["n"],
        "r": fit["r"],
        "p": fit["p"],
        "slope": fit["slope"],
        "rho": sp["rho"],
        "p_spearman": sp["p"],
        "species": ",".join(df["species"].astype(str).tolist()),
    }


def birds_cage(genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    libs, spec = [], []
    for sp in SPECIES_ORDER:
        meta = SPECIES[sp]
        exp = pd.read_csv(EXP / meta["file"], sep="\t")
        exp = annotate_table(exp, peak_to_symbol_de(RUSLAN / meta["de"]))
        cols = lib_cols(exp, "fib")
        sub = exp[exp["symbol"].isin(genes)].copy()
        mat = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        sub = sub.assign(**{c: mat[c] for c in cols})
        per = []
        ngs = []
        for c in cols:
            gmax = sub.groupby("symbol")[c].max()
            mu, n = mean_present(gmax)
            libs.append(
                {
                    "panel": "gse213253_cage",
                    "species": sp.lower(),
                    "sample": c,
                    "mass_g": meta["mass_g"],
                    "n_genes": n,
                    "mean_tpm": mu,
                }
            )
            per.append(mu)
            ngs.append(n)
        spec.append(
            {
                "panel": "gse213253_cage",
                "species": sp.lower(),
                "mass_g": meta["mass_g"],
                "n_libs": len(cols),
                "n_genes_med": int(np.median(ngs)) if ngs else 0,
                "mean_tpm": float(np.nanmean(per)),
                "agg": "max_peak",
            }
        )
    return pd.DataFrame(libs), pd.DataFrame(spec)


def birds_rna(genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    libs, spec = [], []
    for sp in SPECIES_ORDER:
        meta = SPECIES[sp]
        packed = load_geo_rnaseq(sp)
        mat = packed["gene_mat"]
        hit = [g for g in genes if g in mat.index]
        per, ngs = [], []
        for c in packed["norm_cols"]:
            mu, n = mean_present(mat.loc[hit, c] if hit else pd.Series(dtype=float))
            libs.append(
                {
                    "panel": "gse213253_rna",
                    "species": sp.lower(),
                    "sample": c,
                    "mass_g": meta["mass_g"],
                    "n_genes": n,
                    "mean_tpm": mu,
                }
            )
            per.append(mu)
            ngs.append(n)
        spec.append(
            {
                "panel": "gse213253_rna",
                "species": sp.lower(),
                "mass_g": meta["mass_g"],
                "n_libs": len(packed["norm_cols"]),
                "n_genes_med": int(np.median(ngs)) if ngs else 0,
                "mean_tpm": float(np.nanmean(per)),
                "agg": "geo_norm_mean",
                "note": packed["mapped_how"],
            }
        )
    return pd.DataFrame(libs), pd.DataFrame(spec)


def fantom_fib(genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    s38 = _load("s38", Path(__file__).resolve().parent / "38_fantom_msc_oxphos.py")
    FANTOM = s38.FANTOM
    libs, spec = [], []
    blocks = [
        ("human", s38.MASS["human"], FANTOM / "hg19.gene_phase1and2combined_tpm.osc.txt.gz", s38.HUMAN_DERMAL),
        ("mouse", s38.MASS["mouse"], FANTOM / "mm9.gene_phase1and2combined_tpm.osc.txt.gz", s38.MOUSE_MEF),
    ]
    for species, mass, path, ids in blocks:
        rows = s38.score_gene_level(path, ids, genes)
        mus, ngs = [], []
        for r in rows:
            libs.append(
                {
                    "panel": "fantom5_fib",
                    "species": species,
                    "sample": r["cnhs"],
                    "mass_g": mass,
                    "n_genes": r["n_genes_any"],
                    "mean_tpm": r["mean_sum"],
                }
            )
            mus.append(r["mean_sum"])
            ngs.append(r["n_genes_any"])
        if mus:
            spec.append(
                {
                    "panel": "fantom5_fib",
                    "species": species,
                    "mass_g": mass,
                    "n_libs": len(mus),
                    "n_genes_med": int(np.median(ngs)),
                    "mean_tpm": float(np.nanmean(mus)),
                    "agg": "gene_sum",
                }
            )
    return pd.DataFrame(libs), pd.DataFrame(spec)


def wagner(genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    s37 = _load("s37", Path(__file__).resolve().parent / "37_wagner_sf_oxphos.py")
    libs = []
    for sp, sample, _srr in s37.SAMPLES:
        qf = s37.QUANT / sample / "quant.sf"
        if not qf.exists():
            continue
        t2g = pd.read_csv(s37.build_tx2gene(sp), sep="\t")
        q = pd.read_csv(qf, sep="\t")
        q = q.merge(t2g, left_on="Name", right_on="tx", how="left")
        q["symbol"] = q["symbol"].fillna("").map(lambda s: normalize_symbol(s) if s else "")
        if "raw_symbol" in q.columns:
            miss = q["symbol"].eq("")
            q.loc[miss, "symbol"] = q.loc[miss, "raw_symbol"].fillna("").map(normalize_symbol)
        by = q[q["symbol"].isin(genes)].groupby("symbol", observed=True)["TPM"].sum()
        mu, n = mean_present(by)
        libs.append(
            {
                "panel": "wagner_sf",
                "species": sp,
                "sample": sample,
                "mass_g": s37.MASS_G[sp],
                "n_genes": n,
                "mean_tpm": mu,
            }
        )
    lib = pd.DataFrame(libs)
    spec = (
        lib.groupby("species", as_index=False)
        .agg(n_libs=("sample", "size"), n_genes_med=("n_genes", "median"), mean_tpm=("mean_tpm", "mean"), mass_g=("mass_g", "first"))
        .assign(panel="wagner_sf", agg="salmon_gene_sum")
    )
    spec["n_genes_med"] = spec["n_genes_med"].astype(int)
    return lib, spec


def gse167580(genes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    s36 = _load("s36", Path(__file__).resolve().parent / "36_gse167580_fib_oxphos_module.py")
    libs = []
    for sp, fn in s36.FILES.items():
        df = pd.read_csv(s36.FIB / fn, sep="\t")
        id2, _how = s36.id2symbol(sp, genes)
        ids = df["id"].astype(str)
        sample_cols = [c for c in df.columns if c != "id"]
        mat = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        mapped = ids.map(id2).fillna("")
        keep = mapped.isin(genes)
        sub = mat.loc[keep].copy()
        sub["_sym"] = mapped[keep].values
        collapsed = sub.groupby("_sym", observed=True)[sample_cols].max()
        for col in sample_cols:
            klass = s36.SAMPLE_CLASS.get((sp, col), "unknown")
            if klass != "primary_like":
                continue
            mu, n = mean_present(collapsed[col] if col in collapsed.columns else pd.Series(dtype=float))
            libs.append(
                {
                    "panel": "gse167580_primary",
                    "species": sp,
                    "sample": col,
                    "mass_g": s36.MASS_G[sp],
                    "n_genes": n,
                    "mean_tpm": mu,
                }
            )
    lib = pd.DataFrame(libs)
    spec = (
        lib.groupby("species", as_index=False)
        .agg(n_libs=("sample", "size"), n_genes_med=("n_genes", "median"), mean_tpm=("mean_tpm", "mean"), mass_g=("mass_g", "first"))
        .assign(panel="gse167580_primary", agg="kallisto_gene_max")
    )
    spec["n_genes_med"] = spec["n_genes_med"].astype(int)
    return lib, spec


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = oxphos_targets()
    print(f"nuclear OXPHOS after alias collapse: {len(genes)}", flush=True)

    blocks = [
        birds_cage(genes),
        birds_rna(genes),
        fantom_fib(genes),
        wagner(genes),
        gse167580(genes),
    ]
    libs = pd.concat([a for a, _ in blocks], ignore_index=True)
    spec = pd.concat([b for _, b in blocks], ignore_index=True)
    spec = spec.sort_values(["panel", "mass_g"])
    libs.to_csv(OUT / "libraries.tsv", sep="\t", index=False)
    spec.to_csv(OUT / "species.tsv", sep="\t", index=False)

    fits = []
    for panel, g in spec.groupby("panel", sort=False):
        g = g.sort_values("mass_g")
        if len(g) >= 4:
            fits.append(fit_species(g, panel))
        else:
            fits.append(
                {
                    "panel": panel,
                    "n": int(len(g)),
                    "r": np.nan,
                    "p": np.nan,
                    "slope": np.nan,
                    "rho": np.nan,
                    "p_spearman": np.nan,
                    "species": ",".join(g["species"].astype(str)),
                }
            )
    fits = pd.DataFrame(fits)
    fits.to_csv(OUT / "fits.tsv", sep="\t", index=False)
    print(spec.to_string(index=False))
    print(fits.to_string(index=False))


if __name__ == "__main__":
    main()
