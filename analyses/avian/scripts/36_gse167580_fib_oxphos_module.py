#!/usr/bin/env python3
"""Fibroblast OXPHOS module fraction vs mass: GSE167580 + GSE213253 birds.

Same recipe that produced the bird GEO RNA signal: frozen 29-gene CEF
coexpression module, within-table transcriptome fraction
(sum of mapped module TPM / sum of all TPM). Studies are never pooled.
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
    PREFERRED_ALIAS,
    RES,
    SPECIES,
    SPECIES_ORDER,
    load_geo_rnaseq,
    normalize_symbol,
    ols_logmass,
    spearman_logmass,
)

OUT = RES / "gse167580_fib_module"
FIB = DATA / "gse167580_fib"
ORTH = DATA / "ensembl_orthologs"
GTF = DATA / "old_gtf"

# AnAge adult body mass (g) unless noted. Birds from oxphos_thread_lib.
MASS_G = {
    "mouse": 20.5,
    "rat": 300.0,
    "rabbit": 1800.0,
    "cat": 4000.0,
    "rhesus": 8000.0,
    "dog": 40000.0,
    "human": 62035.0,
    "sheep": 80000.0,  # AnAge-typical domestic sheep ~80 kg
    "pig": 130000.0,
    "cow": 750000.0,
    "chicken": 912.5,
}
MASS_G.update({sp.lower(): SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER})

LATIN = {
    "mouse": "M. musculus",
    "rat": "R. norvegicus",
    "rabbit": "O. cuniculus",
    "cat": "F. catus",
    "rhesus": "M. mulatta",
    "dog": "C. familiaris",
    "human": "H. sapiens",
    "sheep": "O. aries",
    "pig": "S. scrofa",
    "cow": "B. taurus",
    "chicken": "G. gallus",
    "zebra": "T. guttata",
    "quail": "C. japonica",
    "duck": "A. platyrhynchos",
    "peacock": "P. cristatus",
    "turkey": "M. gallopavo",
    "emu": "D. novaehollandiae",
}

FILES = {
    "human": "GSE167580_human_quant_TPM.txt.gz",
    "mouse": "GSE167580_mouse_quant_TPM-1.txt.gz",
    "rat": "GSE167580_rat_quant_TPM.txt.gz",
    "rabbit": "GSE167580_rabbit_quant_TPM-1.txt.gz",
    "dog": "GSE167580_dog_quant_TPM.txt.gz",
    "pig": "GSE167580_pig_quant_TPM-1.txt.gz",
    "cow": "GSE167580_cow_quant_TPM.txt.gz",
    "sheep": "GSE167580_sheep_quant_TPM.txt.gz",
    "rhesus": "GSE167580_monkey_quant_TPM.txt.gz",
    "cat": "GSE167580_cat_quant_TPM.txt.gz",
    "chicken": "GSE167580_chicken_quant_TPM.txt.gz",
}

# He et al. sample flags. Rat BRL-3A is a liver epithelial line, not a fibroblast.
SAMPLE_CLASS = {
    ("human", "ESF"): "primary_like",
    ("human", "HPF"): "primary_like",
    ("human", "WI38"): "line",
    ("mouse", "MEF"): "primary_like",
    ("mouse", "3T6"): "line",
    ("mouse", "PA12"): "line",
    ("rat", "BRL_3A"): "not_fibroblast",
    ("rabbit", "rabbit1"): "primary_like",
    ("rabbit", "rabbit3"): "primary_like",
    ("rabbit", "R9ab"): "line",
    ("dog", "BG-1"): "primary_like",
    ("dog", "BG-2"): "primary_like",
    ("dog", "A-72"): "tumor",
    ("pig", "DB-2"): "primary_like",
    ("pig", "DB-3"): "primary_like",
    ("pig", "RC-7"): "primary_like",
    ("pig", "RC-8"): "primary_like",
    ("cow", "HSF-3"): "primary_like",
    ("cow", "HSF-4"): "primary_like",
    ("sheep", "XWHY-1"): "primary_like",
    ("sheep", "XWHY2"): "primary_like",
    ("sheep", "XWHY3"): "primary_like",
    ("rhesus", "DBS"): "primary_like",
    ("cat", "CEF-1"): "primary_like",
    ("cat", "CEF-2"): "primary_like",
    ("chicken", "BLH1"): "primary_like",
    ("chicken", "BLH4"): "primary_like",
    ("chicken", "BLH9"): "primary_like",
    ("chicken", "BLH12"): "primary_like",
    ("chicken", "DF1"): "line",
}

BIOMART = {
    "mouse": ORTH / "mouse_to_human.tsv",
    "rat": ORTH / "rat_to_human.tsv",
    "rabbit": ORTH / "rabbit_to_human.tsv",
    "pig": ORTH / "pig_to_human.tsv",
    "cow": ORTH / "cow_to_human.tsv",
    "rhesus": ORTH / "rhesus_to_human.tsv",
}

GTF_MAP = {
    "dog": GTF / "dog_canfam3_genes.tsv",
    "cat": GTF / "cat_felcat62_genes.tsv",
    "chicken": GTF / "chicken_gg5_genes.tsv",
    "sheep": GTF / "sheep_oar31_genes.tsv",
}


def load_defs() -> dict[str, set[str]]:
    mod = {normalize_symbol(g) for g in pd.read_csv(RES / "definitions" / "oxphos_module.tsv", sep="\t")["gene"]}
    raw = pd.read_csv(RES / "definitions" / "oxphos_set.tsv", sep="\t")["gene"]
    oxset = {normalize_symbol(g) for g in raw}
    oxset |= set(PREFERRED_ALIAS.values()) & {normalize_symbol(g) for g in raw}
    # set file already mixes old/new aliases; keep unique normalized nuclear names
    oxset = {g for g in oxset if g}
    mod = {g for g in mod if g}
    return {"module": mod, "set": oxset}


def _biomart_id2sym(path: Path, wanted: set[str]) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t")
    df["sym"] = df["Human gene name"].map(normalize_symbol)
    df = df[df["sym"] != ""]
    one2 = df["Human homology type"].astype(str).str.contains("one2one", na=False)
    out: dict[str, str] = {}
    for prefer in (True, False):
        sub = df[one2] if prefer else df[~one2]
        for eid, sym in zip(sub["Gene stable ID"].astype(str), sub["sym"]):
            if eid in out:
                continue
            if prefer or sym in wanted:
                out[eid] = sym
    return out


def _gtf_id2sym(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t")
    out = {}
    for eid, name in zip(df["ensembl_gene_id"].astype(str), df["gene_name"].astype(str)):
        sym = normalize_symbol(name)
        if sym:
            out[eid] = sym
    return out


def _human_id2sym() -> dict[str, str]:
    df = pd.read_csv(ORTH / "human_ensg_symbol.tsv", sep="\t")
    out = {}
    for eid, name in zip(df["Gene stable ID"].astype(str), df["Gene name"].astype(str)):
        sym = normalize_symbol(name)
        if sym:
            out[eid] = sym
    return out


def id2symbol(species: str, wanted: set[str]) -> tuple[dict[str, str], str]:
    if species == "human":
        return _human_id2sym(), "biomart_symbol"
    if species in BIOMART:
        return _biomart_id2sym(BIOMART[species], wanted), "biomart_human_ortholog"
    if species in GTF_MAP:
        return _gtf_id2sym(GTF_MAP[species]), "old_gtf_gene_name"
    raise KeyError(species)


def score_table(df: pd.DataFrame, id2sym: dict[str, str], genes: set[str]) -> dict[str, dict]:
    """Per-sample fraction and coverage. Duplicate IDs for one symbol → max TPM."""
    ids = df["id"].astype(str)
    sample_cols = [c for c in df.columns if c != "id"]
    mat = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib_tot = mat.sum(axis=0)
    mapped = ids.map(id2sym).fillna("")
    keep = mapped.isin(genes)
    sub = mat.loc[keep].copy()
    sub["_sym"] = mapped[keep].values
    collapsed = sub.groupby("_sym", observed=True)[sample_cols].max()
    out = {}
    for col in sample_cols:
        den = float(lib_tot[col])
        present = [g for g in genes if g in collapsed.index]
        num = float(collapsed.reindex(present)[col].fillna(0).sum()) if present else 0.0
        vals = collapsed.reindex(present)[col].fillna(0).astype(float) if present else pd.Series(dtype=float)
        out[col] = {
            "module_n_mapped": len(present),
            "module_n_expr": int((vals > 0).sum()) if len(vals) else 0,
            "frac": num / den if den else np.nan,
            "mean_log10": float(np.log10(vals.clip(lower=0) + 1).mean()) if len(vals) else np.nan,
            "lib_total": den,
        }
    return out


def fit_pair(y: pd.Series, masses: pd.Series) -> dict:
    ols = ols_logmass(y, masses)
    spr = spearman_logmass(y, masses)
    return {**{f"ols_{k}": v for k, v in ols.items()}, "spearman_rho": spr["rho"], "spearman_p": spr["p"]}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    defs = load_defs()
    wanted = defs["module"] | defs["set"]

    maps = {}
    map_how = {}
    for sp in FILES:
        maps[sp], map_how[sp] = id2symbol(sp, wanted)

    lib_rows = []
    cov_rows = []
    for sp, fn in FILES.items():
        df = pd.read_csv(FIB / fn, sep="\t")
        file_ids = set(df["id"].astype(str))
        mapped_ids = set(maps[sp]) & file_ids
        hit_mod = {maps[sp][i] for i in mapped_ids if maps[sp][i] in defs["module"]}
        hit_set = {maps[sp][i] for i in mapped_ids if maps[sp][i] in defs["set"]}
        cov_rows.append(
            dict(
                species=sp,
                study="GSE167580",
                map=map_how[sp],
                n_file_genes=len(file_ids),
                n_ids_mapped=len(mapped_ids),
                n_module=len(hit_mod),
                n_set=len(hit_set),
                missing_module=";".join(sorted(defs["module"] - hit_mod)),
            )
        )
        scored = {dname: score_table(df, maps[sp], genes) for dname, genes in defs.items()}
        for col in [c for c in df.columns if c != "id"]:
            klass = SAMPLE_CLASS.get((sp, col), "unknown")
            row = dict(
                study="GSE167580",
                species=sp,
                latin=LATIN[sp],
                sample=col,
                sample_class=klass,
                mass_g=MASS_G[sp],
                map=map_how[sp],
            )
            for dname in defs:
                s = scored[dname][col]
                row[f"{dname}_n_mapped"] = s["module_n_mapped"]
                row[f"{dname}_n_expr"] = s["module_n_expr"]
                row[f"{dname}_frac"] = s["frac"]
                row[f"{dname}_mean_log10"] = s["mean_log10"]
                row["lib_total"] = s["lib_total"]
            lib_rows.append(row)

    # GSE213253 seven birds, same fraction metric, official GEO tables
    for sp in SPECIES_ORDER:
        geo = load_geo_rnaseq(sp)
        gmat = geo["gene_mat"]
        tot = geo["lib_total"]
        key = sp.lower()
        for dname, genes in defs.items():
            present = sorted(set(gmat.index) & genes)
            for col in gmat.columns:
                den = float(tot[col])
                vals = gmat.reindex(present)[col].fillna(0).astype(float)
                lib_rows.append(
                    dict(
                        study="GSE213253",
                        species=key,
                        latin=LATIN[key],
                        sample=col,
                        sample_class="primary_like",
                        mass_g=MASS_G[key],
                        map=geo["mapped_how"],
                        **{
                            f"{dname}_n_mapped": len(present),
                            f"{dname}_n_expr": int((vals > 0).sum()),
                            f"{dname}_frac": float(vals.sum()) / den if den else np.nan,
                            f"{dname}_mean_log10": float(np.log10(vals.clip(lower=0) + 1).mean()) if present else np.nan,
                        },
                        lib_total=den,
                    )
                )
        cov_rows.append(
            dict(
                species=key,
                study="GSE213253",
                map=geo["mapped_how"],
                n_file_genes=geo["n_raw"],
                n_ids_mapped=geo["n_mapped_rows"],
                n_module=len(set(geo["gene_mat"].index) & defs["module"]),
                n_set=len(set(geo["gene_mat"].index) & defs["set"]),
                missing_module=";".join(sorted(defs["module"] - set(geo["gene_mat"].index))),
            )
        )

    lib = pd.DataFrame(lib_rows)
    # GSE213253 loop appended one row per definition per library; collapse
    # if both module_* and set_* landed on different rows.
    if lib.duplicated(["study", "species", "sample"]).any():
        keys = ["study", "species", "latin", "sample", "sample_class", "mass_g", "map", "lib_total"]
        parts = []
        for (study, species, sample), sub in lib.groupby(["study", "species", "sample"], sort=False):
            row = {k: sub.iloc[0][k] for k in keys if k in sub.columns}
            for col in sub.columns:
                if col in keys:
                    continue
                vals = sub[col].dropna()
                row[col] = vals.iloc[0] if len(vals) else np.nan
            parts.append(row)
        lib = pd.DataFrame(parts)

    lib.to_csv(OUT / "phenotype_by_library.tsv", sep="\t", index=False)
    pd.DataFrame(cov_rows).drop_duplicates(["study", "species"]).to_csv(OUT / "gene_coverage.tsv", sep="\t", index=False)

    # species means
    phen_rows = []
    cohorts = [
        ("GSE167580", "all_fib_labelled", lambda r: r.study.eq("GSE167580") & r.sample_class.ne("not_fibroblast")),
        ("GSE167580", "primary_like", lambda r: r.study.eq("GSE167580") & r.sample_class.eq("primary_like")),
        ("GSE167580", "primary_like_no_chicken", lambda r: r.study.eq("GSE167580") & r.sample_class.eq("primary_like") & r.species.ne("chicken")),
        ("GSE213253", "birds7", lambda r: r.study.eq("GSE213253")),
    ]
    for study, unit, filt in cohorts:
        sub = lib[filt(lib)]
        for sp, g in sub.groupby("species"):
            phen_rows.append(
                dict(
                    study=study,
                    unit=unit,
                    species=sp,
                    latin=g["latin"].iloc[0],
                    mass_g=float(g["mass_g"].iloc[0]),
                    n_samples=int(len(g)),
                    samples=";".join(g["sample"].astype(str)),
                    module_n_mapped=int(g["module_n_mapped"].median()),
                    set_n_mapped=int(g["set_n_mapped"].median()) if "set_n_mapped" in g else np.nan,
                    module_frac=float(g["module_frac"].mean()),
                    set_frac=float(g["set_frac"].mean()) if "set_frac" in g else np.nan,
                    module_mean_log10=float(g["module_mean_log10"].mean()),
                )
            )
    phen = pd.DataFrame(phen_rows)
    phen.to_csv(OUT / "phenotype_by_species.tsv", sep="\t", index=False)

    ols_rows = []
    for unit in phen["unit"].unique():
        sub = phen[phen.unit == unit].copy()
        masses = pd.Series(sub["mass_g"].values)
        for metric in ("module_frac", "set_frac", "module_mean_log10"):
            y = pd.Series(sub[metric].values)
            y.index = masses.index = range(len(y))
            fit = fit_pair(y, masses)
            ols_rows.append(dict(unit=unit, metric=metric, n_species=int(len(sub)), **fit))
        # per-library, same unit
        lib_sub = lib[lib.species.isin(set(sub.species)) & lib.study.eq(sub.study.iloc[0])]
        if unit == "primary_like":
            lib_sub = lib_sub[lib_sub.sample_class.eq("primary_like")]
        elif unit == "primary_like_no_chicken":
            lib_sub = lib_sub[lib_sub.sample_class.eq("primary_like") & lib_sub.species.ne("chicken")]
        elif unit == "all_fib_labelled":
            lib_sub = lib_sub[lib_sub.sample_class.ne("not_fibroblast")]
        for metric, col in (("module_frac", "module_frac"), ("set_frac", "set_frac")):
            y = pd.Series(lib_sub[col].values)
            xm = pd.Series(lib_sub["mass_g"].values)
            y.index = xm.index = range(len(y))
            fit = fit_pair(y, xm)
            ols_rows.append(dict(unit=unit + "_libraries", metric=metric, n_species=int(lib_sub.species.nunique()), **fit))
    ols = pd.DataFrame(ols_rows)
    ols.to_csv(OUT / "ols_spearman_vs_log10mass.tsv", sep="\t", index=False)

    # figure: two studies, primary-like GSE167580 + GSE213253 birds
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), constrained_layout=True)
    panels = [
        ("primary_like", "GSE167580 fibroblasts, primary-like", axes[0]),
        ("birds7", "GSE213253 embryonic fibroblasts (birds)", axes[1]),
    ]
    for unit, title, ax in panels:
        sub = phen[phen.unit == unit].sort_values("mass_g")
        x = np.log10(sub["mass_g"].values)
        y = sub["module_frac"].values
        ax.scatter(x, y, c="#1d3557", s=56, zorder=3)
        masses = pd.Series(sub["mass_g"].values)
        yy = pd.Series(y)
        fit = ols_logmass(yy, masses)
        if np.isfinite(fit["slope"]):
            xx = np.linspace(x.min(), x.max(), 40)
            ax.plot(xx, fit["intercept"] + fit["slope"] * xx, color="#1d3557", lw=1.4)
        for name, xi, yi in zip(sub["species"], x, y):
            ax.text(xi, yi, name[:4], fontsize=7)
        ax.set_title(f"{title}\nOLS r={fit['r']:.2f} p={fit['p']:.3g}  n={int(fit['n'])}")
        ax.set_xlabel("log10 adult mass (g)")
        ax.set_ylabel("frozen 29-gene module / library TPM")
    fig.suptitle("Fibroblast OXPHOS module fraction vs mass — studies not pooled", fontsize=11)
    fig.savefig(FIG / "fig_gse167580_fib_module_vs_mass.png", dpi=160)
    fig.savefig(OUT / "fig_gse167580_fib_module_vs_mass.png", dpi=160)
    plt.close()

    primary = ols[ols.unit.isin(["primary_like", "primary_like_no_chicken", "birds7"]) & ols.metric.eq("module_frac")]
    summary = {
        "date": "2026-08-28",
        "module": sorted(defs["module"]),
        "n_module": len(defs["module"]),
        "metric": "sum(mapped module TPM) / sum(all TPM); duplicate IDs per symbol take max",
        "studies_not_pooled": True,
        "dropped": {
            "rat_BRL_3A": "Buffalo rat liver epithelial line, not a fibroblast",
            "dog_A-72": "tumor fibroblasts / tumor line",
            "lines_kept_out_of_primary": ["WI38", "3T6", "PA12", "R9ab", "DF1"],
            "horse": "not in GSE167580",
        },
        "coverage": pd.DataFrame(cov_rows).drop_duplicates(["study", "species"]).to_dict(orient="records"),
        "primary_fits": primary.to_dict(orient="records"),
        "species_primary_like": phen[phen.unit == "primary_like"].sort_values("mass_g").to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(pd.DataFrame(cov_rows).drop_duplicates(["study", "species"]).to_string(index=False))
    print()
    print(phen[phen.unit.isin(["primary_like", "birds7"])].sort_values(["unit", "mass_g"]).to_string(index=False))
    print()
    print(ols.to_string(index=False))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
