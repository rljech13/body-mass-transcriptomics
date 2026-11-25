#!/usr/bin/env python3
"""Paper-A differentiation package (arms 1–4) + corpus-ready tables.

1) Shared / DE genes epi vs fib (n=7 Ruslan matrix) ∩ iES1300
2) iES GEMCAT compartment + consistent metabolite map across spp
3) OXPHOS module logFC(epi/fib) and OXPHOS-linked mets ~ mass (test vs −1/4)
4) CAGE promoter metrics for OXPHOS/mito genes (fib; chicken ESC gene-TPM)
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
DATA = ROOT / "data"
RES = ROOT / "results" / "paperA_diff_pack"
FIG = ROOT / "figures"
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")
GEM = ROOT / "results" / "gemcat_multispp_n7"
ATLAS_TPM = Path(
    "/path/to/user-home/chicken/CAGE_atlases/qc_atlas/results/05_gene_tpm_matrix.tsv"
)

MASS_G = {
    "chicken": 912.5,
    "duck": 1121.0,
    "emu": 36200.0,
    "peacock": 4093.75,
    "quail": 96.6,
    "turkey": 5811.0,
    "zebra": 11.93,
}
SPECIES_ORDER = ["zebra", "quail", "chicken", "duck", "peacock", "turkey", "emu"]

SAMPLE_COLS = [
    "Fibroblast.C07", "Fibroblast.C08", "Fibroblast.C09", "Fibroblast.C10",
    "Fibroblast.C11", "Fibroblast.C12", "PGC.C37", "PGC.C38",
    "Epiblast.C05", "Epiblast.C06", "Epiblast.C07",
    "Fibroblast.C25", "Fibroblast.C26", "Fibroblast.C27", "Fibroblast.C28",
    "Fibroblast.C29", "Fibroblast.C30",
    "Epiblast.C12", "Epiblast.C13", "Epiblast.C14",
    "Fibroblast.C13", "Fibroblast.C14", "Fibroblast.C15", "Fibroblast.C16",
    "Fibroblast.C17", "Fibroblast.C18",
    "Epiblast.C21", "Epiblast.C22", "Epiblast.C23",
    "Fibroblast.C01", "Fibroblast.C02", "Fibroblast.C03", "Fibroblast.C04",
    "Fibroblast.C05", "Fibroblast.C06",
    "Epiblast.C01", "Epiblast.C02", "Epiblast.C03", "Epiblast.C04",
    "Fibroblast.C31", "Fibroblast.C32", "Fibroblast.C33", "Fibroblast.C34",
    "Fibroblast.C35", "Fibroblast.C36",
    "Epiblast.C15", "Epiblast.C16", "Epiblast.C17",
    "Fibroblast.C39", "Fibroblast.C40", "Fibroblast.C41", "Fibroblast.C42",
    "Epiblast.C18", "Epiblast.C19", "Epiblast.C20",
    "Fibroblast.C19", "Fibroblast.C20", "Fibroblast.C21", "Fibroblast.C22",
    "Fibroblast.C23", "Fibroblast.C24",
    "Epiblast.C08", "Epiblast.C09", "Epiblast.C10", "Epiblast.C11",
]

COMP_NAME = {
    "c": "cytosol",
    "m": "mitochondria",
    "e": "extracellular",
    "r": "ER",
    "g": "golgi",
    "l": "lysosome",
    "n": "nucleus",
    "x": "peroxisome",
    "w": "cell wall/other",
    "q": "other_q",
}

CAGE_FIB = {
    "Chicken": ("GSE213253_Chicken_CAGE.xlsx", 2.0, "chicken"),
    "Duck": ("GSE213253_Duck_CAGE.xlsx", 1.1, "duck"),
    "Quail": ("GSE213253_Quail_CAGE.xlsx", 0.12, "quail"),
    "ZebraFinch": ("GSE213253_ZebraFinch_CAGE.xlsx", 0.0157, "zebra"),
    "Turkey": ("GSE213253_Turkey_CAGE.xlsx", 5.811, "turkey"),
    "Peacock": ("GSE213253_Peacock_CAGE.xlsx", 4.094, "peacock"),
    "Emu": ("GSE213253_Emu_CAGE.xlsx", 36.2, "emu"),
}


def load_modules() -> dict[str, list[str]]:
    out = {}
    for p in (DATA / "genesets").glob("*.txt"):
        out[p.stem] = [ln.strip().upper() for ln in p.read_text().splitlines() if ln.strip()]
    return out


def load_matrix() -> pd.DataFrame:
    genes = pd.read_csv(RUSLAN / "cross-birds/counts.csv", index_col=0)
    # counts.csv: index=0..n, single column of symbols OR two cols
    if genes.shape[1] == 1:
        symbols = genes.iloc[:, 0].astype(str).str.upper()
    else:
        symbols = genes.iloc[:, 0].astype(str).str.upper()
    exp = pd.read_csv(RUSLAN / "cross-birds/file.csv", index_col=0)
    exp.columns = SAMPLE_COLS
    exp.index = symbols.values
    # drop duplicate symbols keeping first
    exp = exp[~exp.index.duplicated(keep="first")]
    return exp


def species_of(col: str) -> str:
    # reverse from SAMPLE_COLS blocks
    idx = SAMPLE_COLS.index(col)
    blocks = [
        (0, 11, "chicken"),
        (11, 20, "quail"),
        (20, 29, "turkey"),
        (29, 39, "peacock"),
        (39, 48, "duck"),
        (48, 55, "zebra"),
        (55, 65, "emu"),
    ]
    for a, b, sp in blocks:
        if a <= idx < b:
            return sp
    raise KeyError(col)


def arm1_genes(exp: pd.DataFrame, modules: dict[str, list[str]]) -> dict:
    sym2 = pd.read_csv(DATA / "symbol_to_chicken_entrez.tsv", sep="\t")
    sym2["symbol"] = sym2["symbol"].astype(str).str.upper()
    ies_syms = set(sym2["symbol"])

    per = []
    for sp in SPECIES_ORDER:
        fib = [c for c in exp.columns if c.startswith("Fibroblast") and species_of(c) == sp]
        epi = [c for c in exp.columns if c.startswith("Epiblast") and species_of(c) == sp]
        fib_m = exp[fib].mean(axis=1)
        epi_m = exp[epi].mean(axis=1)
        # matrix already log-like; logFC ≈ epi - fib
        logfc = epi_m - fib_m
        shared_high = (fib_m >= fib_m.quantile(0.75)) & (epi_m >= epi_m.quantile(0.75))
        df = pd.DataFrame(
            {
                "gene": exp.index,
                "species": sp,
                "fib_mean": fib_m.values,
                "epi_mean": epi_m.values,
                "logFC_epi_minus_fib": logfc.values,
                "shared_high_q75": shared_high.values,
                "in_ies1300_map": [g in ies_syms for g in exp.index],
            }
        )
        per.append(df)
    long = pd.concat(per, ignore_index=True)
    long.to_csv(RES / "arm1_gene_epi_fib_per_species.tsv", sep="\t", index=False)

    # wide logFC
    wide = long.pivot(index="gene", columns="species", values="logFC_epi_minus_fib")
    wide["n_up_epi"] = (wide[SPECIES_ORDER] > 0.5).sum(axis=1)
    wide["n_down_epi"] = (wide[SPECIES_ORDER] < -0.5).sum(axis=1)
    wide["n_stable"] = ((wide[SPECIES_ORDER].abs() <= 0.25).sum(axis=1))
    wide["mean_logFC"] = wide[SPECIES_ORDER].mean(axis=1)
    wide["in_ies1300_map"] = wide.index.isin(ies_syms)
    for mod, genes in modules.items():
        wide[f"module_{mod}"] = wide.index.isin(set(genes))
    wide = wide.reset_index()
    wide.to_csv(RES / "arm1_gene_logFC_wide_n7.tsv", sep="\t", index=False)

    cons_up = wide[wide["n_up_epi"] >= 5].sort_values("mean_logFC", ascending=False)
    cons_down = wide[wide["n_down_epi"] >= 5].sort_values("mean_logFC")
    cons_stable = wide[(wide["n_stable"] >= 5) & (wide["mean_logFC"].abs() <= 0.25)]
    shared_all = (
        long.groupby("gene")["shared_high_q75"].sum().reset_index(name="n_spp_shared_high")
    )
    shared_all = shared_all[shared_all["n_spp_shared_high"] >= 5]
    shared_all = shared_all.merge(
        wide[["gene", "mean_logFC", "in_ies1300_map"] + [f"module_{m}" for m in modules]],
        on="gene",
        how="left",
    )
    cons_up.to_csv(RES / "arm1_consistent_up_epi_ge5spp.tsv", sep="\t", index=False)
    cons_down.to_csv(RES / "arm1_consistent_down_epi_ge5spp.tsv", sep="\t", index=False)
    cons_stable.to_csv(RES / "arm1_stable_ge5spp.tsv", sep="\t", index=False)
    shared_all.to_csv(RES / "arm1_shared_high_ge5spp.tsv", sep="\t", index=False)

    summary = {
        "n_genes_matrix": int(exp.shape[0]),
        "n_in_ies_map": int(wide["in_ies1300_map"].sum()),
        "n_cons_up_epi_ge5": int(len(cons_up)),
        "n_cons_down_epi_ge5": int(len(cons_down)),
        "n_stable_ge5": int(len(cons_stable)),
        "n_shared_high_ge5": int(len(shared_all)),
        "n_cons_up_in_ies": int(cons_up["in_ies1300_map"].sum()),
        "n_cons_down_in_ies": int(cons_down["in_ies1300_map"].sum()),
        "module_coverage_in_matrix": {
            m: int(wide[f"module_{m}"].sum()) for m in modules
        },
    }
    (RES / "arm1_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def arm2_gemcat_map() -> dict:
    rows = []
    met_long = []
    for sp in SPECIES_ORDER:
        df = pd.read_csv(GEM / f"{sp}_gemcat_ies1300_sorted.tsv", sep="\t")
        df = df[pd.to_numeric(df["score"], errors="coerce").fillna(0) > 0].copy()
        df["log2"] = np.log2(df["score"].astype(float))
        df["comp"] = df["metabolite"].astype(str).str.extract(r"\[([^\]]+)\]$")[0]
        df["species"] = sp
        met_long.append(df)
        for comp, sub in df.groupby("comp"):
            rows.append(
                {
                    "species": sp,
                    "compartment": comp,
                    "compartment_name": COMP_NAME.get(comp, comp),
                    "n_mets": len(sub),
                    "median_log2": float(sub["log2"].median()),
                    "mean_log2": float(sub["log2"].mean()),
                    "frac_up": float((sub["log2"] > 0).mean()),
                    "log10_mass_g": float(np.log10(MASS_G[sp])),
                }
            )
    comp_tab = pd.DataFrame(rows)
    comp_tab.to_csv(RES / "arm2_ies_compartment_medians.tsv", sep="\t", index=False)

    # OLS compartment ~ mass
    fits = []
    for comp, sub in comp_tab.groupby("compartment"):
        if len(sub) < 5:
            continue
        sl, icept, r, p, se = stats.linregress(sub["log10_mass_g"], sub["median_log2"])
        fits.append(
            {
                "compartment": comp,
                "compartment_name": COMP_NAME.get(comp, comp),
                "slope": sl,
                "intercept": icept,
                "r": r,
                "p": p,
                "se": se,
                "mean_median_log2": float(sub["median_log2"].mean()),
            }
        )
    fits_df = pd.DataFrame(fits).sort_values("p")
    fits_df.to_csv(RES / "arm2_ies_compartment_vs_mass_ols.tsv", sep="\t", index=False)

    big = pd.concat(met_long, ignore_index=True)
    # consistent direction across spp
    piv = big.pivot_table(index="metabolite", columns="species", values="log2", aggfunc="first")
    piv = piv.reindex(columns=SPECIES_ORDER)
    piv["n_up"] = (piv[SPECIES_ORDER] > np.log2(1.5)).sum(axis=1)
    piv["n_down"] = (piv[SPECIES_ORDER] < np.log2(1 / 1.5)).sum(axis=1)
    piv["mean_log2"] = piv[SPECIES_ORDER].mean(axis=1)
    piv["comp"] = piv.index.to_series().astype(str).str.extract(r"\[([^\]]+)\]$")[0].values
    cons_up = piv[piv["n_up"] >= 5].sort_values("mean_log2", ascending=False)
    cons_down = piv[piv["n_down"] >= 5].sort_values("mean_log2")
    cons_up.to_csv(RES / "arm2_consistent_up_mets_ge5spp.tsv", sep="\t")
    cons_down.to_csv(RES / "arm2_consistent_down_mets_ge5spp.tsv", sep="\t")

    # figure compartments
    fig, ax = plt.subplots(figsize=(8, 4.2), constrained_layout=True)
    for comp, sub in comp_tab.groupby("compartment"):
        if comp not in {"c", "m", "e", "r", "g", "n"}:
            continue
        sub = sub.sort_values("log10_mass_g")
        ax.plot(
            sub["log10_mass_g"],
            sub["median_log2"],
            marker="o",
            label=COMP_NAME.get(comp, comp),
        )
    ax.axhline(0, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("log10 adult body mass (g)")
    ax.set_ylabel("median log2 GEMCAT (epi/fib FC)")
    ax.set_title("iES1300 compartments across species")
    ax.legend(fontsize=7, ncol=2)
    fig.savefig(FIG / "fig_arm2_ies_compartments_vs_mass.png", dpi=160)
    plt.close(fig)

    summary = {
        "n_consistent_up_mets_ge5": int(len(cons_up)),
        "n_consistent_down_mets_ge5": int(len(cons_down)),
        "mito_mean_median_log2": float(
            comp_tab.loc[comp_tab.compartment == "m", "median_log2"].mean()
        ),
        "top_comp_by_p": fits_df.head(3).to_dict(orient="records"),
    }
    (RES / "arm2_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def oxphos_mets_from_ies() -> set[str]:
    """Metabolites participating in reactions with OXPHOS-like genes/subsystems."""
    import cobra

    model = cobra.io.read_sbml_model(
        str(DATA / "ies1300/S1_File/S1_File in sbml format.xml")
    )
    ox_syms = set(load_modules()["OXPHOS"])
    # map symbols → chicken entrez present in model
    sym2 = pd.read_csv(DATA / "symbol_to_chicken_entrez.tsv", sep="\t")
    sym2["symbol"] = sym2["symbol"].astype(str).str.upper()
    ox_entrez = set(
        sym2.loc[sym2["symbol"].isin(ox_syms), "chicken_entrez"].astype(int).astype(str)
    )
    mets = set()
    for rxn in model.reactions:
        gids = {g.id for g in rxn.genes}
        name = (rxn.name or "").lower() + " " + rxn.id.lower()
        hit = bool(gids & ox_entrez) or any(
            k in name
            for k in (
                "oxidative phosphorylation",
                "complex i",
                "complex ii",
                "complex iii",
                "complex iv",
                "atp synthase",
                "cytochrome",
                "ubiquinol",
            )
        )
        if hit:
            for m in rxn.metabolites:
                mets.add(m.id)
    return mets


def arm3_oxphos(exp: pd.DataFrame, modules: dict[str, list[str]]) -> dict:
    ox = set(modules["OXPHOS"])
    present = [g for g in ox if g in exp.index]
    rows = []
    gene_fc_rows = []
    for sp in SPECIES_ORDER:
        fib = [c for c in exp.columns if c.startswith("Fibroblast") and species_of(c) == sp]
        epi = [c for c in exp.columns if c.startswith("Epiblast") and species_of(c) == sp]
        fib_s = exp.loc[present, fib].mean(axis=1)
        epi_s = exp.loc[present, epi].mean(axis=1)
        # module score = mean across genes; logFC module = mean(epi)-mean(fib)
        fib_mod = float(fib_s.mean())
        epi_mod = float(epi_s.mean())
        logfc = epi_mod - fib_mod
        # also gene-wise mean logFC
        gfc = float((epi_s - fib_s).mean())
        rows.append(
            {
                "species": sp,
                "mass_g": MASS_G[sp],
                "log10_mass_g": np.log10(MASS_G[sp]),
                "n_oxphos_genes": len(present),
                "fib_module": fib_mod,
                "epi_module": epi_mod,
                "logFC_epi_minus_fib": logfc,
                "mean_gene_logFC": gfc,
            }
        )
        for g in present:
            gene_fc_rows.append(
                {
                    "species": sp,
                    "gene": g,
                    "logFC": float(epi_s[g] - fib_s[g]),
                    "log10_mass_g": np.log10(MASS_G[sp]),
                }
            )
    mod = pd.DataFrame(rows)
    mod.to_csv(RES / "arm3_oxphos_module_epi_fib.tsv", sep="\t", index=False)
    pd.DataFrame(gene_fc_rows).to_csv(
        RES / "arm3_oxphos_gene_logFC_per_species.tsv", sep="\t", index=False
    )

    fits = []
    for col in ["logFC_epi_minus_fib", "fib_module", "epi_module", "mean_gene_logFC"]:
        sl, icept, r, p, se = stats.linregress(mod["log10_mass_g"], mod[col])
        # distance to cellular Kleiber targets for logFC~mass only
        err_m14 = abs(sl - (-0.25)) if "logFC" in col or col == "mean_gene_logFC" else np.nan
        err_m13 = abs(sl - (-1 / 3)) if "logFC" in col or col == "mean_gene_logFC" else np.nan
        fits.append(
            {
                "endpoint": col,
                "slope": sl,
                "intercept": icept,
                "r": r,
                "p": p,
                "se": se,
                "abs_err_vs_m14": err_m14,
                "abs_err_vs_m13": err_m13,
            }
        )

    # OXPHOS-linked metabolite median from iES GEMCAT
    ox_mets = oxphos_mets_from_ies()
    met_rows = []
    for sp in SPECIES_ORDER:
        df = pd.read_csv(GEM / f"{sp}_gemcat_ies1300_sorted.tsv", sep="\t")
        df = df[df["score"] > 0].copy()
        df["log2"] = np.log2(df["score"].astype(float))
        sub = df[df["metabolite"].isin(ox_mets)]
        mito_ox = sub[sub["metabolite"].astype(str).str.endswith("[m]")]
        met_rows.append(
            {
                "species": sp,
                "log10_mass_g": np.log10(MASS_G[sp]),
                "n_oxphos_mets": len(sub),
                "median_log2_oxphos_mets": float(sub["log2"].median()) if len(sub) else np.nan,
                "n_oxphos_mito_mets": len(mito_ox),
                "median_log2_oxphos_mito": float(mito_ox["log2"].median())
                if len(mito_ox)
                else np.nan,
            }
        )
    met = pd.DataFrame(met_rows)
    met.to_csv(RES / "arm3_oxphos_mets_gemcat.tsv", sep="\t", index=False)
    for col in ["median_log2_oxphos_mets", "median_log2_oxphos_mito"]:
        y = met[col].values
        if np.all(np.isfinite(y)):
            sl, icept, r, p, se = stats.linregress(met["log10_mass_g"], y)
            fits.append(
                {
                    "endpoint": col,
                    "slope": sl,
                    "intercept": icept,
                    "r": r,
                    "p": p,
                    "se": se,
                    "abs_err_vs_m14": abs(sl - (-0.25)),
                    "abs_err_vs_m13": abs(sl - (-1 / 3)),
                }
            )
    fits_df = pd.DataFrame(fits)
    fits_df.to_csv(RES / "arm3_oxphos_vs_mass_ols.tsv", sep="\t", index=False)

    # figure
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), constrained_layout=True)
    for ax, col, title in [
        (axes[0], "logFC_epi_minus_fib", "OXPHOS module logFC (epi−fib)"),
        (axes[1], "fib_module", "OXPHOS absolute fib"),
        (axes[2], "epi_module", "OXPHOS absolute epi"),
    ]:
        ax.scatter(mod["log10_mass_g"], mod[col], s=50)
        for _, r in mod.iterrows():
            ax.text(r["log10_mass_g"], r[col], r["species"], fontsize=7)
        sl, icept, _, p, _ = stats.linregress(mod["log10_mass_g"], mod[col])
        xs = np.linspace(mod["log10_mass_g"].min(), mod["log10_mass_g"].max(), 40)
        ax.plot(xs, icept + sl * xs, color="C1")
        if "logFC" in col:
            ax.plot(xs, icept + (-0.25) * (xs - xs.mean()) + mod[col].mean(), "k--", lw=0.8, label="slope −1/4 ref")
            ax.legend(fontsize=7)
        ax.set_title(f"{title}\nslope={sl:.3f}, p={p:.3f}")
        ax.set_xlabel("log10 mass (g)")
    fig.suptitle("Arm3: OXPHOS vs mass (n=7)")
    fig.savefig(FIG / "fig_arm3_oxphos_vs_mass.png", dpi=160)
    plt.close(fig)

    summary = {
        "n_oxphos_in_matrix": len(present),
        "n_oxphos_mets_ies": len(ox_mets),
        "fits": fits_df.to_dict(orient="records"),
    }
    (RES / "arm3_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def promoter_stats(path: Path, genes: set[str]) -> pd.DataFrame:
    df = pd.read_excel(path)
    if "geneId" not in df.columns:
        return pd.DataFrame()
    norms = [c for c in df.columns if str(c).startswith("norm.")]
    if not norms:
        return pd.DataFrame()
    df["geneId"] = df["geneId"].astype(str).str.upper()
    sub = df[df["geneId"].isin(genes)].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["_expr"] = sub[norms].apply(pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
    rows = []
    for g, gdf in sub.groupby("geneId"):
        vals = gdf["_expr"].astype(float).values
        pos = vals[vals > 0]
        if len(pos) == 0:
            continue
        p = pos / pos.sum()
        shannon = float(-(p * np.log2(p + 1e-12)).sum())
        rows.append(
            {
                "gene": g,
                "n_promoters": int(len(gdf)),
                "n_active_promoters": int((gdf["_expr"] > 0).sum()),
                "total_tpm": float(gdf["_expr"].sum()),
                "dominant_frac": float(p.max()),
                "shannon": shannon,
            }
        )
    return pd.DataFrame(rows)


def arm4_cage(modules: dict[str, list[str]]) -> dict:
    oxmito = set(modules["OXPHOS"]) | set(modules["MITO_BIOGENESIS"]) | set(modules["TCA"])
    all_rows = []
    coverage = []
    for sp, (fname, mass_kg, sp_key) in CAGE_FIB.items():
        path = DATA / fname
        if not path.exists():
            continue
        st = promoter_stats(path, oxmito)
        n_hit = 0 if st.empty else st["gene"].nunique()
        coverage.append(
            {
                "species": sp,
                "sp_key": sp_key,
                "mass_kg": mass_kg,
                "n_oxmito_genes_with_promoter": n_hit,
                "usable_symbols": n_hit >= 10,
            }
        )
        if st.empty or n_hit < 5:
            continue
        st["species"] = sp
        st["sp_key"] = sp_key
        st["mass_kg"] = mass_kg
        st["module"] = st["gene"].map(
            lambda g: (
                "OXPHOS"
                if g in modules["OXPHOS"]
                else "MITO_BIOGENESIS"
                if g in modules["MITO_BIOGENESIS"]
                else "TCA"
            )
        )
        all_rows.append(st)
    cov = pd.DataFrame(coverage)
    cov.to_csv(RES / "arm4_cage_species_coverage.tsv", sep="\t", index=False)

    summary = {"n_species_usable_fib_cage": 0}
    if all_rows:
        big = pd.concat(all_rows, ignore_index=True)
        big.to_csv(RES / "arm4_cage_promoter_oxmito_genes.tsv", sep="\t", index=False)
        # overall + OXPHOS-only summaries
        for label, sub in [
            ("oxmito", big),
            ("OXPHOS", big[big["module"] == "OXPHOS"]),
        ]:
            summ = (
                sub.groupby("species")
                .agg(
                    mass_kg=("mass_kg", "first"),
                    n_genes=("gene", "nunique"),
                    mean_n_active=("n_active_promoters", "mean"),
                    mean_shannon=("shannon", "mean"),
                    mean_dominant_frac=("dominant_frac", "mean"),
                    mean_total_tpm=("total_tpm", "mean"),
                )
                .reset_index()
            )
            summ.to_csv(RES / f"arm4_cage_{label}_species_summary.tsv", sep="\t", index=False)
            fit_rows = []
            for col in [
                "mean_n_active",
                "mean_shannon",
                "mean_dominant_frac",
                "mean_total_tpm",
            ]:
                if len(summ) >= 3:
                    sl, icept, r, p, se = stats.linregress(
                        np.log10(summ["mass_kg"]), summ[col]
                    )
                    fit_rows.append(
                        {"module_set": label, "metric": col, "slope": sl, "r": r, "p": p, "se": se}
                    )
            pd.DataFrame(fit_rows).to_csv(
                RES / f"arm4_cage_{label}_vs_mass_ols.tsv", sep="\t", index=False
            )
        summary["n_species_usable_fib_cage"] = int(big["species"].nunique())

        fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), constrained_layout=True)
        summ = pd.read_csv(RES / "arm4_cage_OXPHOS_species_summary.tsv", sep="\t")
        for ax, col, title in zip(
            axes,
            ["mean_n_active", "mean_shannon", "mean_total_tpm"],
            ["Active promoters", "Shannon", "Total TPM"],
        ):
            x = np.log10(summ["mass_kg"])
            y = summ[col]
            ax.scatter(x, y, s=50)
            for sp, xi, yi in zip(summ["species"], x, y):
                ax.text(xi, yi, sp, fontsize=7)
            if len(summ) >= 3:
                sl, icept, _, p, _ = stats.linregress(x, y)
                xs = np.linspace(x.min(), x.max(), 40)
                ax.plot(xs, icept + sl * xs, color="C1")
                ax.set_title(f"{title}\nslope={sl:.3f}, p={p:.3f}")
            ax.set_xlabel("log10 mass (kg)")
        fig.suptitle("Arm4: CAGE OXPHOS promoters (fib)")
        fig.savefig(FIG / "fig_arm4_cage_oxphos_promoters_vs_mass.png", dpi=160)
        plt.close(fig)

    # Chicken ESC vs fib gene-level TPM (atlas) for OXPHOS
    if ATLAS_TPM.exists():
        tpm = pd.read_csv(ATLAS_TPM, sep="\t")
        tpm["gene_name"] = tpm["gene_name"].astype(str).str.upper()
        # chicken fib C07-C12, ESC often E05-E07 or similar — detect
        cols = list(tpm.columns)
        fib_cols = [c for c in cols if c.startswith("C0") and any(f"C0{i}" in c or f"C{i}" in c for i in range(7, 13))]
        # more robust: C07..C12 in name
        fib_cols = [c for c in cols if any(x in c for x in ["C07", "C08", "C09", "C10", "C11", "C12"])]
        esc_cols = [
            c
            for c in cols
            if any(x in c for x in ["E05", "E06", "E07", "ESC", "esc"])
            or c.startswith("E0")
        ]
        # Lizio / atlas: look for known patterns from earlier work
        meta_path = Path(
            "/path/to/user-home/chicken/CAGE_atlases/meta/atlas_metadata_unified.tsv"
        )
        esc_note = {"fib_cols": fib_cols[:20], "esc_cols": esc_cols[:20]}
        if meta_path.exists():
            meta = pd.read_csv(meta_path, sep="\t")
            esc_note["meta_cols"] = meta.columns.tolist()[:30]
            # try find ESC/fib sample ids
            for col in meta.columns:
                if meta[col].dtype == object:
                    vals = meta[col].astype(str).str.lower()
                    if vals.str.contains("fibroblast").any() or vals.str.contains("esc").any():
                        esc_note["annot_col"] = col
                        break
        # Use columns from chicken_esc_vs_fib if available
        prev = ROOT / "results" / "chicken_esc_vs_fib_gene_tpm.tsv"
        if prev.exists():
            prev_df = pd.read_csv(prev, sep="\t")
            prev_df["gene"] = prev_df.iloc[:, 0].astype(str).str.upper()
            # expect fib_mean, esc_mean or similar
            cols_l = {c.lower(): c for c in prev_df.columns}
            fc_col = None
            for cand in ["logfc", "log2fc", "delta", "logfc_esc_minus_fib"]:
                if cand in cols_l:
                    fc_col = cols_l[cand]
                    break
            ox = prev_df[prev_df["gene"].isin(oxmito)].copy()
            ox.to_csv(RES / "arm4_chicken_atlas_oxmito_esc_fib.tsv", sep="\t", index=False)
            summary["chicken_atlas_oxmito_n"] = int(len(ox))
            summary["chicken_atlas_columns"] = prev_df.columns.tolist()
            if fc_col:
                summary["chicken_atlas_oxphos_mean_logFC"] = float(
                    ox.loc[ox["gene"].isin(modules["OXPHOS"]), fc_col].mean()
                )
        summary["atlas_col_probe"] = esc_note

    (RES / "arm4_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    return summary


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    modules = load_modules()
    print("Loading Ruslan matrix...", flush=True)
    exp = load_matrix()
    print(f"matrix {exp.shape}", flush=True)

    print("=== ARM1 genes ===", flush=True)
    s1 = arm1_genes(exp, modules)
    print(json.dumps(s1, indent=2), flush=True)

    print("=== ARM2 GEMCAT map ===", flush=True)
    s2 = arm2_gemcat_map()
    print(json.dumps(s2, indent=2), flush=True)

    print("=== ARM3 OXPHOS ===", flush=True)
    s3 = arm3_oxphos(exp, modules)
    print(json.dumps({k: v for k, v in s3.items() if k != "fits"}, indent=2), flush=True)
    print(pd.DataFrame(s3["fits"]).to_string(index=False), flush=True)

    print("=== ARM4 CAGE ===", flush=True)
    s4 = arm4_cage(modules)
    print(json.dumps(s4, indent=2, default=str)[:2000], flush=True)

    # combined claims stub for corpus
    claims = {
        "date": "2026-08-06",
        "arms": {"1": s1, "2": s2, "3": s3, "4": s4},
        "paths": {
            "results": str(RES),
            "figures": [
                "fig_arm2_ies_compartments_vs_mass.png",
                "fig_arm3_oxphos_vs_mass.png",
                "fig_arm4_cage_oxphos_promoters_vs_mass.png",
            ],
        },
    }
    (RES / "paperA_diff_pack_summary.json").write_text(json.dumps(claims, indent=2, default=str))
    print("DONE", RES, flush=True)


if __name__ == "__main__":
    main()
