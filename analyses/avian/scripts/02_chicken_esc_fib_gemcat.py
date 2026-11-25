#!/usr/bin/env python3
"""A4–A5: chicken ESC vs fibroblast metabolic modules + GEMCAT on fold-change."""
from __future__ import annotations

import json
import os
import re
import subprocess
import warnings
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
DATA = ROOT / "data"
RES = ROOT / "results"
FIG = ROOT / "figures"
ATLAS_TPM = Path("/path/to/user-home/chicken/CAGE_atlases/qc_atlas/results/05_gene_tpm_matrix.tsv")
RECON = DATA / "Recon3D.xml"
GEMCAT_ENV = {"LD_PRELOAD": "/usr/lib/x86_64-linux-gnu/libsqlite3.so.0", **os.environ}

FIB_COLS = [
    "C08_S8_R1_001.tpm",
    "C09_S9_R1_001.tpm",
    "C10_S32_R1_001.tpm",
    "C11_S10_R1_001.tpm",
    "C12_S11_R1_001.tpm",
]
ESC_COLS = [
    "E05_S13_R1_001.tpm",
    "E06_S14_R1_001.tpm",
    "E07_S15_R1_001.tpm",
]


def load_modules() -> dict[str, list[str]]:
    out = {}
    for p in (DATA / "genesets").glob("*.txt"):
        out[p.stem] = [ln.strip().upper() for ln in p.read_text().splitlines() if ln.strip()]
    return out


def module_score(expr: pd.Series, genes: list[str], min_genes: int = 5) -> tuple[float, int]:
    present = [g for g in genes if g in expr.index]
    if len(present) < min_genes:
        return np.nan, len(present)
    return float(np.log2(expr.loc[present].astype(float) + 1.0).mean()), len(present)


def symbol_to_human_entrez(symbols: list[str]) -> dict[str, str]:
    """Map chicken/human-shared symbols to human Entrez via mygene (species=human)."""
    import mygene

    mg = mygene.MyGeneInfo()
    mapping = {}
    for i in range(0, len(symbols), 800):
        batch = symbols[i : i + 800]
        print(f"  mygene batch {i}-{i+len(batch)} / {len(symbols)}")
        hits = mg.querymany(
            batch,
            scopes="symbol",
            fields="entrezgene",
            species=9606,
            returnall=False,
            as_dataframe=False,
            verbose=False,
        )
        for h in hits:
            q = str(h.get("query", "")).upper()
            if "entrezgene" in h and not h.get("notfound"):
                mapping[q] = str(int(h["entrezgene"]))
    return mapping


def expand_entrez_to_recon_ids(entrez_fc: pd.DataFrame, model_path: Path) -> pd.DataFrame:
    import cobra

    model = cobra.io.read_sbml_model(str(model_path))
    entrez_to_at = defaultdict(list)
    for g in model.genes:
        m = re.match(r"^(\d+)_AT\d+$", g.id)
        if m:
            entrez_to_at[m.group(1)].append(g.id)
    rows = []
    for entrez, val in zip(entrez_fc["gene"].astype(str), entrez_fc["logFC"]):
        for gid in entrez_to_at.get(entrez, []):
            rows.append((gid, val))
    return pd.DataFrame(rows, columns=["gene", "logFC"])


def run_gemcat(de_path: Path, out_path: Path) -> None:
    cmd = ["gemcat", str(RECON), str(de_path), "-e", "logFC", "-o", str(out_path)]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, env=GEMCAT_ENV)


def summarize_gemcat(result_csv: Path, top_n: int = 30) -> pd.DataFrame:
    df = pd.read_csv(result_csv, header=None)
    # gemcat output: metabolite, score (sometimes with header)
    if df.shape[1] >= 2:
        df.columns = ["metabolite", "score"] + [f"c{i}" for i in range(df.shape[1] - 2)]
    else:
        df.columns = ["metabolite"]
        df["score"] = np.nan
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna(subset=["score"]).sort_values("score", ascending=False)
    return df


def pathway_proxy_from_gene_fc(gene_fc: pd.Series, modules: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for mod, genes in modules.items():
        present = [g for g in genes if g in gene_fc.index]
        if len(present) < 5:
            rows.append({"module": mod, "n": len(present), "mean_logFC": np.nan, "t_p": np.nan})
            continue
        vals = gene_fc.loc[present].astype(float)
        t, p = stats.ttest_1samp(vals, 0.0)
        rows.append({"module": mod, "n": len(present), "mean_logFC": float(vals.mean()), "t_p": float(p)})
    return pd.DataFrame(rows)


def main():
    modules = load_modules()
    print("Loading atlas gene TPM...")
    tpm = pd.read_csv(ATLAS_TPM, sep="\t", index_col=0)
    tpm.index = tpm.index.astype(str).str.upper()
    for c in FIB_COLS + ESC_COLS:
        if c not in tpm.columns:
            raise SystemExit(f"Missing column {c} in atlas TPM")

    fib = tpm[FIB_COLS].mean(axis=1)
    esc = tpm[ESC_COLS].mean(axis=1)
    # logFC ESC/fib (pseudocount)
    logfc = np.log2((esc + 1.0) / (fib + 1.0))
    logfc.name = "logFC_ESC_over_fib"

    pd.DataFrame({"fib_tpm_mean": fib, "esc_tpm_mean": esc, "logFC_ESC_over_fib": logfc}).to_csv(
        RES / "chicken_esc_vs_fib_gene_tpm.tsv", sep="\t"
    )

    # module scores
    rows = []
    for mod, genes in modules.items():
        sf, nf = module_score(fib, genes)
        se, ne = module_score(esc, genes)
        rows.append(
            {
                "module": mod,
                "fib_score": sf,
                "esc_score": se,
                "delta_esc_minus_fib": se - sf if np.isfinite(se) and np.isfinite(sf) else np.nan,
                "n_fib": nf,
                "n_esc": ne,
            }
        )
    mod_tab = pd.DataFrame(rows)
    mod_tab.to_csv(RES / "chicken_esc_fib_module_scores.tsv", sep="\t", index=False)
    print(mod_tab.to_string(index=False))

    # pathway proxy from gene FC
    path_proxy = pathway_proxy_from_gene_fc(logfc, modules)
    path_proxy.to_csv(RES / "chicken_esc_fib_module_logFC_ttest.tsv", sep="\t", index=False)

    # overlap with fib allometry gene slopes (A2/A3)
    slope_path = RES / "per_gene_logtpm_vs_logmass_slopes.tsv"
    overlap_note = {"available": False}
    if slope_path.exists():
        slopes = pd.read_csv(slope_path, sep="\t")
        # genes with negative allometry slope (mass-up → expr-down)
        neg = set(slopes.loc[slopes.slope < 0, "gene"].astype(str).str.upper())
        # genes downregulated in ESC vs fib
        down_esc = set(logfc[logfc < -0.5].index)
        up_esc = set(logfc[logfc > 0.5].index)
        overlap_note = {
            "available": True,
            "n_neg_allometry": len(neg),
            "n_down_in_esc": len(down_esc & set(slopes.gene.str.upper())),
            "overlap_neg_allometry_and_down_esc": sorted(neg & down_esc)[:100],
            "n_overlap_neg_down": len(neg & down_esc),
            "overlap_neg_allometry_and_up_esc": sorted(neg & up_esc)[:100],
            "n_overlap_neg_up": len(neg & up_esc),
            "interpretation": "Q4: do fib mass-scaling genes also change at ESC differentiation? Intersection sizes reported; multi-spp ESC unavailable.",
        }
        with open(RES / "A4_gene_overlap_allometry_vs_escfib.json", "w") as f:
            json.dump(overlap_note, f, indent=2)

    # figure module delta
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    m = mod_tab.dropna(subset=["delta_esc_minus_fib"])
    colors = ["#b33" if v < 0 else "#33b" for v in m.delta_esc_minus_fib]
    ax.barh(m.module, m.delta_esc_minus_fib, color=colors)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("Δ module score (ESC − fibroblast)")
    ax.set_title("Chicken CAGE atlas: metabolic modules ESC vs fibroblasts")
    fig.savefig(FIG / "fig_A4_chicken_esc_fib_modules.png", dpi=160)
    plt.close(fig)

    # GEMCAT: map symbols → human Entrez → Recon3D gene IDs
    # Restrict to informative DE genes to keep mygene/GEMCAT tractable
    keep = (fib + esc) > 1.0
    de = logfc[keep].replace([np.inf, -np.inf], np.nan).dropna()
    # prioritize metabolic modules + strongest |logFC|
    meta_genes = set()
    for gs in modules.values():
        meta_genes.update(gs)
    strong = set(de.index[de.abs() >= 0.75].astype(str))
    symbols = sorted((set(de.index.astype(str)) & meta_genes) | strong)
    # hard cap
    if len(symbols) > 2500:
        top = de.abs().sort_values(ascending=False).head(2000).index.astype(str)
        symbols = sorted(set(top) | (set(de.index.astype(str)) & meta_genes))
    print(f"Mapping {len(symbols)} symbols to human Entrez...")
    sym2ent = symbol_to_human_entrez(symbols)
    mapped = [(sym2ent[s], float(de.loc[s])) for s in symbols if s in sym2ent]
    entrez_fc = pd.DataFrame(mapped, columns=["gene", "logFC"]).drop_duplicates("gene", keep="first")
    entrez_fc.to_csv(RES / "chicken_esc_fib_logFC_human_entrez.csv", index=False)
    print(f"Mapped to Entrez: {len(entrez_fc)} / {len(symbols)}")

    print("Expanding to Recon3D _AT ids...")
    fixed = expand_entrez_to_recon_ids(entrez_fc, RECON)
    fixed_path = RES / "chicken_esc_fib_logFC_recon3d_fixedID.csv"
    fixed.to_csv(fixed_path, index=False)
    print(f"Recon3D gene rows: {len(fixed)}")

    gemcat_out = RES / "gemcat_esc_over_fib_result.csv"
    if len(fixed) >= 50:
        run_gemcat(fixed_path, gemcat_out)
        gem = summarize_gemcat(gemcat_out)
        gem.to_csv(RES / "gemcat_esc_over_fib_sorted.tsv", sep="\t", index=False)
        gem.head(50).to_csv(RES / "gemcat_esc_over_fib_top50.tsv", sep="\t", index=False)
        gem.tail(50).to_csv(RES / "gemcat_esc_over_fib_bottom50.tsv", sep="\t", index=False)

        fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
        top = pd.concat([gem.head(15), gem.tail(15)])
        ax.barh(top.metabolite.astype(str), top.score)
        ax.set_xlabel("GEMCAT centrality score (ESC/fib logFC weights)")
        ax.set_title("GEMCAT Recon3D: metabolites shifted ESC vs fibroblast")
        fig.savefig(FIG / "fig_A5_gemcat_esc_fib.png", dpi=160)
        plt.close(fig)
    else:
        print("Too few mapped genes for GEMCAT; skipped.")

    # Q6–Q7 note: multi-species paired ESC/fib unavailable
    q67 = {
        "Q6_Q7_multi_spp_OXPHOS_delta_vs_mass": "NOT AVAILABLE",
        "reason": "GSE213253 is fibroblasts-only across 7 species; chicken atlas has ESC+fib but no matched ESC for other birds.",
        "chicken_only_OXPHOS_delta": None,
    }
    ox = mod_tab[mod_tab.module == "OXPHOS"]
    if len(ox):
        q67["chicken_only_OXPHOS_delta"] = float(ox.iloc[0]["delta_esc_minus_fib"])
        q67["chicken_OXPHOS_mean_gene_logFC"] = float(
            path_proxy.loc[path_proxy.module == "OXPHOS", "mean_logFC"].iloc[0]
        ) if (path_proxy.module == "OXPHOS").any() else None
    with open(RES / "A5_Q6_Q7_multispecies_unavailable.json", "w") as f:
        json.dump(q67, f, indent=2)

    print("Done A4–A5.")


if __name__ == "__main__":
    main()
