#!/usr/bin/env python3
"""Is gene length the thing that decides which gene set 'works' against body mass?

The deck shows three sets: mitochondrial and ribosomal decline with mass, cell
cycle does not. Those three differ in an obvious non-biological way — the first
two are short genes, the third is long. The exp_tab files carry `geneLength`, so
this can be tested directly.

Checks:
1. Per gene, does length predict whether the three Augustus-annotated species
   detect it at all?
2. Per set, does median gene length predict the set's correlation with mass on
   the deck's summed metric?
3. Does a length-matched null reproduce the mitochondrial and ribosomal result?
4. Where does my own symbol universe come from, and how much of the per-species
   gene count is my mapping choice rather than the genome annotation?
"""
from __future__ import annotations

import json
import re
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
    EXP, RES, RUSLAN, SPECIES, SPECIES_ORDER, annotate_table, normalize_symbol,
    nuclear_oxphos, peak_to_symbol_de)

RMD = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/sizes.Rmd")
OUT = RES / "colleagues_matrix_check"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
RNG = np.random.default_rng(20260821)
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))
AUG = ["Peacock", "Turkey", "Emu"]


def fit(y):
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan
    r = stats.linregress(LOGM[ok], y[ok])
    return float(r.rvalue), float(r.pvalue)


def rmd_list(var: str) -> list[str]:
    txt = RMD.read_text(errors="ignore")
    m = re.search(rf"^{var}\s*<-\s*c\((.*?)\)\s*$", txt, re.S | re.M)
    return sorted({normalize_symbol(t) for t in re.findall(r"[A-Za-z][A-Za-z0-9\-]*", m.group(1))} - {""})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    per_species, lengths = {}, {}
    for sp in SPECIES_ORDER:
        meta = SPECIES[sp]
        df = pd.read_csv(EXP / meta["file"], sep="\t")
        df = annotate_table(df, peak_to_symbol_de(RUSLAN / meta["de"]))
        df["symbol"] = df["symbol"].map(normalize_symbol)
        prom = df[df["class"].astype(str).str.lower().eq("promoter") & df["symbol"].ne("")]
        cols = [c for c in prom.columns if str(c).startswith("norm.C")]
        g = prom.groupby("symbol", observed=True)
        per_species[sp] = g[cols].sum().mean(axis=1)
        lengths[sp] = g["geneLength"].median()
    wide = pd.DataFrame(per_species)
    length_chicken = lengths["Chicken"]
    print(f"universe {wide.shape}, chicken gene lengths for {length_chicken.notna().sum()} symbols",
          flush=True)

    det = wide.notna()
    report: dict = {}

    # ---------- 1. does length predict detection in the Augustus species ----------
    found_aug = det[AUG].sum(axis=1)
    common = length_chicken.dropna().index.intersection(found_aug.index)
    L = np.log10(length_chicken.reindex(common).astype(float).clip(lower=1))
    F = found_aug.reindex(common).astype(float)
    rho = stats.spearmanr(L, F)
    bins = pd.qcut(L.rank(method="first"), 10, labels=False)
    by_bin = pd.DataFrame({"len": L, "found": F, "bin": bins}).groupby("bin").agg(
        median_length_bp=("len", lambda v: float(10 ** np.median(v))),
        mean_found_of_3=("found", "mean"), n=("found", "size"))
    print("\n1. detection in the three Augustus genomes by gene-length decile")
    print(by_bin.to_string())
    print(f"   spearman(log10 gene length, n Augustus species detecting) = "
          f"{rho.statistic:+.3f} p={rho.pvalue:.3g}")
    report["length_vs_detection"] = dict(spearman=float(rho.statistic), p=float(rho.pvalue),
                                         by_decile=by_bin.reset_index().to_dict("records"))

    # ---------- 2. set-level: median length vs the summed-metric correlation ----------
    sets = {"mitochondrial (deck)": rmd_list("mt_genes"),
            "cell cycle (deck)": rmd_list("cell_cycle_genes"),
            "OXPHOS (canon)": sorted({normalize_symbol(g) for g in nuclear_oxphos()} - {""}),
            "cytosolic ribosome": [g for g in wide.index if str(g).startswith(("RPL", "RPS"))],
            "mitoribosome": [g for g in wide.index if str(g).startswith(("MRPL", "MRPS"))],
            "proteasome": [g for g in wide.index if str(g).startswith("PSM")],
            "collagen": [g for g in wide.index if str(g).startswith("COL")],
            "spliceosome": [g for g in wide.index if str(g).startswith(("SNRP", "SF3"))],
            "kinases (long, control)": [g for g in wide.index if str(g).endswith("K1")]}
    rows = []
    for name, genes in sets.items():
        hit = wide.reindex(genes)
        n_present = int(hit.notna().any(axis=1).sum())
        if n_present < 8:
            continue
        r_sum, p_sum = fit(np.log10(hit.sum(axis=0).values))
        core = hit.dropna()
        r_core, p_core = fit(np.log10(core.sum(axis=0).values)) if len(core) >= 8 else (np.nan, np.nan)
        med_len = float(length_chicken.reindex(genes).median())
        rows.append(dict(gene_set=name, n_present=n_present, median_gene_length_bp=med_len,
                         r_variable_list=r_sum, p_variable_list=p_sum,
                         n_core=len(core), r_fixed_core=r_core, p_fixed_core=p_core))
        print(f"2. {name:24s} n={n_present:4d} medLen={med_len:8.0f}bp  "
              f"variable r={r_sum:+.3f} p={p_sum:.4f}  |  fixed core n={len(core):4d} "
              f"r={r_core:+.3f} p={p_core:.4f}")
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "gene_length_vs_setresult.tsv", sep="\t", index=False)
    sr = stats.spearmanr(tab["median_gene_length_bp"], tab["r_variable_list"])
    print(f"   spearman(median gene length, r on the summed metric) = "
          f"{sr.statistic:+.3f} p={sr.pvalue:.4f} over {len(tab)} sets")
    report["set_level"] = dict(rows=rows, spearman_length_vs_r=float(sr.statistic),
                               p=float(sr.pvalue))

    # ---------- 3. length-matched null for the mitochondrial set ----------
    mito = [g for g in sets["mitochondrial (deck)"] if g in wide.index]
    prof = pd.DataFrame({"len": np.log10(length_chicken.astype(float).clip(lower=1))}).dropna()
    prof["bin"] = pd.qcut(prof["len"].rank(method="first"), 20, labels=False)
    tgt = prof.reindex([g for g in mito if g in prof.index]).dropna()
    pool = prof.drop(index=[g for g in mito if g in prof.index], errors="ignore")
    obs, _ = fit(np.log10(wide.reindex(mito).sum(axis=0).values))
    nulls = []
    for _ in range(3000):
        out = []
        for b, grp in tgt.groupby("bin", observed=True):
            cand = pool.index.values[pool["bin"].values == b]
            if len(cand) == 0:
                out = None
                break
            out.append(RNG.choice(cand, size=len(grp), replace=len(cand) < len(grp)))
        if out is None:
            continue
        r, _ = fit(np.log10(wide.reindex(np.concatenate(out)).sum(axis=0).values))
        if np.isfinite(r):
            nulls.append(r)
    nulls = np.array(nulls)
    emp = float((nulls <= obs).mean())
    report["length_matched_null"] = dict(observed_r=obs, null_median=float(np.median(nulls)),
                                         empirical_p=emp, n=int(len(nulls)))
    print(f"\n3. length-matched null: observed r={obs:+.3f}, null median={np.median(nulls):+.3f}, "
          f"empirical p={emp:.4f} (n={len(nulls)})")

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    ax = axes[0]
    ax.plot(by_bin["median_length_bp"], by_bin["mean_found_of_3"], "o-", color="#c1121f")
    ax.set_xscale("log")
    ax.set(xlabel="median gene length of the decile (bp, log scale)",
           ylabel="Augustus species detecting the gene (of 3)",
           title=f"1. short genes are the ones the Augustus\ngenomes miss "
                 f"(spearman {rho.statistic:+.2f})")

    ax = axes[1]
    ax.scatter(tab["median_gene_length_bp"], tab["r_variable_list"], s=70, color="#1d3557")
    for _, r_ in tab.iterrows():
        ax.annotate(r_.gene_set, (r_.median_gene_length_bp, r_.r_variable_list), fontsize=6.5,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(0, color="#adb5bd", lw=1)
    ax.set_xscale("log")
    ax.set(xlabel="median gene length of the set (bp, log scale)",
           ylabel="r vs log10 mass, summed metric",
           title=f"2. the shorter the set, the stronger\nthe 'scaling' (spearman {sr.statistic:+.2f})")

    ax = axes[2]
    ax.hist(nulls, bins=45, color="#adb5bd", edgecolor="none")
    ax.axvline(obs, color="#c1121f", lw=2)
    ax.axvline(float(np.median(nulls)), color="#264653", lw=1.2, ls="--")
    ax.set(xlabel="r vs log10 mass", ylabel="length-matched null sets",
           title=f"3. random sets of the same gene lengths\nempirical p={emp:.3f}")

    fig.suptitle("Gene length decides which set 'scales with body mass'", fontsize=12)
    for p in (FIG / "fig_gene_length_discriminator.png", OUT / "fig_gene_length_discriminator.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "gene_length_discriminator.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
