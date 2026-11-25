#!/usr/bin/env python3
"""The deck's own metric and the deck's own gene lists, on a gene list that cannot move.

Script 32 turned up the fact that settles the argument: on the deck's summed
metric, a deliberately arbitrary long-gene set ('symbols ending in K1') gives
r=+0.866 p=0.012 against body mass, as strong as the mitochondrial set gives in
the other direction. A metric that manufactures significance of both signs from
arbitrary sets is not evidence about mitochondria.

This script closes the loop on the substrate the colleagues themselves built:
the cross-birds 1-to-1 ortholog matrix, 3055 genes x 65 libraries, zero missing
values. There the gene list is identical in all seven species by construction, so
the count term of the summed metric is frozen. If the decline is biological it
survives; if it is the count term it disappears.

Also reports what the mtDNA-encoded genes can and cannot say on this panel.
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
CROSS = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/cross-birds")
OUT = RES / "colleagues_matrix_check"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
RNG = np.random.default_rng(20260821)
N_RAND = 4000
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))

BLOCKS = [("Chicken", "fib", 0, 6), ("Chicken", "pgc", 6, 8), ("Chicken", "epi", 8, 11),
          ("Quail", "fib", 11, 17), ("Quail", "epi", 17, 20),
          ("Turkey", "fib", 20, 26), ("Turkey", "epi", 26, 29),
          ("Peacock", "fib", 29, 35), ("Peacock", "epi", 35, 39),
          ("Duck", "fib", 39, 45), ("Duck", "epi", 45, 48),
          ("Zebra", "fib", 48, 52), ("Zebra", "epi", 52, 55),
          ("Emu", "fib", 55, 61), ("Emu", "epi", 61, 65)]


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
    if not m:
        return []
    return sorted({normalize_symbol(t) for t in re.findall(r"[A-Za-z][A-Za-z0-9\-]*", m.group(1))} - {""})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {}

    # ---------- our CAGE assembly, for the decomposition and the mtDNA count ----------
    per, mtdna = {}, {}
    for sp in SPECIES_ORDER:
        meta = SPECIES[sp]
        df = pd.read_csv(EXP / meta["file"], sep="\t")
        df = annotate_table(df, peak_to_symbol_de(RUSLAN / meta["de"]))
        df["symbol"] = df["symbol"].map(normalize_symbol)
        prom = df[df["class"].astype(str).str.lower().eq("promoter") & df["symbol"].ne("")]
        cols = [c for c in prom.columns if str(c).startswith("norm.C")]
        per[sp] = prom.groupby("symbol", observed=True)[cols].sum().mean(axis=1)
        mtdna[sp] = sorted({s for s in per[sp].index if str(s).startswith("MT-")})
    wide = pd.DataFrame(per)
    report["mtdna_genes_per_species"] = {sp: v for sp, v in mtdna.items()}
    print("mtDNA-encoded genes present in the CAGE promoter tables:")
    for sp in SPECIES_ORDER:
        print(f"  {sp:8s} n={len(mtdna[sp])}  {mtdna[sp][:6]}")

    # ---------- the see-saw: sum metric on arbitrary sets, both directions ----------
    mito = rmd_list("mt_genes")
    cc = rmd_list("cell_cycle_genes")
    canon = sorted({normalize_symbol(g) for g in nuclear_oxphos()} - {""})
    arb = {"mitochondrial (deck)": mito, "cell cycle (deck)": cc, "OXPHOS (canon)": canon,
           "cytosolic ribosome": [g for g in wide.index if str(g).startswith(("RPL", "RPS"))],
           "arbitrary: symbols ending in K1": [g for g in wide.index if str(g).endswith("K1")],
           "arbitrary: symbols starting with ZN": [g for g in wide.index if str(g).startswith("ZN")],
           "arbitrary: symbols starting with SLC": [g for g in wide.index if str(g).startswith("SLC")]}
    rows = []
    for name, genes in arb.items():
        hit = wide.reindex([g for g in genes if g in wide.index])
        if len(hit) < 10:
            continue
        n_det = hit.notna().sum(axis=0).values.astype(float)
        mean_det = hit.mean(axis=0).values
        r_sum, p_sum = fit(np.log10(hit.sum(axis=0).values))
        r_n, p_n = fit(np.log10(n_det))
        r_m, p_m = fit(np.log10(mean_det))
        rows.append(dict(gene_set=name, n_list=len(hit), r_sum=r_sum, p_sum=p_sum,
                         r_count=r_n, p_count=p_n, r_mean=r_m, p_mean=p_m,
                         n_detected=n_det.astype(int).tolist()))
        print(f"\n{name:34s} n={len(hit):4d}")
        print(f"   log10 sum   r={r_sum:+.3f} p={p_sum:.4f}"
              f"   = count r={r_n:+.3f} p={p_n:.4f}  +  mean r={r_m:+.3f} p={p_m:.4f}")
    seesaw = pd.DataFrame(rows)
    seesaw.to_csv(OUT / "seesaw_arbitrary_sets.tsv", sep="\t", index=False)
    report["seesaw"] = rows

    # how often does the summed metric hand out p<0.05 to a random set of 123 genes?
    pool = wide.index.values
    hits, rs = 0, []
    for _ in range(N_RAND):
        draw = RNG.choice(pool, size=len(mito), replace=False)
        r, p = fit(np.log10(wide.reindex(draw).sum(axis=0).values))
        if np.isfinite(r):
            rs.append(r)
            hits += int(p < 0.05)
    rs = np.array(rs)
    report["random_set_false_positive_rate"] = dict(
        n=len(rs), frac_p_below_005=hits / max(len(rs), 1), median_r=float(np.median(rs)),
        frac_r_below_minus_08=float((rs <= -0.8).mean()), set_size=len(mito))
    print(f"\nrandom sets of {len(mito)} genes on the summed metric: "
          f"{100*hits/len(rs):.1f}% reach p<0.05, median r={np.median(rs):+.3f}")

    # ---------- the frozen list: colleagues' own ortholog matrix ----------
    genes = pd.read_csv(CROSS / "counts.csv")["0"].map(normalize_symbol).values
    mat = pd.read_csv(CROSS / "file.csv", index_col=0)
    mat.index = genes
    mat = mat[~mat.index.duplicated()]
    per_t: dict[str, dict[str, pd.Series]] = {}
    for sp, tis, a, b in BLOCKS:
        per_t.setdefault(tis, {})[sp] = mat.iloc[:, a:b].mean(axis=1)
    cfib = pd.DataFrame({sp: per_t["fib"][sp] for sp in SPECIES_ORDER})
    cepi = pd.DataFrame({sp: per_t["epi"][sp] for sp in SPECIES_ORDER})
    print(f"\northolog matrix {cfib.shape}, missing {int(cfib.isna().sum().sum())}")

    frows = []
    for name, gl in (("mitochondrial (deck)", mito), ("cell cycle (deck)", cc),
                     ("OXPHOS (canon)", canon),
                     ("cytosolic ribosome", [g for g in cfib.index if str(g).startswith(("RPL", "RPS"))])):
        present = [g for g in gl if g in cfib.index]
        if len(present) < 8:
            print(f"  {name:24s} only {len(present)} orthologs, skipped")
            continue
        for tis, m in (("fib", cfib), ("epi", cepi)):
            sub = m.loc[present]
            r_sum, p_sum = fit(np.log10(sub.sum(axis=0).values))
            r_mean, p_mean = fit(sub.mean(axis=0).values)
            frows.append(dict(gene_set=name, tissue=tis, n_orthologs=len(present),
                              r_sum=r_sum, p_sum=p_sum, r_mean=r_mean, p_mean=p_mean))
            print(f"  {name:24s} {tis} n={len(present):3d}  sum r={r_sum:+.3f} p={p_sum:.4f}"
                  f"   mean r={r_mean:+.3f} p={p_mean:.4f}")
    frozen = pd.DataFrame(frows)
    frozen.to_csv(OUT / "frozen_list_same_metric.tsv", sep="\t", index=False)
    report["frozen_list"] = frows

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)

    ax = axes[0]
    s2 = seesaw.sort_values("r_sum")
    cols = ["#c1121f" if v < 0 else "#2a9d8f" for v in s2.r_sum]
    ax.barh(range(len(s2)), s2.r_sum.values, color=cols)
    for i, (rv, pv) in enumerate(zip(s2.r_sum, s2.p_sum)):
        ax.annotate(f"p={pv:.3f}", (rv, i), fontsize=7, va="center",
                    ha="right" if rv < 0 else "left",
                    xytext=(-4 if rv < 0 else 4, 0), textcoords="offset points")
    ax.set_yticks(range(len(s2)))
    ax.set_yticklabels(s2.gene_set, fontsize=7)
    ax.axvline(0, color="#333", lw=1)
    ax.set(xlim=(-1.15, 1.15), xlabel="r vs log10 mass, the deck's summed metric",
           title="the metric gives strong correlations\nof both signs to arbitrary sets")

    ax = axes[1]
    ax.hist(rs, bins=50, color="#adb5bd", edgecolor="none")
    obs_r = float(seesaw.loc[seesaw.gene_set == "mitochondrial (deck)", "r_sum"].iloc[0])
    ax.axvline(obs_r, color="#c1121f", lw=2)
    ax.axvline(0, color="#333", lw=1)
    ax.set(xlabel=f"r vs log10 mass, random sets of {len(mito)} genes",
           ylabel="random sets",
           title=f"{100*hits/len(rs):.0f}% of random sets reach p<0.05\n"
                 f"null median r={np.median(rs):+.2f}, not 0")

    ax = axes[2]
    f2 = frozen[frozen.tissue == "fib"]
    x = np.arange(len(f2))
    w = 0.38
    ax.bar(x - w / 2, f2.r_sum, w, color="#c1121f", label="summed metric")
    ax.bar(x + w / 2, f2.r_mean, w, color="#457b9d", label="mean per gene")
    ax.axhline(0, color="#333", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{g}\n(n={n})" for g, n in zip(f2.gene_set, f2.n_orthologs)], fontsize=7)
    ax.set(ylim=(-1.05, 1.05), ylabel="r vs log10 mass (fibroblasts)",
           title="on the colleagues' 1-to-1 ortholog matrix\nthe gene list cannot move")
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("Freeze the gene list and the scaling goes away", fontsize=12)
    for p in (FIG / "fig_same_metric_fixed_list.png", OUT / "fig_same_metric_fixed_list.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "same_metric_fixed_list.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
