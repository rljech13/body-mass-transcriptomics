#!/usr/bin/env python3
"""The mechanism behind the CAGE fibroblast OXPHOS decline, and why the null missed it.

Findings that this script establishes:

1. Breadth gradient. The fit strengthens monotonically as genes with more dropout
   are admitted: breadth>=7 p=0.58, >=6 p=0.019, >=5 p=0.0035. Significance is
   manufactured by dropout, not damaged by it.
2. Inverted dropout. Genome wide, the genes a species fails to detect are its low
   expressed ones. Among OXPHOS the missing genes are the high expressed ones, so
   removing them drags the mean of the surviving set down.
3. Gene property, not species biology. The same nine subunits are absent from all
   three Augustus annotated genomes independently.
4. Why the script 25 null missed it. For a random gene, breadth and abundance are
   positively coupled, so a low-breadth draw is also low-expressed and its loss
   barely moves the mean. OXPHOS breaks that coupling. Matching on breadth and a
   chicken abundance DECILE cannot separate the two, because every OXPHOS gene
   that matters sits inside the top decile. Rematching on breadth by fine
   abundance rank inside the top decile removes the result.
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
from oxphos_thread_lib import RES, SPECIES, SPECIES_ORDER, normalize_symbol, nuclear_oxphos  # noqa: E402

OUT = RES / "explain_cage_decline"
CACHE = RES / "reverify_confounds" / "gene_matrices.pkl"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
RNG = np.random.default_rng(20260820)
N_NULL = 3000
AUGUSTUS = ("Peacock", "Turkey", "Emu")
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))


def fit(y: np.ndarray) -> tuple[float, float]:
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3 or np.ptp(y[ok]) == 0:
        return np.nan, np.nan
    res = stats.linregress(LOGM[ok], y[ok])
    return float(res.rvalue), float(res.pvalue)


def y1(wide: pd.DataFrame, genes) -> np.ndarray:
    return wide.reindex(genes).mean(axis=0).values


def draw_matched(prof: pd.DataFrame, target, keys: list[str]):
    tgt = prof.reindex(target).dropna(subset=keys)
    pool = prof.drop(index=[g for g in target if g in prof.index], errors="ignore")
    out = []
    for _, grp in tgt.groupby(keys, observed=True):
        mask = np.ones(len(pool), dtype=bool)
        for k in keys:
            mask &= (pool[k].values == grp.iloc[0][k])
        cand = pool.index.values[mask]
        if len(cand) == 0:
            return None
        out.append(RNG.choice(cand, size=len(grp), replace=len(cand) < len(grp)))
    return np.concatenate(out)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with CACHE.open("rb") as fh:
        fib = pickle.load(fh)["CAGE_fib"]
    wide = pd.DataFrame({sp: fib[sp][0] for sp in SPECIES_ORDER})
    ox_all = {normalize_symbol(g) for g in nuclear_oxphos()} - {""}
    ox = [g for g in wide.index if g in ox_all]
    det = wide.notna()
    breadth = det.sum(axis=1)
    chick = wide["Chicken"]
    report: dict = {}

    # ---------- 1. breadth gradient ----------
    grad = []
    for k in range(7, 2, -1):
        keep = [g for g in ox if breadth[g] >= k]
        r, p = fit(y1(wide, keep))
        grad.append(dict(min_breadth=k, n_genes=len(keep), r=r, p=p))
        print(f"1. breadth>={k}: n={len(keep):3d} r={r:+.3f} p={p:.4f}")
    report["breadth_gradient"] = grad

    # ---------- 2. inverted dropout ----------
    inv = []
    for sp in SPECIES_ORDER:
        kept = chick.reindex([g for g in ox if det.loc[g, sp]]).dropna()
        lost = chick.reindex([g for g in ox if not det.loc[g, sp]]).dropna()
        gk = chick[det[sp] & chick.notna()]
        gl = chick[~det[sp] & chick.notna()]
        inv.append(dict(species=sp, n_kept=len(kept), n_lost=len(lost),
                        oxphos_kept_median=float(kept.median()),
                        oxphos_lost_median=float(lost.median()) if len(lost) else np.nan,
                        genome_kept_median=float(gk.median()),
                        genome_lost_median=float(gl.median()) if len(gl) else np.nan))
    inv_df = pd.DataFrame(inv)
    inv_df.to_csv(OUT / "dropout_inversion.tsv", sep="\t", index=False)
    report["dropout_inversion"] = inv
    print("\n2. dropout inversion (chicken level of kept vs lost genes)")
    print(inv_df.to_string(index=False))

    # coupling of breadth and abundance: genome wide vs inside OXPHOS
    ok = chick.notna()
    rho_all = stats.spearmanr(breadth[ok], chick[ok])
    ox_ok = [g for g in ox if pd.notna(chick.get(g))]
    rho_ox = stats.spearmanr(breadth[ox_ok], chick[ox_ok])
    report["breadth_abundance_coupling"] = dict(
        genome_rho=float(rho_all.statistic), genome_p=float(rho_all.pvalue),
        oxphos_rho=float(rho_ox.statistic), oxphos_p=float(rho_ox.pvalue), n_oxphos=len(ox_ok))
    print(f"\n   spearman(breadth, chicken level): genome {rho_all.statistic:+.3f} "
          f"(p={rho_all.pvalue:.3g}) vs OXPHOS {rho_ox.statistic:+.3f} (p={rho_ox.pvalue:.3g})")

    # ---------- 3. same genes lost in independent genomes ----------
    lost_sets = {sp: set(det.index[~det[sp]]) & set(ox) for sp in AUGUSTUS}
    shared_lost = sorted(set.intersection(*lost_sets.values()))
    exp_by_chance = np.prod([len(lost_sets[sp]) / len(ox) for sp in AUGUSTUS]) * len(ox)
    report["shared_dropout"] = dict(
        genes=shared_lost, n_observed=len(shared_lost),
        n_expected_if_independent=float(exp_by_chance),
        per_species_lost={sp: len(v) for sp, v in lost_sets.items()}, n_oxphos=len(ox))
    print(f"\n3. lost in all three Augustus genomes: {len(shared_lost)} vs "
          f"{exp_by_chance:.1f} expected if independent -> {shared_lost}")

    # ---------- 4. the null, at decile vs fine abundance resolution ----------
    prof = pd.DataFrame({"breadth": breadth, "abundance": chick})
    seen = prof["abundance"].notna()
    prof["abund_decile"] = -1
    prof.loc[seen, "abund_decile"] = pd.qcut(
        prof.loc[seen, "abundance"].rank(method="first"), 10, labels=False).astype(int)
    prof["abund_fine"] = -1
    prof.loc[seen, "abund_fine"] = pd.qcut(
        prof.loc[seen, "abundance"].rank(method="first"), 100, labels=False).astype(int)

    obs_r, obs_p = fit(y1(wide, ox))
    ox_arr = np.array(ox)
    res_null = {}
    for label, keys in (("breadth + abundance decile", ["breadth", "abund_decile"]),
                        ("breadth + fine abundance percentile", ["breadth", "abund_fine"])):
        vals = []
        for _ in range(N_NULL):
            d = draw_matched(prof, ox_arr, keys)
            if d is None:
                continue
            r, _ = fit(y1(wide, d))
            if np.isfinite(r):
                vals.append(r)
        vals = np.array(vals)
        emp = float((vals <= obs_r).mean())
        res_null[label] = dict(empirical_p=emp, n_null=int(len(vals)),
                               null_median_r=float(np.median(vals)),
                               null_q05_r=float(np.quantile(vals, 0.05)))
        print(f"\n4. null on {label}: n={len(vals)} median r={np.median(vals):+.3f} "
              f"observed r={obs_r:+.3f} empirical p={emp:.4f}")
    report["nulls"] = dict(observed_r=obs_r, naive_p=obs_p, **res_null)

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    g = pd.DataFrame(grad)
    ax = axes[0]
    ax.plot(g.min_breadth, -np.log10(g.p), "o-", color="#c1121f")
    ax.axhline(-np.log10(0.05), color="#adb5bd", ls="--", lw=1)
    ax.invert_xaxis()
    for _, row in g.iterrows():
        ax.annotate(f"n={row.n_genes:.0f}", (row.min_breadth, -np.log10(row.p)),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set(xlabel="genes admitted: detected in at least k of 7 species",
           ylabel="-log10 p", title="1. the result is created by dropout\n"
                                    "strict gene list = no signal")

    ax = axes[1]
    x = np.arange(len(inv_df))
    ax.plot(x, inv_df.oxphos_kept_median, "o-", color="#c1121f", label="OXPHOS kept")
    ax.plot(x, inv_df.oxphos_lost_median, "s--", color="#e07a5f", label="OXPHOS lost")
    ax.plot(x, inv_df.genome_kept_median, "o-", color="#264653", label="all genes kept")
    ax.plot(x, inv_df.genome_lost_median, "s--", color="#2a9d8f", label="all genes lost")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_ORDER, rotation=45, ha="right", fontsize=8)
    ax.set(ylabel="chicken CAGE level (log)",
           title="2. genome wide the lost genes are weak,\nfor OXPHOS they are the strong ones")
    ax.legend(fontsize=7, frameon=False)

    ax = axes[2]
    labels = list(res_null)
    ax.bar(range(2), [res_null[k]["empirical_p"] for k in labels],
           color=["#e07a5f", "#2a9d8f"])
    ax.axhline(0.05, color="#c1121f", ls="--", lw=1)
    ax.set_xticks(range(2))
    ax.set_xticklabels(["decile\nmatching", "fine percentile\nmatching"], fontsize=8)
    ax.set(ylabel="empirical p", title="4. the null was too coarse")
    for i, k in enumerate(labels):
        ax.annotate(f"{res_null[k]['empirical_p']:.3f}", (i, res_null[k]["empirical_p"]),
                    ha="center", va="bottom", fontsize=9)

    fig.suptitle("The CAGE fibroblast OXPHOS decline is annotation dropout", fontsize=12)
    for p in (FIG / "fig_dropout_mechanism.png", OUT / "fig_dropout_mechanism.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "dropout_mechanism.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
