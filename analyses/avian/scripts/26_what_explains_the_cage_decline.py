#!/usr/bin/env python3
"""What, if anything, explains the surviving CAGE fibroblast OXPHOS decline?

Script 25 left the CAGE fibroblast decline standing against nulls matched on
detection breadth and chicken abundance. This script tests the explanations that
were never put to the data:

A. Annotation platform. The symbol universe collapses from ~10k to ~6.4-7.6k in
   exactly the three heaviest species (Peacock, Turkey, Emu use Augustus models),
   so "large" and "poorly annotated" are the same binary variable here.
A2. The paired fibroblast-minus-epiblast contrast. Both tissues come from the same
   peak tables, the same symbol universe and the same normalisation, so any
   explanation living in the annotation or the denominator cancels in the
   difference.
B. Module coherence. The script 25 null matched breadth and abundance but not
   co-regulation. A tightly co-expressed module produces a cleaner species trend
   than a random mixture, so it can beat that null on coherence alone.
C. Specificity. Whether other conserved housekeeping modules fall the same way.
D. What rises. If the OXPHOS share falls, something has to take its place.
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
    RES, SPECIES, SPECIES_ORDER, load_geneset, normalize_symbol, nuclear_oxphos)

OUT = RES / "explain_cage_decline"
CACHE = RES / "reverify_confounds" / "gene_matrices.pkl"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
RNG = np.random.default_rng(20260820)
N_NULL = 2000
N_NULL_FAMILY = 600

# Peacock/Turkey/Emu carry Augustus gene models; symbols exist only where the DE
# table could map them. The other four have native symbol annotation.
AUGUSTUS = {"Peacock", "Turkey", "Emu"}

LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))
IS_AUG = np.array([sp in AUGUSTUS for sp in SPECIES_ORDER])

FAMILIES = {
    "cytosolic ribosome (RPL/RPS)": ("prefix", ("RPL", "RPS")),
    "mitoribosome (MRPL/MRPS)": ("prefix", ("MRPL", "MRPS")),
    "proteasome (PSM)": ("prefix", ("PSM",)),
    "translation initiation (EIF)": ("prefix", ("EIF",)),
    "spliceosome (SNRP/SF3)": ("prefix", ("SNRP", "SF3")),
    "heat shock (HSP)": ("prefix", ("HSP",)),
    "collagen (COL)": ("prefix", ("COL",)),
    "complex I (NDUF)": ("prefix", ("NDUF",)),
    "complex IV (COX)": ("prefix", ("COX",)),
    "ATP synthase (ATP5)": ("prefix", ("ATP5",)),
}


def fit(y: np.ndarray, x: np.ndarray = LOGM) -> tuple[float, float, float]:
    ok = np.isfinite(y) & np.isfinite(x)
    if ok.sum() < 3 or np.ptp(y[ok]) == 0 or np.ptp(x[ok]) == 0:
        return np.nan, np.nan, np.nan
    res = stats.linregress(x[ok], y[ok])
    return float(res.rvalue), float(res.pvalue), float(res.slope)


def y1(layer: dict, genes: np.ndarray) -> np.ndarray:
    """Mean level over the genes of the set that are detected in that species."""
    out = []
    for sp in SPECIES_ORDER:
        hit = layer[sp][0].reindex(genes).dropna()
        out.append(float(hit.mean()) if len(hit) else np.nan)
    return np.array(out)


def build_gene_profile(layer: dict) -> pd.DataFrame:
    breadth: dict[str, int] = {}
    for sp in SPECIES_ORDER:
        for g in layer[sp][0].index:
            breadth[g] = breadth.get(g, 0) + 1
    prof = pd.DataFrame({"symbol": list(breadth), "breadth": list(breadth.values())})
    prof["abundance"] = prof["symbol"].map(layer["Chicken"][0]).astype(float)
    seen = prof["abundance"].notna()
    prof["abund_decile"] = -1
    prof.loc[seen, "abund_decile"] = pd.qcut(
        prof.loc[seen, "abundance"].rank(method="first"), 10, labels=False).astype(int)
    return prof.set_index("symbol")


def matched_draw(prof: pd.DataFrame, target: np.ndarray, keys: list[str]) -> np.ndarray | None:
    tgt = prof.reindex(target).dropna(subset=keys)
    if tgt.empty:
        return None
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


def wide_matrix(layer: dict) -> pd.DataFrame:
    return pd.DataFrame({sp: layer[sp][0] for sp in SPECIES_ORDER})


def coherence(zmat: pd.DataFrame, genes: np.ndarray) -> float:
    """Mean pairwise correlation of across-species profiles, genes seen everywhere."""
    sub = zmat.reindex(genes).dropna()
    if len(sub) < 5:
        return np.nan
    c = np.corrcoef(sub.values)
    iu = np.triu_indices_from(c, k=1)
    return float(np.nanmean(c[iu]))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with CACHE.open("rb") as fh:
        layers = pickle.load(fh)

    fib, epi = layers["CAGE_fib"], layers["CAGE_epi"]
    oxphos = np.array(sorted({normalize_symbol(g) for g in nuclear_oxphos()} - {""}))
    module = np.array(sorted({normalize_symbol(g) for g in pd.read_csv(
        RES / "definitions" / "oxphos_module.tsv", sep="\t")["gene"]} - {""}))

    prof = build_gene_profile(fib)
    present = np.array([g for g in oxphos if g in prof.index])
    mod_present = np.array([g for g in module if g in prof.index])

    wide_fib, wide_epi = wide_matrix(fib), wide_matrix(epi)
    lg = np.log10(wide_fib + 1.0)
    zmat = lg.dropna().pipe(lambda d: d.sub(d.mean(axis=1), axis=0).div(d.std(axis=1).replace(0, np.nan), axis=0)).dropna()
    print(f"universe: {len(prof)} symbols, {len(zmat)} detected in all 7 species; "
          f"OXPHOS present {len(present)}, module {len(mod_present)}", flush=True)

    report: dict = {}

    # ---------- A. annotation platform ----------
    n_sym = np.array([len(fib[sp][0]) for sp in SPECIES_ORDER], dtype=float)
    r_ann, p_ann, _ = fit(n_sym)
    obs_fib = y1(fib, present)
    obs_epi = y1(epi, present)

    nat = ~IS_AUG
    r_all, p_all, _ = fit(obs_fib)
    r_nat, p_nat, _ = fit(obs_fib[nat], LOGM[nat])
    r_aug, p_aug, _ = fit(obs_fib[IS_AUG], LOGM[IS_AUG])

    # within-group centring removes the platform step entirely
    yc, xc = obs_fib.copy(), LOGM.copy()
    for grp in (nat, IS_AUG):
        yc[grp] -= yc[grp].mean()
        xc[grp] -= xc[grp].mean()
    r_wg, p_wg, _ = fit(yc, xc)
    between = float(np.mean(obs_fib[nat]) - np.mean(obs_fib[IS_AUG]))

    report["A_annotation"] = dict(
        n_symbols=dict(zip(SPECIES_ORDER, n_sym.tolist())),
        n_symbols_vs_mass=dict(r=r_ann, p=p_ann),
        platform_split_is_perfect=bool(LOGM[nat].max() < LOGM[IS_AUG].min()),
        full_r=r_all, full_p=p_all,
        native4_r=r_nat, native4_p=p_nat, native4_mass_range_fold=float(
            10 ** (LOGM[nat].max() - LOGM[nat].min())),
        augustus3_r=r_aug, augustus3_p=p_aug,
        within_group_centred_r=r_wg, within_group_centred_p=p_wg,
        between_group_gap=between)
    print(f"\nA. symbol universe vs mass: r={r_ann:+.3f} p={p_ann:.4f}")
    print(f"   full n=7      r={r_all:+.3f} p={p_all:.4g}")
    print(f"   native4 only  r={r_nat:+.3f} p={p_nat:.4g}  ({report['A_annotation']['native4_mass_range_fold']:.0f}x mass)")
    print(f"   augustus3     r={r_aug:+.3f} p={p_aug:.4g}")
    print(f"   within-group  r={r_wg:+.3f} p={p_wg:.4g}   between-group gap={between:+.2f}")

    # ---------- A2. paired fib - epi contrast ----------
    rows = []
    for name, genes in (("set", present), ("module", mod_present)):
        f, e = y1(fib, genes), y1(epi, genes)
        for label, vec in (("fib", f), ("epi", e), ("fib-epi", f - e),
                           ("log2 fib/epi", np.log2(f / e))):
            r, p, s = fit(vec)
            rows.append(dict(definition=name, contrast=label, r=r, p=p, slope=s,
                             values=np.round(vec, 3).tolist()))
            print(f"A2 {name:7s} {label:13s} r={r:+.3f} p={p:.4g}")
    paired = pd.DataFrame(rows)
    paired.to_csv(OUT / "paired_fib_epi_contrast.tsv", sep="\t", index=False)
    report["A2_paired"] = rows

    # ---------- B. coherence-matched null ----------
    obs_r = r_all
    obs_coh = coherence(zmat, present)
    null_r, null_coh = [], []
    for i in range(N_NULL):
        draw = matched_draw(prof, present, ["breadth", "abund_decile"])
        if draw is None:
            continue
        r, _, _ = fit(y1(fib, draw))
        c = coherence(zmat, draw)
        if np.isfinite(r) and np.isfinite(c):
            null_r.append(r)
            null_coh.append(c)
        if (i + 1) % 500 == 0:
            print(f"   null {i+1}/{N_NULL}", flush=True)
    null_r, null_coh = np.array(null_r), np.array(null_coh)
    p_flat = float((null_r <= obs_r).mean())
    # conditional on comparable coherence
    win = null_coh >= np.quantile(null_coh, 0.90)
    p_coh = float((null_r[win] <= obs_r).mean()) if win.sum() >= 50 else np.nan
    rho_coh_r = float(stats.spearmanr(null_coh, np.abs(null_r)).statistic)
    report["B_coherence"] = dict(
        observed_r=obs_r, observed_coherence=obs_coh,
        null_coherence_median=float(np.median(null_coh)),
        null_coherence_q90=float(np.quantile(null_coh, 0.90)),
        observed_coherence_percentile=float((null_coh <= obs_coh).mean()),
        spearman_coherence_vs_abs_r=rho_coh_r,
        p_breadth_abundance=p_flat, p_top_decile_coherence=p_coh,
        n_null=int(len(null_r)), n_in_window=int(win.sum()))
    print(f"\nB. observed coherence {obs_coh:+.3f} (null median {np.median(null_coh):+.3f}, "
          f"percentile {(null_coh <= obs_coh).mean():.3f})")
    print(f"   spearman(coherence, |r|) in nulls = {rho_coh_r:+.3f}")
    print(f"   p vs breadth+abundance null = {p_flat:.4f}; "
          f"p vs top-coherence-decile null = {p_coh:.4f}")

    # ---------- C. specificity across conserved modules ----------
    sets: dict[str, np.ndarray] = {"OXPHOS set": present, "OXPHOS module": mod_present}
    for nm in ("TCA", "GLYCOLYSIS", "MITO_BIOGENESIS"):
        sets[nm] = np.array([g for g in sorted({normalize_symbol(x) for x in load_geneset(nm)} - {""})
                             if g in prof.index])
    uni = prof.index.astype(str)
    for label, (_, prefixes) in FAMILIES.items():
        sets[label] = uni[[s.startswith(prefixes) for s in uni]].values

    frows = []
    for label, genes in sets.items():
        if len(genes) < 8:
            continue
        r, p, _ = fit(y1(fib, genes))
        re_, pe_, _ = fit(y1(epi, genes))
        rd, pd_, _ = fit(y1(fib, genes) - y1(epi, genes))
        coh = coherence(zmat, genes)
        nulls = []
        for _ in range(N_NULL_FAMILY):
            draw = matched_draw(prof, genes, ["breadth", "abund_decile"])
            if draw is None:
                continue
            rr, _, _ = fit(y1(fib, draw))
            if np.isfinite(rr):
                nulls.append(rr)
        emp = float((np.array(nulls) <= r).mean()) if len(nulls) >= 100 else np.nan
        frows.append(dict(gene_set=label, n_genes=int(len(genes)), coherence=coh,
                          fib_r=r, fib_p=p, matched_p=emp, epi_r=re_, epi_p=pe_,
                          fib_minus_epi_r=rd, fib_minus_epi_p=pd_))
        print(f"C {label:30s} n={len(genes):4d} coh={coh:+.2f} fib r={r:+.3f} "
              f"p={p:.4f} matched={emp:.3f} epi r={re_:+.3f} diff p={pd_:.4f}", flush=True)
    fam = pd.DataFrame(frows).sort_values("fib_r")
    fam.to_csv(OUT / "gene_family_sweep.tsv", sep="\t", index=False)
    report["C_families"] = frows

    # ---------- D. what rises ----------
    keep = wide_fib.notna().sum(axis=1) >= 6
    sub = np.log10(wide_fib[keep] + 1.0)
    slopes, pvals = [], []
    for vals in sub.values:
        r, p, s = fit(vals)
        slopes.append(s)
        pvals.append(p)
    genes_tab = pd.DataFrame({"symbol": sub.index, "slope_log10_per_log10mass": slopes,
                              "p": pvals}).sort_values("slope_log10_per_log10mass")
    genes_tab["is_oxphos"] = genes_tab["symbol"].isin(set(present))
    genes_tab.to_csv(OUT / "per_gene_slopes_fib.tsv", sep="\t", index=False)
    rank = stats.mannwhitneyu(
        genes_tab.loc[genes_tab.is_oxphos, "slope_log10_per_log10mass"],
        genes_tab.loc[~genes_tab.is_oxphos, "slope_log10_per_log10mass"],
        alternative="less")
    report["D_ranks"] = dict(
        n_genes_tested=int(len(genes_tab)),
        oxphos_n=int(genes_tab.is_oxphos.sum()),
        oxphos_median_slope=float(genes_tab.loc[genes_tab.is_oxphos, "slope_log10_per_log10mass"].median()),
        other_median_slope=float(genes_tab.loc[~genes_tab.is_oxphos, "slope_log10_per_log10mass"].median()),
        mannwhitney_p=float(rank.pvalue),
        top_risers=genes_tab.tail(25)["symbol"].tolist()[::-1],
        top_fallers=genes_tab.head(25)["symbol"].tolist())
    print(f"\nD. per-gene slopes on {len(genes_tab)} genes; OXPHOS median "
          f"{report['D_ranks']['oxphos_median_slope']:+.3f} vs rest "
          f"{report['D_ranks']['other_median_slope']:+.3f}, MW p={rank.pvalue:.3g}")
    print("   top risers:", ", ".join(report["D_ranks"]["top_risers"][:12]))
    print("   top fallers:", ", ".join(report["D_ranks"]["top_fallers"][:12]))

    # ---------- figure ----------
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)

    ax = axes[0, 0]
    for grp, col, lab in ((nat, "#1d3557", "native symbols"), (IS_AUG, "#e07a5f", "Augustus models")):
        ax.scatter(LOGM[grp], obs_fib[grp], s=70, color=col, label=lab, zorder=3)
    for sp, x, y in zip(SPECIES_ORDER, LOGM, obs_fib):
        ax.annotate(sp, (x, y), fontsize=7, xytext=(3, 4), textcoords="offset points")
    xs = np.linspace(LOGM.min(), LOGM.max(), 10)
    sl = stats.linregress(LOGM, obs_fib)
    ax.plot(xs, sl.intercept + sl.slope * xs, color="#adb5bd", lw=1.2, zorder=1)
    ax.set(xlabel="log10 body mass (g)", ylabel="OXPHOS mean CAGE level (fib)",
           title=f"A. platform and mass are the same split\nall 7 p={p_all:.4g}, "
                 f"within-platform p={p_wg:.3g}")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[0, 1]
    d = y1(fib, present) - y1(epi, present)
    ax.axhline(0, color="#adb5bd", lw=1)
    ax.scatter(LOGM, d, s=70, color="#2a9d8f", zorder=3)
    rd, pdv, _ = fit(d)
    sl = stats.linregress(LOGM, d)
    ax.plot(xs, sl.intercept + sl.slope * xs, color="#adb5bd", lw=1.2)
    for sp, x, y in zip(SPECIES_ORDER, LOGM, d):
        ax.annotate(sp, (x, y), fontsize=7, xytext=(3, 4), textcoords="offset points")
    ax.set(xlabel="log10 body mass (g)", ylabel="fibroblast - epiblast",
           title=f"A2. paired contrast cancels annotation\nr={rd:+.3f} p={pdv:.4g}")

    ax = axes[1, 0]
    ax.scatter(null_coh, null_r, s=6, alpha=0.25, color="#adb5bd")
    ax.axhline(obs_r, color="#c1121f", lw=1.5)
    ax.axvline(obs_coh, color="#c1121f", lw=1.5, ls="--")
    ax.set(xlabel="set coherence (mean pairwise r across species)", ylabel="r vs log10 mass",
           title=f"B. coherence-matched null\np(flat)={p_flat:.3f}, p(top-coherence)={p_coh:.3f}")

    ax = axes[1, 1]
    f2 = fam.dropna(subset=["fib_r"])
    cols = ["#c1121f" if "OXPHOS" in s or s.startswith(("complex", "ATP")) else "#457b9d"
            for s in f2.gene_set]
    ax.barh(range(len(f2)), f2.fib_r.values, color=cols)
    ax.set_yticks(range(len(f2)))
    ax.set_yticklabels(f2.gene_set, fontsize=7)
    ax.axvline(0, color="#333", lw=1)
    ax.set(xlabel="r of set mean vs log10 mass (fibroblasts)",
           title="C. is the fall specific to OXPHOS?")

    fig.suptitle("What explains the CAGE fibroblast OXPHOS decline", fontsize=12)
    for p in (FIG / "fig_explain_cage_decline.png", OUT / "fig_explain_cage_decline.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "summary.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
