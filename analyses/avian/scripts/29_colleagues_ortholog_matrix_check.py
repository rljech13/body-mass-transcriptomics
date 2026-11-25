#!/usr/bin/env python3
"""Re-run the OXPHOS-vs-mass question on the colleagues' 1-to-1 ortholog matrix.

Scripts 26-28 concluded that mass and annotation platform cannot be separated on
this panel, because the symbol universe collapses in exactly the three heaviest
species. That conclusion overlooked a substrate that does not have the problem:
Ruslan's cross-birds matrix, 3055 genes x 65 libraries, one row per 1-to-1
ortholog, present in all seven species with zero missing values.

On that matrix there is no dropout, no per-species gene universe, no annotation
asymmetry and no compositional denominator. If the fibroblast OXPHOS decline is
real it must appear here; if it is an artefact of the variable gene list it must
vanish. This script runs the same battery: level, within-species rank, matched
null, gene family sweep, and the fibroblast-versus-epiblast contrast.
"""
from __future__ import annotations

import json
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

CROSS = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/cross-birds")
OUT = RES / "colleagues_matrix_check"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
RNG = np.random.default_rng(20260820)
N_NULL = 4000

COLUMNS = (
    ["Fibroblast.C07", "Fibroblast.C08", "Fibroblast.C09", "Fibroblast.C10",
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
     "Epiblast.C08", "Epiblast.C09", "Epiblast.C10", "Epiblast.C11"])

# column blocks in file order
BLOCKS = [("Chicken", "fib", 0, 6), ("Chicken", "pgc", 6, 8), ("Chicken", "epi", 8, 11),
          ("Quail", "fib", 11, 17), ("Quail", "epi", 17, 20),
          ("Turkey", "fib", 20, 26), ("Turkey", "epi", 26, 29),
          ("Peacock", "fib", 29, 35), ("Peacock", "epi", 35, 39),
          ("Duck", "fib", 39, 45), ("Duck", "epi", 45, 48),
          ("Zebra", "fib", 48, 52), ("Zebra", "epi", 52, 55),
          ("Emu", "fib", 55, 61), ("Emu", "epi", 61, 65)]

LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))

FAMILIES = {
    "cytosolic ribosome (RPL/RPS)": ("RPL", "RPS"),
    "mitoribosome (MRPL/MRPS)": ("MRPL", "MRPS"),
    "proteasome (PSM)": ("PSM",),
    "translation initiation (EIF)": ("EIF",),
    "heat shock (HSP)": ("HSP",),
    "collagen (COL)": ("COL",),
    "complex I (NDUF)": ("NDUF",),
    "complex IV (COX)": ("COX",),
    "ATP synthase (ATP5)": ("ATP5",),
}


def fit(y: np.ndarray) -> tuple[float, float]:
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3 or np.ptp(y[ok]) == 0:
        return np.nan, np.nan
    res = stats.linregress(LOGM[ok], y[ok])
    return float(res.rvalue), float(res.pvalue)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(CROSS / "counts.csv")["0"].map(normalize_symbol).values
    mat = pd.read_csv(CROSS / "file.csv", index_col=0)
    assert mat.shape == (3055, 65), mat.shape
    mat.index = genes
    mat.columns = COLUMNS
    mat = mat[~mat.index.duplicated()]
    print(f"matrix {mat.shape}, missing values {int(mat.isna().sum().sum())}", flush=True)

    per = {}
    for sp, tis, a, b in BLOCKS:
        per.setdefault(tis, {})[sp] = mat.iloc[:, a:b].mean(axis=1)
    fib = pd.DataFrame({sp: per["fib"][sp] for sp in SPECIES_ORDER})
    epi = pd.DataFrame({sp: per["epi"][sp] for sp in SPECIES_ORDER})

    ox_all = {normalize_symbol(g) for g in nuclear_oxphos()} - {""}
    ox = [g for g in fib.index if g in ox_all]
    report: dict = {"n_genes_matrix": int(len(fib)), "n_oxphos_in_matrix": len(ox),
                    "oxphos_genes": sorted(ox)}
    print(f"OXPHOS genes present: {len(ox)} -> {sorted(ox)}", flush=True)

    # ---------- level and rank ----------
    rows = []
    for tis, m in (("fib", fib), ("epi", epi)):
        pct = m.rank(axis=0, pct=True)
        lvl = m.loc[ox].mean(axis=0).values
        rk = pct.loc[ox].mean(axis=0).values
        for nm, v in (("mean level", lvl), ("mean percentile rank", rk)):
            r, p = fit(v)
            rows.append(dict(tissue=tis, metric=nm, r=r, p=p, values=np.round(v, 4).tolist()))
            print(f"  {tis} {nm:22s} r={r:+.3f} p={p:.4f}  {np.round(v,3).tolist()}")
    diff = fib.loc[ox].mean(axis=0).values - epi.loc[ox].mean(axis=0).values
    rd, pd_ = fit(diff)
    rows.append(dict(tissue="fib-epi", metric="mean level", r=rd, p=pd_,
                     values=np.round(diff, 4).tolist()))
    print(f"  fib-epi mean level          r={rd:+.3f} p={pd_:.4f}")
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "oxphos_vs_mass.tsv", sep="\t", index=False)
    report["oxphos_vs_mass"] = rows

    # ---------- matched null on the fixed ortholog set ----------
    chick = fib["Chicken"]
    prof = pd.DataFrame({"abund": chick})
    prof["decile"] = pd.qcut(prof["abund"].rank(method="first"), 10, labels=False)
    tgt = prof.reindex(ox).dropna()
    pool = prof.drop(index=ox, errors="ignore")
    obs_r, obs_p = fit(fib.loc[ox].mean(axis=0).values)
    nulls = []
    for _ in range(N_NULL):
        out = []
        for d, grp in tgt.groupby("decile", observed=True):
            cand = pool.index.values[pool["decile"].values == d]
            if len(cand) == 0:
                out = None
                break
            out.append(RNG.choice(cand, size=len(grp), replace=len(cand) < len(grp)))
        if out is None:
            continue
        r, _ = fit(fib.loc[np.concatenate(out)].mean(axis=0).values)
        if np.isfinite(r):
            nulls.append(r)
    nulls = np.array(nulls)
    emp = float((nulls <= obs_r).mean())
    report["matched_null"] = dict(observed_r=obs_r, naive_p=obs_p, empirical_p=emp,
                                  n_null=int(len(nulls)), null_median_r=float(np.median(nulls)))
    print(f"\nmatched null (abundance decile): obs r={obs_r:+.3f} naive p={obs_p:.4f} "
          f"null median={np.median(nulls):+.3f} empirical p={emp:.4f}")

    # ---------- gene family sweep ----------
    frows = []
    sets = {"OXPHOS": ox}
    for nm in ("TCA", "GLYCOLYSIS", "MITO_BIOGENESIS"):
        sets[nm] = [g for g in fib.index if g in {normalize_symbol(x) for x in load_geneset(nm)}]
    for label, pref in FAMILIES.items():
        sets[label] = [g for g in fib.index if str(g).startswith(pref)]
    for label, gl in sets.items():
        if len(gl) < 5:
            print(f"  {label:30s} n={len(gl)} skipped")
            continue
        r, p = fit(fib.loc[gl].mean(axis=0).values)
        re_, pe_ = fit(epi.loc[gl].mean(axis=0).values)
        frows.append(dict(gene_set=label, n_genes=len(gl), fib_r=r, fib_p=p,
                          epi_r=re_, epi_p=pe_))
        print(f"  {label:30s} n={len(gl):4d} fib r={r:+.3f} p={p:.4f}  epi r={re_:+.3f}")
    fam = pd.DataFrame(frows).sort_values("fib_r")
    fam.to_csv(OUT / "gene_family_sweep.tsv", sep="\t", index=False)
    report["families"] = frows

    # ---------- per-gene slopes ----------
    slopes = []
    for g, vals in zip(fib.index, fib.values):
        r, p = fit(vals)
        sl = stats.linregress(LOGM, vals).slope
        slopes.append(dict(symbol=g, slope=sl, r=r, p=p, is_oxphos=g in set(ox)))
    gt = pd.DataFrame(slopes).sort_values("slope")
    gt.to_csv(OUT / "per_gene_slopes.tsv", sep="\t", index=False)
    mw = stats.mannwhitneyu(gt.loc[gt.is_oxphos, "slope"], gt.loc[~gt.is_oxphos, "slope"],
                            alternative="less")
    report["per_gene"] = dict(oxphos_median_slope=float(gt.loc[gt.is_oxphos, "slope"].median()),
                              other_median_slope=float(gt.loc[~gt.is_oxphos, "slope"].median()),
                              mannwhitney_p=float(mw.pvalue))
    print(f"\nper-gene: OXPHOS median slope {gt.loc[gt.is_oxphos,'slope'].median():+.4f} vs "
          f"rest {gt.loc[~gt.is_oxphos,'slope'].median():+.4f}, MW p={mw.pvalue:.4f}")

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    ax = axes[0]
    lvl_f = fib.loc[ox].mean(axis=0).values
    lvl_e = epi.loc[ox].mean(axis=0).values
    ax.scatter(LOGM, lvl_f, s=70, color="#c1121f", label="fibroblasts", zorder=3)
    ax.scatter(LOGM, lvl_e, s=70, color="#2a9d8f", label="epiblast", zorder=3)
    xs = np.linspace(LOGM.min(), LOGM.max(), 10)
    for v, col in ((lvl_f, "#c1121f"), (lvl_e, "#2a9d8f")):
        s = stats.linregress(LOGM, v)
        ax.plot(xs, s.intercept + s.slope * xs, color=col, lw=1.1, alpha=0.5)
    for sp, x, y in zip(SPECIES_ORDER, LOGM, lvl_f):
        ax.annotate(sp, (x, y), fontsize=7, xytext=(3, 4), textcoords="offset points")
    rf, pf = fit(lvl_f)
    ax.set(xlabel="log10 body mass (g)", ylabel=f"mean of {len(ox)} OXPHOS orthologs (log units)",
           title=f"colleagues' 1-to-1 ortholog matrix\nfib r={rf:+.3f} p={pf:.3f}")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    ax.hist(nulls, bins=45, color="#adb5bd", edgecolor="none")
    ax.axvline(obs_r, color="#c1121f", lw=2)
    ax.axvline(float(np.median(nulls)), color="#264653", lw=1.2, ls="--")
    ax.set(xlabel="r vs log10 mass", ylabel="matched null gene sets",
           title=f"abundance-matched null\nempirical p={emp:.3f}")

    ax = axes[2]
    f2 = fam.dropna(subset=["fib_r"])
    cols = ["#c1121f" if s in ("OXPHOS",) or s.startswith(("complex", "ATP")) else "#457b9d"
            for s in f2.gene_set]
    ax.barh(range(len(f2)), f2.fib_r.values, color=cols)
    ax.set_yticks(range(len(f2)))
    ax.set_yticklabels(f2.gene_set, fontsize=7)
    ax.axvline(0, color="#333", lw=1)
    ax.set(xlabel="r vs log10 mass (fibroblasts)", title="gene families on the same matrix")

    fig.suptitle("OXPHOS vs body mass on the colleagues' cross-birds ortholog matrix", fontsize=12)
    for p in (FIG / "fig_colleagues_matrix_check.png", OUT / "fig_colleagues_matrix_check.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "summary.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
