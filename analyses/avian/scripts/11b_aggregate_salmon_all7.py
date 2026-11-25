#!/usr/bin/env python3
"""Aggregate salmon→chicken quants (7 spp) into gene TPM matrix; OXPHOS coverage vs Ruslan 3055."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
WORKDIR = Path("/path/to/storage/kleiber_birds_gemcat")
RES = ROOT / "results" / "salmon_all7_chicken"
FIG = ROOT / "figures"
DATA = ROOT / "data"
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")

# Ruslan adult body mass (g) — same as PaperA
MASS_G = {
    "ZebraFinch": 11.93,
    "Quail": 96.6,
    "Chicken": 912.5,
    "Duck": 1121.0,
    "Peacock": 4093.75,
    "Turkey": 5811.0,
    "Emu": 36200.0,
}

SPECIES_ORDER = list(MASS_G.keys())


def parse_gene(name: str) -> str:
    n = str(name)
    if "|" in n:
        return n.split("|", 1)[0].upper()
    if "gene_symbol:" in n:
        return n.split("gene_symbol:")[1].split()[0].upper()
    return n.upper()


def load_modules() -> dict[str, list[str]]:
    out = {}
    for p in (DATA / "genesets").glob("*.txt"):
        out[p.stem] = [ln.strip().upper() for ln in p.read_text().splitlines() if ln.strip()]
    return out


def load_library_tpms() -> pd.DataFrame:
    """genes × libraries TPM (sum of isoforms)."""
    cols = {}
    for qdir in sorted((WORKDIR / "quant").glob("*")):
        qf = qdir / "quant.sf"
        if not qf.exists():
            continue
        sf = pd.read_csv(qf, sep="\t")
        genes = sf["Name"].map(parse_gene)
        s = pd.Series(sf["TPM"].astype(float).values, index=genes).groupby(level=0).sum()
        cols[qdir.name] = s
    if not cols:
        raise SystemExit(f"No quant.sf under {WORKDIR}/quant")
    mat = pd.DataFrame(cols).fillna(0.0)
    mat.index.name = "gene"
    return mat


def species_means(lib_mat: pd.DataFrame) -> pd.DataFrame:
    """Mean TPM across sex libraries per species."""
    sp_cols = {}
    for sp in SPECIES_ORDER:
        libs = [c for c in lib_mat.columns if c.startswith(sp + "_")]
        if not libs:
            continue
        sp_cols[sp] = lib_mat[libs].mean(axis=1)
    return pd.DataFrame(sp_cols)


def z_module(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Gene-wise z across species, then mean — only genes present with TPM>0 in ≥3 spp."""
    g = [x for x in genes if x in expr.index]
    sub = expr.loc[g].astype(float)
    # keep genes detected in ≥3 species
    keep = (sub > 0).sum(axis=1) >= 3
    sub = sub.loc[keep]
    if sub.empty:
        return pd.Series({c: np.nan for c in expr.columns})
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    modules = load_modules()

    lib = load_library_tpms()
    lib.to_csv(RES / "gene_tpm_by_library.tsv", sep="\t")
    sp = species_means(lib)
    sp.to_csv(RES / "gene_tpm_by_species.tsv", sep="\t")

    present_libs = sorted(lib.columns.tolist())
    present_spp = [s for s in SPECIES_ORDER if s in sp.columns]

    # Coverage vs Ruslan matrix
    ruslan_genes = set(
        pd.read_csv(RUSLAN / "cross-birds/counts.csv").iloc[:, 1].astype(str).str.upper()
    )
    salmon_genes = set(sp.index.astype(str))
    # genes with TPM>0 in ≥1 / ≥5 / all present spp
    det1 = set(sp.index[(sp > 0).sum(axis=1) >= 1])
    det5 = set(sp.index[(sp > 0).sum(axis=1) >= min(5, len(present_spp))])
    det_all = set(sp.index[(sp > 0).sum(axis=1) >= len(present_spp)]) if present_spp else set()

    cov_rows = []
    for mod, genes in modules.items():
        gset = set(genes)
        cov_rows.append(
            {
                "module": mod,
                "list_n": len(gset),
                "in_ruslan_3055": len(gset & ruslan_genes),
                "in_salmon_index_symbols": len(gset & salmon_genes),
                "detected_ge1_spp": len(gset & det1),
                "detected_ge5_spp": len(gset & det5),
                "detected_all_present_spp": len(gset & det_all),
                "gained_vs_ruslan_ge1": len((gset & det1) - ruslan_genes),
            }
        )
    cov = pd.DataFrame(cov_rows)
    cov.to_csv(RES / "module_coverage_salmon_vs_ruslan.tsv", sep="\t", index=False)

    ox = set(modules["OXPHOS"])
    ox_present = sorted(ox & det1)
    ox_gained = sorted((ox & det1) - ruslan_genes)
    pd.Series(ox_present, name="gene").to_csv(RES / "OXPHOS_detected_salmon.tsv", sep="\t", index=False)
    pd.Series(ox_gained, name="gene").to_csv(RES / "OXPHOS_gained_vs_ruslan3055.tsv", sep="\t", index=False)

    # Module scores + OLS vs mass (only if ≥5 species done)
    score_rows = []
    fit_rows = []
    for mod, genes in modules.items():
        sc = z_module(sp[present_spp], genes) if present_spp else pd.Series(dtype=float)
        for s in present_spp:
            score_rows.append(
                {
                    "module": mod,
                    "species": s,
                    "mass_g": MASS_G[s],
                    "zmodule": float(sc.get(s, np.nan)),
                    "n_genes_used": int(((sp.loc[[g for g in genes if g in sp.index], s] > 0).sum()) if any(g in sp.index for g in genes) else 0),
                }
            )
        if len(present_spp) >= 5:
            y = sc.loc[present_spp].astype(float).values
            x = np.log10([MASS_G[s] for s in present_spp])
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() >= 5:
                slope, icept, r, p, se = stats.linregress(x[m], y[m])
                fit_rows.append(
                    {
                        "module": mod,
                        "n_spp": int(m.sum()),
                        "slope": float(slope),
                        "r": float(r),
                        "p": float(p),
                        "se": float(se),
                    }
                )
    scores = pd.DataFrame(score_rows)
    scores.to_csv(RES / "module_zscores_by_species.tsv", sep="\t", index=False)
    fits = pd.DataFrame(fit_rows)
    if len(fits):
        fits.to_csv(RES / "module_zscores_vs_mass_ols.tsv", sep="\t", index=False)

    # Figure: coverage bars + OXPHOS~mass if ready
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    ax = axes[0]
    sub = cov[cov.module == "OXPHOS"].iloc[0]
    labels = ["list", "Ruslan\n3055", "salmon\n≥1 spp", "salmon\n≥5 spp"]
    vals = [sub.list_n, sub.in_ruslan_3055, sub.detected_ge1_spp, sub.detected_ge5_spp]
    ax.bar(labels, vals, color=["#bbb", "#e76f51", "#1d3557", "#2a9d8f"])
    ax.set_ylabel("OXPHOS genes")
    ax.set_title(f"OXPHOS coverage (libs done: {len(present_libs)})")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.5, str(int(v)), ha="center", fontsize=9)

    ax = axes[1]
    if len(fits) and "OXPHOS" in set(fits.module):
        f = fits[fits.module == "OXPHOS"].iloc[0]
        subsc = scores[scores.module == "OXPHOS"].dropna(subset=["zmodule"])
        x = np.log10(subsc.mass_g)
        y = subsc.zmodule
        ax.scatter(x, y, s=70, c="#1d3557")
        for _, r in subsc.iterrows():
            ax.text(np.log10(r.mass_g), r.zmodule, r.species, fontsize=7)
        xs = np.linspace(x.min(), x.max(), 40)
        ax.plot(xs, f.slope * (xs - x.mean()) + y.mean(), color="#e76f51", lw=2)
        ax.set_title(f"OXPHOS z-module ~ mass\nslope={f.slope:.3f} p={f.p:.3g} n={f.n_spp}")
    else:
        ax.text(0.5, 0.5, f"Waiting for ≥5 species\n(now {len(present_spp)}: {present_spp})", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    ax.set_xlabel("log10 mass (g)")
    ax.set_ylabel("OXPHOS z-module")
    fig.suptitle("salmon → chicken cDNA (GSE213253 fib High 1F+1M)", fontsize=11)
    fig.savefig(FIG / "fig_salmon_all7_oxphos_coverage.png", dpi=160)
    fig.savefig(RES / "fig_salmon_all7_oxphos_coverage.png", dpi=160)
    plt.close(fig)

    summary = {
        "n_libraries_quanted": len(present_libs),
        "libraries": present_libs,
        "n_species": len(present_spp),
        "species": present_spp,
        "n_genes_in_matrix": int(lib.shape[0]),
        "oxphos_coverage": cov[cov.module == "OXPHOS"].iloc[0].to_dict(),
        "fits": fits.to_dict(orient="records") if len(fits) else [],
        "status": "complete" if len(present_spp) == 7 else "partial",
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(cov.to_string(index=False))


if __name__ == "__main__":
    main()
