#!/usr/bin/env python3
"""Recompute fibroblast module allometry with n=7 after salmon→chicken quant for Emu/Turkey/Peacock."""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
WORKDIR = Path("/path/to/storage/kleiber_birds_gemcat")
RES = ROOT / "results"
FIG = ROOT / "figures"
DATA = ROOT / "data"

MASS = {
    "ZebraFinch": 0.0157,
    "Quail": 0.12,
    "Duck": 1.1,
    "Chicken": 2.0,
    "Peacock": 4.5,
    "Turkey": 8.0,
    "Emu": 35.0,
}
ORG2SP = {
    "Taeniopygia guttata": "ZebraFinch",
    "Coturnix japonica": "Quail",
    "Anas platyrhynchos": "Duck",
    "Gallus gallus": "Chicken",
    "Pavo cristatus": "Peacock",
    "Meleagris gallopavo": "Turkey",
    "Dromaius novaehollandiae": "Emu",
}


def load_modules():
    out = {}
    for p in (DATA / "genesets").glob("*.txt"):
        out[p.stem] = [ln.strip().upper() for ln in p.read_text().splitlines() if ln.strip()]
    return out


def module_score(expr: pd.Series, genes, min_genes=5):
    present = [g for g in genes if g in expr.index]
    if len(present) < min_genes:
        return np.nan, len(present)
    return float(np.log2(expr.loc[present].astype(float) + 1).mean()), len(present)


def load_geo_symbol_species():
    """Existing GEO tables with gene symbols."""
    files = {
        "ZebraFinch": "GSE213253_Zebrafinch_RNAseq.xlsx",
        "Quail": "GSE213253_Quail_RNAseq.xlsx",
        "Duck": "GSE213253_Duck_RNAseq.xlsx",
        "Chicken": "GSE213253_Chicken_RNAseq.xlsx",
    }
    expr = {}
    for sp, fn in files.items():
        df = pd.read_excel(DATA / fn)
        norms = [c for c in df.columns if str(c).startswith("norm.")]
        g = df["gene"].astype(str).str.upper()
        s = df[norms].apply(pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
        expr[sp] = pd.Series(s.values, index=g.values).groupby(level=0).mean()
    return expr


def load_salmon_quant():
    """Mean TPM per species from salmon quant dirs named SPECIES_SEX."""
    quant_root = WORKDIR / "quant"
    by_sp = {}
    for d in sorted(quant_root.glob("*")):
        if not d.is_dir():
            continue
        qf = d / "quant.sf"
        if not qf.exists():
            continue
        # dirname like Emu_female
        parts = d.name.split("_")
        if len(parts) < 2:
            continue
        sp = parts[0]
        sf = pd.read_csv(qf, sep="\t")
        # Name column may be ENSGALT... or gene symbol|...
        names = sf["Name"].astype(str)
        # Ensembl cdna headers often gene_symbol in desc; our fasta uses gene=
        genes = []
        for n in names:
            if "gene_symbol:" in n:
                genes.append(n.split("gene_symbol:")[1].split()[0].upper())
            elif "|" in n:
                genes.append(n.split("|")[0].upper())
            else:
                genes.append(n.upper())
        s = pd.Series(sf["TPM"].values, index=genes).groupby(level=0).sum()
        by_sp.setdefault(sp, []).append(s)
    out = {}
    for sp, lst in by_sp.items():
        mat = pd.concat(lst, axis=1).fillna(0)
        out[sp] = mat.mean(axis=1)
    return out


def fit(mass, score):
    x = np.log10(np.asarray(mass, float))
    y = np.asarray(score, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return dict(n=len(x), slope=np.nan, r=np.nan, p=np.nan, se=np.nan)
    slope, intercept, r, p, se = stats.linregress(x, y)
    return dict(n=int(len(x)), slope=float(slope), intercept=float(intercept), r=float(r), p=float(p), se=float(se))


def main():
    modules = load_modules()
    geo = load_geo_symbol_species()
    sal = load_salmon_quant()
    # merge: prefer salmon for recovered spp; keep GEO for the four
    expr = dict(geo)
    for sp in ("Emu", "Turkey", "Peacock"):
        if sp not in sal:
            raise SystemExit(f"Missing salmon quant for {sp}. Check {WORKDIR}/quant")
        expr[sp] = sal[sp]

    # optional: also overwrite all if salmon has them
    note = {"geo_species": list(geo), "salmon_species": list(sal), "final_species": list(expr)}

    rows = []
    for mod, genes in modules.items():
        for sp, e in expr.items():
            sc, n = module_score(e, genes)
            rows.append({"module": mod, "species": sp, "mass_kg": MASS[sp], "score_log2": sc, "n_genes": n})
    scores = pd.DataFrame(rows)
    scores.to_csv(RES / "module_scores_by_species_n7.tsv", sep="\t", index=False)

    fits = []
    for mod in modules:
        sub = scores[scores.module == mod].dropna(subset=["score_log2"])
        f = fit(sub.mass_kg, sub.score_log2)
        fits.append(
            {
                "module": mod,
                **f,
                "abs_err_m13": abs(f["slope"] + 1 / 3) if np.isfinite(f["slope"]) else np.nan,
                "abs_err_m14": abs(f["slope"] + 1 / 4) if np.isfinite(f["slope"]) else np.nan,
            }
        )
    fits_df = pd.DataFrame(fits)
    fits_df.to_csv(RES / "allometry_ols_fits_n7.tsv", sep="\t", index=False)
    print(fits_df.to_string(index=False))

    fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
    for ax, mod in zip(axes.ravel(), modules):
        sub = scores[scores.module == mod].dropna(subset=["score_log2"])
        x = np.log10(sub.mass_kg.values)
        y = sub.score_log2.values
        ax.scatter(x, y, s=60)
        for sp, xi, yi in zip(sub.species, x, y):
            ax.text(xi, yi, sp, fontsize=7)
        if len(x) >= 3 and np.all(np.isfinite(y)):
            slope, intercept, r, p, _ = stats.linregress(x, y)
            xs = np.linspace(x.min() - 0.1, x.max() + 0.1, 40)
            ax.plot(xs, intercept + slope * xs, color="black", lw=1.5)
            y0, x0 = np.mean(y), np.mean(x)
            ax.plot(xs, y0 - (1 / 3) * (xs - x0), "--", label="M^{-1/3}")
            ax.plot(xs, y0 - (1 / 4) * (xs - x0), "--", label="M^{-1/4}")
            ax.set_title(f"{mod}: b={slope:.2f}, r={r:.2f}, p={p:.3f}")
            ax.legend(fontsize=7)
        ax.set_xlabel("log10 mass (kg)")
        ax.set_ylabel("module score")
    fig.suptitle("Fibroblast modules vs mass — n=7 (GEO symbols + salmon→chicken for Emu/Turkey/Peacock)")
    fig.savefig(FIG / "fig_A2_fibroblast_module_allometry_n7.png", dpi=160)
    plt.close(fig)

    with open(RES / "recover3_notes.json", "w") as f:
        json.dump(note, f, indent=2)
    print("Wrote n=7 fits and figure.")


if __name__ == "__main__":
    main()
