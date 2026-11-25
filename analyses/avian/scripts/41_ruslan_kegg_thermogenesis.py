#!/usr/bin/env python3
"""Ruslan metric on his OXPHOS list and on KEGG thermogenesis (hsa04714).

Same score as birds - ESC.R: mean(colSums) of the gene set, then
log10(score) ~ log10(mass). RNA / Wagner / GSE167580 use the TPM-sum
equivalent (set sum; library is 1e6 so this is fraction × 1e6).
"""
from __future__ import annotations

import importlib.util
import re
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    EXP,
    FIG,
    PREFERRED_ALIAS,
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    load_geo_rnaseq,
    normalize_symbol,
)

OUT = RES / "ruslan_kegg_thermo"
GS = DATA / "genesets"
LINE = "#3366FF"
POINT = "#1a1a1a"

MT_ETC = {
    "ATP6", "ATP8", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "COX1", "COX2", "COX3", "CYC1",
}
ETC_PREFIX = ("NDUF", "COX", "UQCR", "ATP5", "SDH", "COA", "MT-")


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ruslan_oxphos_133() -> list[str]:
    text = (RUSLAN / "birds - ESC.R").read_text(errors="replace")
    m = re.search(r"^oxphos_genes <- c\((.*)\)\s*$", text, re.M)
    return [x.strip().strip('"') for x in re.findall(r'"([^"]+)"', m.group(1))]


def ruslan_thermo_229() -> list[str]:
    """Second oxphos_genes <-c( ... ) in birds - ESC.R: his KEGG-like list."""
    text = (RUSLAN / "birds - ESC.R").read_text(errors="replace")
    m = re.search(r"^oxphos_genes <-c\(\n(.*?)\n\)\s*$", text, re.M | re.S)
    return [x.strip().strip('"') for x in re.findall(r'"([^"]+)"', m.group(1))]


def kegg_04714() -> list[str]:
    src = Path("/tmp/kegg_hsa04714_get.txt")
    dest = GS / "kegg_hsa04714_get.txt"
    if src.exists():
        shutil.copy(src, dest)
    text = dest.read_text()
    genes = []
    in_gene = False
    for ln in text.splitlines():
        if ln.startswith("GENE"):
            in_gene = True
            ln = ln[4:]
        elif in_gene and ln[:12].strip() and not ln.startswith(" "):
            break
        if not in_gene:
            continue
        m = re.match(r"\s*\d+\s+(\S+);", ln)
        if m:
            genes.append(m.group(1))
    return genes


def is_etc(sym: str) -> bool:
    u = str(sym).upper()
    if u in MT_ETC or u.startswith("MT-"):
        return True
    return u.startswith(ETC_PREFIX)


def write_list(name: str, genes: list[str]) -> Path:
    path = GS / name
    path.write_text("\n".join(sorted(set(genes))) + "\n")
    return path


def expand(genes: set[str]) -> set[str]:
    """Match names as written, plus ATP5 aliases and MT- / COX1 variants."""
    out = {g.upper() for g in genes}
    rev = {}
    for old, new in PREFERRED_ALIAS.items():
        rev.setdefault(new, set()).add(old)
    extra = set()
    for g in list(out):
        extra.add(PREFERRED_ALIAS.get(g, g))
        extra.update(rev.get(g, ()))
        if g.startswith("MT-"):
            rest = g[3:]
            extra.add(rest)
            if rest.startswith("CO"):
                extra.add("COX" + rest[2:])
            if rest == "CYB":
                extra.add("CYTB")
        if g in MT_ETC:
            extra.add("MT-" + g)
            if g.startswith("COX"):
                extra.add("MT-CO" + g[3:])
            if g == "CYTB":
                extra.add("MT-CYB")
    return {x for x in out | extra if x}


def cage_score(path: Path, gene_col: str, genes: set[str]) -> dict:
    df = pd.read_csv(path, sep="\t")
    norm_all = [c for c in df.columns if str(c).startswith("norm.")]
    fib = [c for c in df.columns if str(c).startswith("norm.C")]
    df["_score"] = df[norm_all].apply(pd.to_numeric, errors="coerce").fillna(0).sum(axis=1)
    df = df.sort_values("_score", ascending=False)
    raw = df[gene_col].astype(str).str.upper()
    hit = raw.isin(genes)
    sub = df.loc[hit].copy()
    sub["_g"] = raw[hit]
    sub = sub.drop_duplicates("_g", keep="first")
    if sub.empty or not fib:
        return {"n_mapped": 0, "sum_tpm": float("nan"), "n_expr_med": 0}
    mat = sub[fib].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib_sums = mat.sum(axis=0)
    n_expr = [(mat[c] > 0).sum() for c in fib]
    return {
        "n_mapped": int(len(sub)),
        "sum_tpm": float(lib_sums.mean()),
        "n_expr_med": float(np.median(n_expr)),
    }


def theme_light(ax) -> None:
    ax.set_facecolor("white")
    ax.grid(True, which="major", color="#D0D0D0", lw=0.6)
    ax.grid(True, which="minor", color="#ECECEC", lw=0.4)
    for sp in ax.spines.values():
        sp.set_color("#B0B0B0")
        sp.set_linewidth(0.8)
    ax.tick_params(colors="#4A4A4A", labelsize=8)
    ax.xaxis.label.set_color("#4A4A4A")
    ax.yaxis.label.set_color("#4A4A4A")
    ax.set_axisbelow(True)


def draw_panel(ax, names, masses, scores, ylabel: str = "log10_fib") -> dict:
    x = np.log10(np.asarray(masses, float))
    y = np.log10(np.asarray(scores, float))
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    names = [n for n, k in zip(names, ok) if k]
    fit = stats.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 80)
    ax.plot(xs, fit.intercept + fit.slope * xs, color=LINE, lw=1.6, zorder=2)
    ax.scatter(x, y, s=22, c=POINT, zorder=3, linewidths=0)
    texts = [ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4) for name, xi, yi in zip(names, x, y)]
    adjust_text(
        texts, ax=ax, x=x, y=y, expand=(1.15, 1.25),
        force_text=(0.25, 0.4), force_static=(0.15, 0.2),
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
    )
    theme_light(ax)
    ax.set_xlabel("log10_body_mass")
    ax.set_ylabel(ylabel)
    ptxt = f"{fit.pvalue:.3f}" if fit.pvalue >= 0.001 else f"{fit.pvalue:.3g}"
    ax.set_title(f"cor. {fit.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
    return {"r": float(fit.rvalue), "p": float(fit.pvalue), "n": int(len(x)), "slope": float(fit.slope)}


def save_fig(fig, stem: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for folder in (FIG, OUT):
        fig.savefig(folder / f"{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
        fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fit_xy(masses, scores) -> dict:
    x = np.log10(np.asarray(masses, float))
    y = np.log10(np.asarray(scores, float))
    ok = np.isfinite(x) & np.isfinite(y)
    fit = stats.linregress(x[ok], y[ok])
    return {"r": float(fit.rvalue), "p": float(fit.pvalue), "n": int(ok.sum()), "slope": float(fit.slope)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    GS.mkdir(parents=True, exist_ok=True)

    ox = ruslan_oxphos_133()
    kegg = kegg_04714()
    sets = {
        "ruslan_oxphos_133": ox,
        "kegg_04714": kegg,
    }
    for name, genes in sets.items():
        write_list(f"{name}.txt", genes)
        print(f"SET {name:24s} n={len(set(genes))}")
    print("overlap 133∩kegg", len(set(g.upper() for g in ox) & set(g.upper() for g in kegg)))

    # Exact names as written. No homemade subsets, no alias expansion.
    match = {k: {g.upper() for g in v} for k, v in sets.items()}
    s36 = _load("s36", Path(__file__).resolve().parent / "36_gse167580_fib_oxphos_module.py")
    s37 = _load("s37", Path(__file__).resolve().parent / "37_wagner_sf_oxphos.py")

    rows = []

    # --- birds CAGE: Ruslan exact recipe ---
    cage_gene_col = {
        "Chicken": "geneId", "Duck": "geneId", "Quail": "geneId", "Zebra": "geneId",
        "Emu": "chicken_genes", "Peacock": "chicken_genes", "Turkey": "chicken_genes",
    }
    for set_name, genes in match.items():
        for sp in SPECIES_ORDER:
            meta = SPECIES[sp]
            sc = cage_score(EXP / meta["file"], cage_gene_col[sp], genes)
            rows.append(dict(
                panel="gse213253_cage", set=set_name, species=sp.lower(),
                mass_g=meta["mass_g"], **sc,
            ))
            print(f"CAGE {set_name:22s} {sp:8s} n={sc['n_mapped']:3d} sum={sc['sum_tpm']:.1f}", flush=True)

    # --- birds RNA: TPM-sum = fraction × 1e6 ---
    for set_name, genes in sets.items():
        wanted = {normalize_symbol(g) for g in genes if normalize_symbol(g)}
        wanted |= {PREFERRED_ALIAS.get(g, g) for g in wanted}
        for sp in SPECIES_ORDER:
            geo = load_geo_rnaseq(sp)
            gmat, tot = geo["gene_mat"], geo["lib_total"]
            idx = [g for g in gmat.index if g in wanted]
            lib_scores = []
            n_expr = []
            for col in gmat.columns:
                den = float(tot[col])
                num = float(gmat.loc[idx, col].sum()) if idx and den else 0.0
                lib_scores.append(1e6 * num / den if den else np.nan)
                n_expr.append(int((gmat.loc[idx, col] > 0).sum()) if idx else 0)
            rows.append(dict(
                panel="gse213253_rna", set=set_name, species=sp.lower(),
                mass_g=SPECIES[sp]["mass_g"],
                n_mapped=len(idx),
                sum_tpm=float(np.nanmean(lib_scores)),
                n_expr_med=float(np.median(n_expr)) if n_expr else 0,
            ))
            print(f"RNA  {set_name:22s} {sp:8s} n={len(idx):3d} sum={np.nanmean(lib_scores):.1f}", flush=True)

    # --- Wagner SF ---
    for set_name, genes in sets.items():
        wanted = {normalize_symbol(g) for g in genes if normalize_symbol(g)}
        wanted |= {PREFERRED_ALIAS.get(g, g) for g in wanted}
        per = {}
        nmap = {}
        for sp, sample, _srr in s37.SAMPLES:
            qf = s37.QUANT / sample / "quant.sf"
            if not qf.exists():
                continue
            t2g = pd.read_csv(s37.build_tx2gene(sp), sep="\t")
            q = pd.read_csv(qf, sep="\t")
            q = q.merge(t2g, left_on="Name", right_on="tx", how="left")
            q["symbol"] = q["symbol"].fillna("").map(lambda s: normalize_symbol(s) if s else "")
            if "raw_symbol" in q.columns:
                miss = q["symbol"].eq("")
                q.loc[miss, "symbol"] = q.loc[miss, "raw_symbol"].fillna("").map(normalize_symbol)
            lib = float(q["TPM"].sum())
            by = q.groupby("symbol", observed=True)["TPM"].sum()
            by = by[by.index != ""]
            present = [g for g in wanted if g in by.index]
            num = float(by.reindex(present).fillna(0).sum())
            per.setdefault(sp, []).append(1e6 * num / lib if lib else np.nan)
            nmap.setdefault(sp, []).append(len(present))
        for sp, vals in per.items():
            rows.append(dict(
                panel="wagner_sf", set=set_name, species=sp.replace("_", " "),
                mass_g=s37.MASS_G[sp],
                n_mapped=int(np.median(nmap[sp])),
                sum_tpm=float(np.nanmean(vals)),
                n_expr_med=int(np.median(nmap[sp])),
            ))
            print(f"WAG  {set_name:22s} {sp:12s} n={int(np.median(nmap[sp])):3d} sum={np.nanmean(vals):.1f}", flush=True)

    # --- GSE167580 primary-like ---
    all_wanted = set()
    for genes in sets.values():
        all_wanted |= {normalize_symbol(g) for g in genes if normalize_symbol(g)}
    maps = {}
    for sp in s36.FILES:
        maps[sp], _how = s36.id2symbol(sp, all_wanted)
    for set_name, genes in sets.items():
        wanted = {normalize_symbol(g) for g in genes if normalize_symbol(g)}
        for sp, fn in s36.FILES.items():
            df = pd.read_csv(s36.FIB / fn, sep="\t")
            scored = s36.score_table(df, maps[sp], wanted)
            vals, ns = [], []
            for col, rec in scored.items():
                if s36.SAMPLE_CLASS.get((sp, col), "unknown") != "primary_like":
                    continue
                vals.append(1e6 * rec["frac"])
                ns.append(rec["module_n_mapped"])
            if not vals:
                continue
            rows.append(dict(
                panel="gse167580", set=set_name, species=sp,
                mass_g=s36.MASS_G[sp],
                n_mapped=int(np.median(ns)),
                sum_tpm=float(np.nanmean(vals)),
                n_expr_med=int(np.median(ns)),
            ))
            print(f"GSE  {set_name:22s} {sp:8s} n={int(np.median(ns)):3d} sum={np.nanmean(vals):.1f}", flush=True)

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "species_scores.tsv", sep="\t", index=False)

    fit_rows = []
    for (panel, set_name), g in tab.groupby(["panel", "set"]):
        g = g[np.isfinite(g["sum_tpm"]) & (g["sum_tpm"] > 0)]
        if panel == "gse167580":
            g = g[g.species != "chicken"]
        if len(g) < 4:
            continue
        fit = fit_xy(g["mass_g"], g["sum_tpm"])
        fit_rows.append(dict(panel=panel, set=set_name, **fit, n_species=len(g)))
        print(f"FIT {panel:16s} {set_name:22s} n={fit['n']} r={fit['r']:+.3f} p={fit['p']:.3g}")
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(OUT / "fits.tsv", sep="\t", index=False)

    # figures: one 2×2 per set (CAGE, RNA, Wagner, GSE)
    titles = {
        "ruslan_oxphos_133": "Ruslan OXPHOS n=133",
        "kegg_04714": "KEGG hsa04714",
    }
    panel_order = [
        ("gse213253_cage", "birds CAGE"),
        ("gse213253_rna", "birds RNA"),
        ("wagner_sf", "Wagner SF"),
        ("gse167580", "GSE167580 mammals"),
    ]
    for set_name, title in titles.items():
        fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.6))
        for ax, (panel, lab) in zip(axes.ravel(), panel_order):
            g = tab[(tab.set == set_name) & (tab.panel == panel)]
            if panel == "gse167580":
                g = g[g.species != "chicken"]
            g = g[np.isfinite(g["sum_tpm"]) & (g["sum_tpm"] > 0)]
            fit = draw_panel(ax, g.species.tolist(), g.mass_g.tolist(), g.sum_tpm.tolist())
            ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
            ax.set_title(f"{lab}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
        fig.suptitle(title, fontsize=12, y=0.98)
        fig.subplots_adjust(wspace=0.32, hspace=0.40, left=0.08, right=0.99, top=0.92, bottom=0.08)
        save_fig(fig, f"fig_ruslan_kegg_{set_name}")

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.6))
    for ax, set_name in zip(axes, titles):
        g = tab[(tab.set == set_name) & (tab.panel == "gse213253_cage")]
        fit = draw_panel(ax, g.species.tolist(), g.mass_g.tolist(), g.sum_tpm.tolist())
        ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
        ax.set_title(f"{titles[set_name]}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=9, pad=7)
    fig.subplots_adjust(wspace=0.28, left=0.08, right=0.99, top=0.86, bottom=0.16)
    save_fig(fig, "fig_ruslan_kegg_cage_two")

    print("wrote", OUT)


if __name__ == "__main__":
    main()
