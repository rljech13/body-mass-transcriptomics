#!/usr/bin/env python3
"""Ruslan metric on Ma 2016 and GSE105160 mock fibroblasts.

Two KEGG pathways only:
  hsa00190 Oxidative phosphorylation
  hsa04714 Thermogenesis
RNA score = 1e6 × set-counts / library size, missing gene = 0.
Species score = mean of libraries. Fit is log10(score) ~ log10(mass);
the figure swaps axes: x = log10(score), y = log10(mass).
Panels are never pooled.
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from adjustText import adjust_text
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import DATA, FIG, PREFERRED_ALIAS, RES, normalize_symbol  # noqa: E402

ROOT = DATA / "file_panels"
MA_XLSX = ROOT / "ma2016" / "elife19130_supp1.xlsx"
MA_FIG1 = ROOT / "ma2016" / "elife19130_fig1.xlsx"
GSE = ROOT / "gse105160"
COUNTS = GSE / "counts_mock"
SOFT = GSE / "GSE105160_family.soft.gz"
OUT = RES / "file_panels_ruslan"
LINE = "#3366FF"
POINT = "#1a1a1a"

# AnAge Adult weight (g), dump 2023-07-03. P. anubis and P. abelii have no row.
# Pongo: AnAge lists Pongo pygmaeus only. Papio: Strum 1991 wild ♂24 + ♀14.7 kg mean.
GSE_MASS = {
    "Mus musculus": (20.5, "AnAge Adult weight"),
    "Saimiri sciureus": (925.0, "AnAge Adult weight"),
    "Macaca nemestrina": (7912.5, "AnAge Adult weight"),
    "Macaca mulatta": (8235.0, "AnAge Adult weight"),
    "Papio anubis": (19350.0, "no AnAge row; Strum 1991 wild 24+14.7 kg / 2 (PIN)"),
    "Pan paniscus": (39925.0, "AnAge Adult weight"),
    "Pan troglodytes": (44983.5, "AnAge Adult weight"),
    "Homo sapiens": (62035.0, "AnAge Adult weight"),
    "Pongo abelii": (64475.0, "AnAge has Pongo pygmaeus only; same AW"),
    "Gorilla gorilla": (139842.0, "AnAge Adult weight"),
}
GSE_LABEL = {
    "Mus musculus": "mouse",
    "Saimiri sciureus": "squirrel monkey",
    "Macaca nemestrina": "pigtail macaque",
    "Macaca mulatta": "rhesus",
    "Papio anubis": "olive baboon",
    "Pan paniscus": "bonobo",
    "Pan troglodytes": "chimp",
    "Homo sapiens": "human",
    "Pongo abelii": "orangutan",
    "Gorilla gorilla": "gorilla",
}
MA_LABEL = {
    "Shrew": "shrew",
    "BigBrownBat": "big brown bat",
    "LittleBrownBat": "little brown bat",
    "GuineaPig": "guinea pig",
    "Porcupine": "porcupine",
    "Chinchilla": "chinchilla",
    "Chipmunk": "chipmunk",
    "RedSquirrel": "red squirrel",
    "FoxSquirrel": "fox squirrel",
    "Beaver": "beaver",
    "AfricanGrassRat": "grass rat",
    "MeadowVole": "meadow vole",
    "CottonRat": "cotton rat",
    "WhiteFootedMouse": "white-footed mouse",
    "DeerMouse": "deer mouse",
}
MA_PAPER_NAME = {
    "Shrew": "Shrew",
    "BigBrownBat": "Big brown bat",
    "LittleBrownBat": "Little brown bat",
    "GuineaPig": "Guinea pig",
    "Porcupine": "Porcupine",
    "Chinchilla": "Chinchilla",
    "Chipmunk": "Chipmunk",
    "RedSquirrel": "Red squirrel",
    "FoxSquirrel": "Fox squirrel",
    "Beaver": "Beaver",
    "AfricanGrassRat": "African grass rat",
    "MeadowVole": "Meadow vole",
    "CottonRat": "Cotton rat",
    "WhiteFootedMouse": "White footed mouse",
    "DeerMouse": "Deer mouse",
}


# KEGG writes mt genes as ATP6/COX1/ND1; Ensembl/NCBI use MT-ATP6/MT-CO1/MT-ND1.
KEGG_MT = {
    "ATP6": "MT-ATP6", "ATP8": "MT-ATP8",
    "COX1": "MT-CO1", "COX2": "MT-CO2", "COX3": "MT-CO3",
    "CYTB": "MT-CYB",
    "ND1": "MT-ND1", "ND2": "MT-ND2", "ND3": "MT-ND3",
    "ND4": "MT-ND4", "ND4L": "MT-ND4L", "ND5": "MT-ND5", "ND6": "MT-ND6",
}
KEGG_MT_REV = {v: k for k, v in KEGG_MT.items()}


def load_sets() -> dict[str, set[str]]:
    ox = {ln.strip().upper() for ln in (DATA / "genesets" / "kegg_00190.txt").read_text().splitlines() if ln.strip()}
    th = {ln.strip().upper() for ln in (DATA / "genesets" / "kegg_04714.txt").read_text().splitlines() if ln.strip()}
    return {"kegg_00190": ox, "kegg_04714": th}


def wanted_symbols(genes: set[str]) -> set[str]:
    """KEGG names as written, ATP5 rename table, and MT- dual names from the same KEGG entries."""
    out = set()
    for g in genes:
        u = g.upper().strip()
        if not u:
            continue
        out.add(u)
        n = normalize_symbol(u)
        if n:
            out.add(n)
        if u in KEGG_MT:
            out.add(KEGG_MT[u])
        if u in KEGG_MT_REV:
            out.add(KEGG_MT_REV[u])
    return out


def mouse_ens_to_human() -> dict[str, str]:
    df = pd.read_csv(DATA / "ensembl_orthologs" / "mouse_to_human.tsv", sep="\t")
    one = df[df["Human homology type"].astype(str).eq("ortholog_one2one")]
    out: dict[str, str] = {}
    for mid, hname in zip(one["Gene stable ID"], one["Human gene name"]):
        sym = normalize_symbol(str(hname))
        if not sym:
            continue
        out[str(mid)] = sym
    return out


def human_ens_to_symbol() -> dict[str, str]:
    df = pd.read_csv(DATA / "ensembl_orthologs" / "human_ensg_symbol.tsv", sep="\t")
    out: dict[str, str] = {}
    for eid, name in zip(df["Gene stable ID"], df["Gene name"]):
        sym = normalize_symbol(str(name))
        if sym:
            out[str(eid)] = sym
    return out


def ma_paper_aw() -> dict[str, float]:
    wb = openpyxl.load_workbook(MA_FIG1, data_only=True)
    ws = wb["A. Species and Traits"]
    aw = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        if not row or not row[3] or row[4] is None:
            continue
        aw[str(row[3]).strip()] = float(row[4])
    return aw


def parse_soft_mock() -> pd.DataFrame:
    text = gzip.open(SOFT, "rt").read()
    blocks = re.split(r"\n\^SAMPLE = ", text)
    rows = []
    for b in blocks[1:]:
        lines = b.splitlines()
        gsm = lines[0].strip()
        title = org = ""
        chars = []
        for ln in lines:
            if ln.startswith("!Sample_title = "):
                title = ln.split(" = ", 1)[1]
            elif ln.startswith("!Sample_organism_ch1 = "):
                org = ln.split(" = ", 1)[1]
            elif ln.startswith("!Sample_characteristics_ch1 = "):
                chars.append(ln.split(" = ", 1)[1])
        treat = " ".join(c for c in chars if "treatment" in c)
        if "unmatched" in title.lower():
            continue
        if "mock" not in title.lower() and "mock" not in treat.lower():
            continue
        rows.append({"gsm": gsm, "title": title, "organism": org})
    return pd.DataFrame(rows)


def read_count_file(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t", compression="gzip")
    idcol = df.columns[0]
    valcol = df.columns[1]
    s = pd.to_numeric(df[valcol], errors="coerce").fillna(0.0)
    s.index = df[idcol].astype(str)
    return s.groupby(level=0).sum()


def score_counts(counts: pd.Series, id2sym: dict[str, str], wanted: set[str]) -> dict:
    mapped = counts.groupby(counts.index.map(lambda i: id2sym.get(i, ""))).sum()
    mapped = mapped[mapped.index != ""]
    lib = float(counts.sum())
    present = [g for g in wanted if g in mapped.index]
    num = float(mapped.reindex(present).fillna(0).sum()) if present else 0.0
    n_expr = int((mapped.reindex(present).fillna(0) > 0).sum()) if present else 0
    return {
        "lib_total": lib,
        "n_mapped": len(present),
        "n_expr": n_expr,
        "sum_tpm": (1e6 * num / lib) if lib else float("nan"),
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


def draw_panel(ax, names, masses, scores) -> dict:
    mass = np.log10(np.asarray(masses, float))
    score = np.log10(np.asarray(scores, float))
    ok = np.isfinite(mass) & np.isfinite(score)
    mass, score = mass[ok], score[ok]
    names = [n for n, k in zip(names, ok) if k]
    # Report the same r/p as the stored fit (score ~ mass). Draw that line with axes swapped.
    fit = stats.linregress(mass, score)
    x, y = score, mass
    mass_line = np.linspace(mass.min(), mass.max(), 80)
    ax.plot(fit.intercept + fit.slope * mass_line, mass_line, color=LINE, lw=1.6, zorder=2)
    ax.scatter(x, y, s=22, c=POINT, zorder=3, linewidths=0)
    texts = [ax.text(xi, yi, name, fontsize=7.5, color="#222222", zorder=4) for name, xi, yi in zip(names, x, y)]
    adjust_text(
        texts, ax=ax, x=x, y=y, expand=(1.15, 1.25),
        force_text=(0.25, 0.4), force_static=(0.15, 0.2),
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
    )
    theme_light(ax)
    ax.set_xlabel("log10_fib")
    ax.set_ylabel("log10_body_mass")
    return {"r": float(fit.rvalue), "p": float(fit.pvalue), "n": int(len(x)), "slope": float(fit.slope), "intercept": float(fit.intercept)}


def save_fig(fig, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    for folder in (FIG, OUT):
        fig.savefig(folder / f"{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
        fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def score_ma(sets: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    aw_paper = ma_paper_aw()
    id2h = mouse_ens_to_human()
    raw = pd.read_excel(MA_XLSX, sheet_name="A. raw counts")
    raw = raw.rename(columns={raw.columns[0]: "gene"})
    lib_rows = []
    for col in raw.columns[1:]:
        m = re.match(r"^(.+)\.\d", str(col))
        if not m:
            continue
        key = m.group(1)
        paper = MA_PAPER_NAME[key]
        counts = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
        counts.index = raw["gene"].astype(str)
        counts = counts.groupby(level=0).sum()
        for set_name, genes in sets.items():
            sc = score_counts(counts, id2h, wanted_symbols(genes))
            lib_rows.append(dict(
                panel="ma2016", set=set_name, sample=col, species=MA_LABEL[key],
                organism=paper, mass_g=aw_paper[paper], mass_source="Ma eLife fig1 AW",
                **sc,
            ))
    libs = pd.DataFrame(lib_rows)
    sp = libs.groupby(["panel", "set", "species", "organism", "mass_g", "mass_source"], as_index=False).agg(
        n_libs=("sample", "nunique"),
        n_mapped=("n_mapped", "median"),
        n_expr_med=("n_expr", "median"),
        sum_tpm=("sum_tpm", "mean"),
    )
    return libs, sp


def score_gse(sets: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    mock = parse_soft_mock()
    ens_h = human_ens_to_symbol()
    ens_m = mouse_ens_to_human()
    files = {re.search(r"GSM\d+", p.name).group(): p for p in COUNTS.glob("*")}
    lib_rows = []
    for rec in mock.itertuples(index=False):
        path = files.get(rec.gsm)
        if path is None:
            print("MISSING", rec.gsm, rec.organism)
            continue
        counts = read_count_file(path)
        id2sym = ens_m if rec.organism == "Mus musculus" else ens_h
        mass, src = GSE_MASS[rec.organism]
        wanted = {k: wanted_symbols(v) for k, v in sets.items()}
        for set_name, genes in wanted.items():
            sc = score_counts(counts, id2sym, genes)
            lib_rows.append(dict(
                panel="gse105160_mock", set=set_name, sample=rec.gsm, species=GSE_LABEL[rec.organism],
                organism=rec.organism, title=rec.title, mass_g=mass, mass_source=src,
                **sc,
            ))
    libs = pd.DataFrame(lib_rows)
    sp = libs.groupby(["panel", "set", "species", "organism", "mass_g", "mass_source"], as_index=False).agg(
        n_libs=("sample", "nunique"),
        n_mapped=("n_mapped", "median"),
        n_expr_med=("n_expr", "median"),
        sum_tpm=("sum_tpm", "mean"),
    )
    return libs, sp


def fit_xy(masses, scores) -> dict:
    x = np.log10(np.asarray(masses, float))
    y = np.log10(np.asarray(scores, float))
    ok = np.isfinite(x) & np.isfinite(y)
    fit = stats.linregress(x[ok], y[ok])
    return {"r": float(fit.rvalue), "p": float(fit.pvalue), "n": int(ok.sum()),
            "slope": float(fit.slope), "intercept": float(fit.intercept)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sets = load_sets()
    print("sets", {k: len(v) for k, v in sets.items()}, "overlap", len(sets["kegg_00190"] & sets["kegg_04714"]))

    ma_lib, ma_sp = score_ma(sets)
    gse_lib, gse_sp = score_gse(sets)
    libs = pd.concat([ma_lib, gse_lib], ignore_index=True)
    tab = pd.concat([ma_sp, gse_sp], ignore_index=True)
    libs.to_csv(OUT / "libraries.tsv", sep="\t", index=False)
    tab.to_csv(OUT / "species_scores.tsv", sep="\t", index=False)

    fit_rows = []
    for (panel, set_name), g in tab.groupby(["panel", "set"]):
        g = g[np.isfinite(g["sum_tpm"]) & (g["sum_tpm"] > 0)]
        fit = fit_xy(g["mass_g"], g["sum_tpm"])
        fit_rows.append(dict(panel=panel, set=set_name, **fit, n_species=len(g)))
        print(f"FIT {panel:16s} {set_name:22s} n={fit['n']} r={fit['r']:+.3f} p={fit['p']:.3g}")
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(OUT / "fits.tsv", sep="\t", index=False)

    titles = {"kegg_00190": "KEGG hsa00190 OXPHOS", "kegg_04714": "KEGG hsa04714 thermogenesis"}
    panels = [("ma2016", "Ma 2016 · abdominal skin"), ("gse105160_mock", "GSE105160 · dermal mock")]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2))
    for i, set_name in enumerate(titles):
        for j, (panel, lab) in enumerate(panels):
            ax = axes[i, j]
            g = tab[(tab.set == set_name) & (tab.panel == panel)].sort_values("mass_g")
            fit = draw_panel(ax, g.species.tolist(), g.mass_g.tolist(), g.sum_tpm.tolist())
            ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
            ax.set_title(f"{lab}\n{titles[set_name]}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=9, pad=6)
    fig.subplots_adjust(wspace=0.28, hspace=0.42, left=0.07, right=0.99, top=0.90, bottom=0.08)
    save_fig(fig, "fig_file_panels_ruslan_2x2")

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0))
    for ax, set_name in zip(axes, titles):
        g = tab[(tab.set == set_name) & (tab.panel == "ma2016")].sort_values("mass_g")
        fit = draw_panel(ax, g.species.tolist(), g.mass_g.tolist(), g.sum_tpm.tolist())
        ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
        ax.set_title(f"{titles[set_name]}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
    fig.suptitle("Ma 2016 · adult abdominal skin fibroblasts", fontsize=12, y=0.98)
    fig.subplots_adjust(wspace=0.28, left=0.08, right=0.99, top=0.84, bottom=0.14)
    save_fig(fig, "fig_file_panels_ma2016")

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0))
    for ax, set_name in zip(axes, titles):
        g = tab[(tab.set == set_name) & (tab.panel == "gse105160_mock")].sort_values("mass_g")
        fit = draw_panel(ax, g.species.tolist(), g.mass_g.tolist(), g.sum_tpm.tolist())
        ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
        ax.set_title(f"{titles[set_name]}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
    fig.suptitle("GSE105160 · primary dermal fibroblasts, mock only", fontsize=12, y=0.98)
    fig.subplots_adjust(wspace=0.28, left=0.08, right=0.99, top=0.84, bottom=0.14)
    save_fig(fig, "fig_file_panels_gse105160")

    print("wrote", OUT)
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
