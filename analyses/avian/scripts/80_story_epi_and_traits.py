#!/usr/bin/env python3
"""Story sets in epiblast, and the same sums vs Ruslan's other three traits.

Same homologous-CAGE assignment as scripts 75/79 (nearest chicken-protein CDS
≤ 2 kb, peak kept if fibroblast mean > 0). Epiblast is the mean of those same
peaks' norm.E libraries, not a new promoter call. Missing ≠ 0 at assignment;
a gene already in the fibroblast matrix may have epi = 0 at that promoter.

y = log10(sum TPM of the set). Predictors from birds-ESC.R species_traits:
log10 mass, incubation_d, max_longevity_y, log10 genome_size_gb. Tissues
and traits are separate panels. Do not glue sets or studies.
"""
from __future__ import annotations

import importlib.util
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
s74 = importlib.import_module("74_nineteen_via_cage")
s75 = importlib.import_module("75_honest_gene_mass_scan")
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    FIG,
    RES,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    peak_to_symbol_de,
    RUSLAN,
)

OUT = RES / "story_epi_traits"
M_FIB_PATH = RES / "honest_gene_scan" / "M_fib_nearest2kb.tsv"
THR = 1.0
SETS = [
    ("oxphos_16", "OXPHOS 16", "fig_story_oxphos_16"),
    ("folate_3", "folate 3", "fig_story_folate_3"),
    ("ribo_5", "ribo 5", "fig_story_ribo_5"),
    ("green_hubs_7", "green hubs 7", "fig_story_green_hubs_7"),
    ("green_sig_51", "green 51", "fig_story_green_51"),
]
# Ruslan tribble, birds-ESC.R ~5504. Peacock genome is the k-mer estimate.
TRAITS = {
    "Zebra": dict(mass_g=11.93, incubation_d=14.15, max_longevity_y=14.5, genome_size_gb=1.223),
    "Quail": dict(mass_g=96.6, incubation_d=17.0, max_longevity_y=14.58, genome_size_gb=1.379),
    "Chicken": dict(mass_g=912.5, incubation_d=20.5, max_longevity_y=30.0, genome_size_gb=1.223),
    "Duck": dict(mass_g=1121.0, incubation_d=26.5, max_longevity_y=29.1, genome_size_gb=1.428),
    "Peacock": dict(mass_g=4093.75, incubation_d=28.0, max_longevity_y=25.0, genome_size_gb=1.130),
    "Turkey": dict(mass_g=5811.0, incubation_d=28.0, max_longevity_y=13.0, genome_size_gb=1.281),
    "Emu": dict(mass_g=36200.0, incubation_d=56.0, max_longevity_y=16.6, genome_size_gb=1.594),
}
PREDICTORS = [
    ("mass", "log10_body_mass", lambda t: np.log10(t["mass_g"])),
    ("incubation", "incubation_d", lambda t: float(t["incubation_d"])),
    ("longevity", "max_longevity_y", lambda t: float(t["max_longevity_y"])),
    ("genome", "log10_genome_size", lambda t: np.log10(t["genome_size_gb"])),
]


def load40():
    spec = importlib.util.spec_from_file_location(
        "s40", Path(__file__).resolve().parent / "40_ruslan_style_scatters.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def genes_of(stem: str) -> list[str]:
    return [ln.strip() for ln in (DATA / "genesets" / f"story_{stem}.txt").read_text().splitlines() if ln.strip()]


def load_peaks_both(species: str) -> pd.DataFrame:
    """Same fibroblast filter as script 74, plus mean of E columns on those peaks."""
    d = pd.read_csv(s74.s62.EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
    fcols = lib_cols(d, "fib")
    ecols = lib_cols(d, "epi")
    fnum = d[fcols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    fnum = fnum.where(fnum >= THR, 0.0)
    if ecols:
        enum = d[ecols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        enum = enum.where(enum >= THR, 0.0)
        epi = enum.mean(axis=1)
    else:
        epi = pd.Series(0.0, index=d.index)
    d = d.assign(
        mid=(d["start"] + d["end"]) / 2.0,
        fib=fnum.mean(axis=1),
        epi=epi,
        seq=d["seqnames"].map(s74.contig),
        n_fib=len(fcols),
        n_epi=len(ecols),
    )
    return d[d["fib"] > 0].copy()


def assign_both(peaks: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
    """Identical catchment to script 75; also sums epiblast on the same peaks."""
    if loci.empty or peaks.empty:
        return pd.DataFrame(columns=["gene", "fib", "epi", "n_peaks", "dist"])
    index: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
    for (seq, st), g in loci.groupby(["seqname", "strand"], sort=False):
        g = g.sort_values("five")
        index[(seq, st)] = (g["five"].to_numpy(float), g["symbol"].tolist())

    acc: dict[str, list] = {}
    for r in peaks.itertuples():
        packed = index.get((r.seq, r.strand))
        if packed is None:
            continue
        fives, syms = packed
        mid = float(r.mid)
        i = int(np.searchsorted(fives, mid))
        cands: list[int] = []
        if i > 0:
            cands.append(i - 1)
        if i < len(fives):
            cands.append(i)
        if not cands:
            continue
        ds = [abs(mid - fives[j]) for j in cands]
        jbest = cands[int(np.argmin(ds))]
        d = abs(mid - fives[jbest])
        if d > s75.MAX_BP:
            continue
        if len(cands) == 2 and abs(ds[0] - ds[1]) < 1e-6:
            continue
        gene = syms[jbest]
        acc.setdefault(gene, []).append((float(r.fib), float(r.epi), d))

    rows = []
    for gene, hits in acc.items():
        fib = float(sum(h[0] for h in hits))
        if fib <= 0:
            continue
        rows.append(dict(
            gene=gene,
            fib=fib,
            epi=float(sum(h[1] for h in hits)),
            n_peaks=len(hits),
            dist=min(h[2] for h in hits),
        ))
    return pd.DataFrame(rows)


def scatter(ax, names, x, y, xlabel: str, ylabel: str, s40) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    fit = stats.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 80)
    ax.plot(xs, fit.intercept + fit.slope * xs, color=s40.LINE, lw=1.6, zorder=2)
    ax.scatter(x, y, s=22, c=s40.POINT, zorder=3, linewidths=0)
    texts = [
        ax.text(xi, yi, name, fontsize=8, color="#222222", zorder=4)
        for name, xi, yi in zip(names, x, y)
    ]
    s40.adjust_text(
        texts, ax=ax, x=x, y=y,
        expand=(1.15, 1.25), force_text=(0.25, 0.4), force_static=(0.15, 0.2),
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
    )
    s40.theme_light(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ptxt = f"{fit.pvalue:.3f}" if fit.pvalue >= 0.001 else f"{fit.pvalue:.3g}"
    ax.set_title(f"cor. {fit.rvalue:.3f}; p = {ptxt}", fontsize=11, pad=8)
    return {
        "r": float(fit.rvalue), "p": float(fit.pvalue),
        "slope": float(fit.slope), "intercept": float(fit.intercept),
        "n": int(len(names)),
    }


def save(fig, stem: str, s40) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for folder in (FIG, OUT):
        fig.savefig(folder / f"{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
        fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight", facecolor="white")


def set_sums(wide: pd.DataFrame, genes: list[str]) -> np.ndarray:
    hit = [g for g in genes if g in wide.index]
    if len(hit) != len(genes):
        raise SystemExit(f"missing from matrix: {set(genes) - set(hit)}")
    return wide.loc[hit, SPECIES_ORDER].sum(axis=0).to_numpy(float)


def build_matrices() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    print("локусы куриного белка + те же пики, сумма E", flush=True)
    fib_s, epi_s, lib_rows = {}, {}, []
    for sp in SPECIES_ORDER:
        print(f"  {sp}", flush=True)
        loci = s75.all_chicken_loci(s75.ASM / sp / "reanno_named.gff")
        peaks = load_peaks_both(sp)
        n_fib = int(peaks["n_fib"].iloc[0]) if len(peaks) else 0
        n_epi = int(peaks["n_epi"].iloc[0]) if len(peaks) else 0
        assigned = assign_both(peaks, loci)
        fib_s[sp] = assigned.set_index("gene")["fib"] if len(assigned) else pd.Series(dtype=float)
        epi_s[sp] = assigned.set_index("gene")["epi"] if len(assigned) else pd.Series(dtype=float)
        lib_rows.append(dict(species=sp, n_fib_libs=n_fib, n_epi_libs=n_epi,
                             n_genes=len(assigned)))
        print(f"    библиотек C={n_fib} E={n_epi}  генов {len(assigned)}", flush=True)
    Mf = pd.DataFrame(fib_s).reindex(columns=SPECIES_ORDER)
    Me = pd.DataFrame(epi_s).reindex(index=Mf.index, columns=SPECIES_ORDER)
    return Mf, Me, pd.DataFrame(lib_rows)


def five_grid(panels, ylabel: str, xlabel: str, stem: str, s40) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(12.4, 7.2))
    axes = axes.ravel()
    for i, (title, names, x, y, fit) in enumerate(panels):
        scatter(axes[i], names, x, y, xlabel, ylabel, s40)
        ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
        axes[i].set_title(f"{title}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=10, pad=7)
    axes[5].set_visible(False)
    fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.06, right=0.99, top=0.94, bottom=0.08)
    save(fig, stem, s40)
    plt.close(fig)


def main() -> None:
    s40 = load40()
    OUT.mkdir(parents=True, exist_ok=True)
    names = [sp.lower() for sp in SPECIES_ORDER]
    Mf, Me, libs = build_matrices()
    old = pd.read_csv(M_FIB_PATH, sep="\t", index_col=0)[SPECIES_ORDER]
    story_genes = []
    for stem, _, _ in SETS:
        story_genes.extend(genes_of(stem))
    story_genes = sorted(set(story_genes))
    rels = []
    for g in story_genes:
        a = Mf.loc[g].to_numpy(float)
        b = old.loc[g].to_numpy(float)
        rels.append(np.nanmax(np.abs(a - b) / np.clip(b, 1e-6, None)))
    max_rel = float(np.max(rels))
    print(f"сверка M_fib с скриптом 75 по {len(story_genes)} генам истории: max отн. {max_rel:.4g}")
    if max_rel > 1e-6:
        raise SystemExit("fibroblast assignment drifted; stop")
    Mf.to_csv(OUT / "M_fib_nearest2kb_check.tsv", sep="\t")
    Me.to_csv(OUT / "M_epi_nearest2kb.tsv", sep="\t")
    libs.to_csv(OUT / "library_counts.tsv", sep="\t", index=False)

    X = {key: np.array([fn(TRAITS[sp]) for sp in SPECIES_ORDER], float)
         for key, _, fn in PREDICTORS}
    xlab = {key: lab for key, lab, _ in PREDICTORS}

    pred_r = []
    keys = [k for k, _, _ in PREDICTORS]
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            rr = stats.linregress(X[a], X[b])
            pred_r.append(dict(a=a, b=b, r=float(rr.rvalue), p=float(rr.pvalue)))
            print(f"коллинеарность {a} vs {b}: r={rr.rvalue:+.3f} p={rr.pvalue:.3f}")
    pd.DataFrame(pred_r).to_csv(OUT / "predictor_collinearity.tsv", sep="\t", index=False)

    recs, canvas = [], []
    tissues = [("fib", Mf, "log10_fib"), ("epi", Me, "log10_epi")]
    for stem, title, figstem in SETS:
        genes = genes_of(stem)
        for tiss, wide, ylab in tissues:
            sums = set_sums(wide, genes)
            n_zero = int((wide.loc[genes] == 0).to_numpy().sum())
            n_na = int(wide.loc[genes].isna().to_numpy().sum())
            for key, xlab_i, _ in PREDICTORS:
                y = np.log10(sums)
                if not np.isfinite(y).all():
                    raise SystemExit(f"non-finite log10 sum: {stem} {tiss} {key} {sums}")
                fit = stats.linregress(X[key], y)
                recs.append(dict(
                    set=stem, title=title, n=len(genes), tissue=tiss,
                    predictor=key, xlabel=xlab_i, ylabel=ylab,
                    r=float(fit.rvalue), p=float(fit.pvalue),
                    slope=float(fit.slope), intercept=float(fit.intercept),
                    n_zero_cells=n_zero, n_na_cells=n_na,
                    **{f"sum_{sp.lower()}": float(sums[i]) for i, sp in enumerate(SPECIES_ORDER)},
                    **{f"y_{sp.lower()}": float(y[i]) for i, sp in enumerate(SPECIES_ORDER)},
                ))
                canvas.append(dict(
                    set=stem, title=title, tissue=tiss, predictor=key,
                    xlabel=xlab_i, ylabel=ylab,
                    r=float(fit.rvalue), p=float(fit.pvalue),
                    slope=float(fit.slope), intercept=float(fit.intercept),
                    x=X[key].tolist(), y=y.tolist(), names=names,
                ))
            print(f"{stem:16s} {tiss} zeros={n_zero} NA={n_na}", flush=True)

    fits = pd.DataFrame(recs)
    fits.to_csv(OUT / "story_set_trait_fits.tsv", sep="\t", index=False)
    (OUT / "canvas_panels.json").write_text(json.dumps(canvas, indent=2))

    # --- figures: epi vs mass (the locked recipe) ---
    epi_mass = []
    for stem, title, figstem in SETS:
        genes = genes_of(stem)
        sums = set_sums(Me, genes)
        y = np.log10(sums)
        x = X["mass"]
        fig, ax = plt.subplots(figsize=(4.6, 3.8))
        fit = scatter(ax, names, x, y, xlab["mass"], "log10_epi", s40)
        fig.subplots_adjust(left=0.16, right=0.96, top=0.88, bottom=0.16)
        save(fig, f"{figstem}_epi", s40)
        plt.close(fig)
        epi_mass.append((title, names, x, y, fit))
        print(f"{stem:16s} epi vs mass  r={fit['r']:+.3f} p={fit['p']:.4f}  "
              f"sums={[round(s, 1) for s in sums]}")
    five_grid(epi_mass, "log10_epi", xlab["mass"], "fig_story_five_epi", s40)

    # two-up fib | epi vs mass, one row per set
    fig, axes = plt.subplots(5, 2, figsize=(8.4, 16.2))
    for i, (stem, title, figstem) in enumerate(SETS):
        genes = genes_of(stem)
        for j, (tiss, wide, ylab) in enumerate(tissues):
            y = np.log10(set_sums(wide, genes))
            fit = scatter(axes[i, j], names, X["mass"], y, xlab["mass"], ylab, s40)
            ptxt = f"{fit['p']:.3f}" if fit["p"] >= 0.001 else f"{fit['p']:.3g}"
            axes[i, j].set_title(
                f"{title} {tiss}   cor. {fit['r']:.3f}; p = {ptxt}", fontsize=9, pad=6)
    fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.08, right=0.99, top=0.98, bottom=0.03)
    save(fig, "fig_story_five_fib_epi", s40)
    plt.close(fig)

    # three non-mass traits, five panels each tissue
    for key in ("incubation", "longevity", "genome"):
        for tiss, wide, ylab in tissues:
            panels = []
            for stem, title, _ in SETS:
                genes = genes_of(stem)
                y = np.log10(set_sums(wide, genes))
                fit = dict(r=float(stats.linregress(X[key], y).rvalue),
                           p=float(stats.linregress(X[key], y).pvalue))
                panels.append((title, names, X[key], y, fit))
            five_grid(panels, ylab, xlab[key], f"fig_story_five_{tiss}_{key}", s40)

    # compact r heatmap: 5 sets × 4 predictors, two tissues
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 3.6))
    pred_labs = [lab for _, lab, _ in PREDICTORS]
    set_labs = [t for _, t, _ in SETS]
    for ax, tiss, ttl in zip(axes, ("fib", "epi"), ("фибробласты", "эпибласт")):
        mat = np.array([
            [float(fits[(fits.set == stem) & (fits.tissue == tiss) & (fits.predictor == key)].r.iloc[0])
             for key, _, _ in PREDICTORS]
            for stem, _, _ in SETS
        ])
        im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(4), pred_labs, rotation=25, ha="right", fontsize=8)
        ax.set_yticks(range(5), set_labs, fontsize=9)
        ax.set_title(ttl, fontsize=11)
        for i in range(5):
            for j in range(4):
                ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=8,
                        color="#111" if abs(mat[i, j]) < 0.55 else "white")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.subplots_adjust(wspace=0.35, left=0.12, right=0.98, top=0.88, bottom=0.22)
    save(fig, "fig_story_r_heatmap_traits", s40)
    plt.close(fig)

    print(f"wrote {OUT}")
    show = fits[fits.predictor == "mass"][["set", "tissue", "r", "p"]]
    print(show.to_string(index=False))


if __name__ == "__main__":
    main()
