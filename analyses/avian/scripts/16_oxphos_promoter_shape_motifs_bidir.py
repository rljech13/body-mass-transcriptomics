#!/usr/bin/env python3
"""OXPHOS promoter shape: width, GC/motifs (chicken), bidirectionality (n=7).

Uses Dropbox fib+epi peak tables + galGal6 for chicken sequence windows.
"""
from __future__ import annotations

import json
import re
import subprocess
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
EXP = Path("/path/to/storage/kleiber_birds_gemcat/dropbox_20260808/unpacked/exp_tab")
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")
GENOME = Path("/path/to/user-home/chicken/CAGE_atlases/qc_reviewer/bpnet/galGal6.main.fa")
ALIAS = Path("/path/to/storage/kleiber_birds_gemcat/ref/galGal6.chromAlias.txt")
RES = ROOT / "results" / "oxphos_promoter_shape"
FIG = ROOT / "figures"

SPECIES = {
    "Zebra": {"file": "Zebra_fib_epi_exp.txt", "de": "Zebrafinch_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 11.93},
    "Quail": {"file": "Quail_fib_epi_exp.txt", "de": "Quail_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 96.6},
    "Chicken": {"file": "Chicken_fib_epi_exp.txt", "de": "Chicken_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 912.5},
    "Duck": {"file": "Duck_fib_epi_exp.txt", "de": "Duck_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 1121.0},
    "Peacock": {"file": "Peacock_fib_epi_exp.txt", "de": "Peacock_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 4093.75},
    "Turkey": {"file": "Turkey_fib_epi_exp.txt", "de": "Turkey_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 5811.0},
    "Emu": {"file": "Emu_fib_epi_exp.txt", "de": "Emu_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 36200.0},
}

PREFERRED_ALIAS = {
    "ATP5A1": "ATP5F1A", "ATP5B": "ATP5F1B", "ATP5C1": "ATP5F1C", "ATP5D": "ATP5F1D",
    "ATP5E": "ATP5F1E", "ATP5F1": "ATP5PB", "ATP5G1": "ATP5MC1", "ATP5G3": "ATP5MC3",
    "ATP5H": "ATP5PD", "ATP5J": "ATP5PF", "ATP5J2": "ATP5MF", "ATP5L": "ATP5MG",
    "ATP5O": "ATP5PO", "COX4": "COX4I1",
}
MTDNA = {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8"}


def load_ox() -> set[str]:
    return {ln.strip().upper() for ln in (ROOT / "data/genesets/OXPHOS.txt").read_text().splitlines() if ln.strip()}


def is_aug(s: str) -> bool:
    u = str(s).upper().strip()
    return bool(re.match(r"^G\d+(\.T\d+)?$", u)) or u in {"NAN", "NA", "", "NONE", "."}


def norm_sym(s: str, ox: set[str]) -> str:
    s = str(s).upper().strip()
    if not s or is_aug(s):
        return ""
    s2 = PREFERRED_ALIAS.get(s, s)
    return s2 if s2 in ox else (s if s in ox else "")


def annotate(df: pd.DataFrame, de_path: Path, ox: set[str]) -> pd.DataFrame:
    de = pd.read_excel(de_path, usecols=["id", "genes", "geneId"])
    de["sym"] = de["genes"].astype(str).where(~de["genes"].astype(str).str.upper().isin(["NAN", "NA"]), de["geneId"])
    de["sym"] = de["sym"].map(lambda x: norm_sym(x, ox))
    id2 = de[de.sym != ""].drop_duplicates("id").set_index("id")["sym"]
    out = df.copy()
    gid = out["geneId"].astype(str)
    cg = out["chicken_genes"].astype(str) if "chicken_genes" in out.columns else pd.Series([""] * len(out))
    sym = []
    for g, c in zip(gid, cg):
        pick = c if str(c) and not is_aug(c) and str(c).upper() not in {"NAN", "NA"} else (g if not is_aug(g) else "")
        sym.append(norm_sym(pick, ox))
    out["symbol"] = sym
    fill = out["symbol"].eq("")
    out.loc[fill, "symbol"] = out.loc[fill, "id"].map(id2).fillna("").map(lambda x: norm_sym(x, ox))
    return out


def lib_mean(df: pd.DataFrame, tissue: str) -> pd.Series:
    pref = "norm.C" if tissue == "fib" else "norm.E"
    cols = [c for c in df.columns if str(c).startswith(pref)]
    return df[cols].apply(pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)


def width_class(w: float) -> str:
    if w <= 100:
        return "sharp"
    if w <= 300:
        return "intermediate"
    return "broad"


def load_alias() -> dict[str, str]:
    out: dict[str, str] = {}
    if not ALIAS.exists():
        return out
    df = pd.read_csv(ALIAS, sep="\t", header=0)
    cols = [str(c).lstrip("# ").strip() for c in df.columns]
    df.columns = cols
    ucsc = "ucsc" if "ucsc" in cols else cols[0]
    for c in cols:
        if c == ucsc:
            continue
        for u, a in zip(df[ucsc], df[c]):
            if pd.notna(a) and pd.notna(u):
                out[str(a)] = str(u)
                out[str(u)] = str(u)
    return out


def fetch_seq(chrom: str, start: int, end: int) -> str:
    start = max(1, start)
    if end <= start:
        return ""
    try:
        out = subprocess.check_output(
            ["samtools", "faidx", str(GENOME), f"{chrom}:{start}-{end}"],
            text=True,
            stderr=subprocess.DEVNULL,
        )
        return "".join(out.splitlines()[1:]).upper()
    except Exception:
        return ""


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def seq_stats(seq: str) -> dict:
    if not seq:
        return {"gc": np.nan, "cpg": np.nan, "n": 0, "tata": 0, "inr": 0, "gcbox": 0, "caat": 0}
    s = seq.upper().replace("N", "")
    n = len(s)
    if n == 0:
        return {"gc": np.nan, "cpg": np.nan, "n": 0, "tata": 0, "inr": 0, "gcbox": 0, "caat": 0}
    gc = (s.count("G") + s.count("C")) / n
    cpg = len(re.findall("CG", s)) / max(n - 1, 1)
    return {
        "gc": gc,
        "cpg": cpg,
        "n": n,
        "tata": int(bool(re.search("TATA[AT]A", s))),
        "inr": int(bool(re.search("TCA[GTC]T", s))),  # loose Initiator
        "gcbox": int(bool(re.search("GGGCGG|CCGCCC", s))),
        "caat": int(bool(re.search("CCAAT|ATTGG", s))),
    }


def windows_for_tss(tss: int, strand: str) -> dict[str, tuple[int, int]]:
    """Relative to TSS, sense orientation after revcomp for - strand."""
    # genomic coords; for - strand we'll revcomp the fetched seq
    if strand == "+":
        return {
            "upstream_200_50": (tss - 200, tss - 50),
            "core_50_plus50": (tss - 50, tss + 50),
            "downstream_50_200": (tss + 50, tss + 200),
            "full_200_plus100": (tss - 200, tss + 100),
        }
    return {
        "upstream_200_50": (tss + 50, tss + 200),
        "core_50_plus50": (tss - 50, tss + 50),
        "downstream_50_200": (tss - 200, tss - 50),
        "full_200_plus100": (tss - 100, tss + 200),
    }


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    ox = load_ox()
    nuclear = ox - MTDNA
    alias = load_alias()

    peak_rows = []
    bidir_rows = []
    species_sum = []

    for sp, meta in SPECIES.items():
        print(f"=== {sp} ===", flush=True)
        df = pd.read_csv(EXP / meta["file"], sep="\t")
        df = annotate(df, RUSLAN / meta["de"], ox)
        df["width"] = pd.to_numeric(df["width"], errors="coerce")
        df["tss"] = pd.to_numeric(df["tss"], errors="coerce")
        df["strand"] = df["strand"].astype(str)
        df["is_promoter"] = df["class"].astype(str).str.lower().eq("promoter")
        for tissue in ("fib", "epi"):
            df[f"expr_{tissue}"] = lib_mean(df, tissue)

        # OXPHOS promoter peaks
        oxp = df[df["symbol"].isin(nuclear) & df["is_promoter"]].copy()
        bg = df[df["is_promoter"] & ~df["symbol"].isin(ox)].copy()  # non-OXPHOS promoters as background

        for tissue in ("fib", "epi"):
            ecol = f"expr_{tissue}"
            sub = oxp[oxp[ecol] >= 1.0].copy()
            if sub.empty:
                continue
            sub["width_class"] = sub["width"].map(width_class)
            for _, r in sub.iterrows():
                peak_rows.append(
                    {
                        "species": sp,
                        "mass_g": meta["mass_g"],
                        "tissue": tissue,
                        "gene": r["symbol"],
                        "peak_id": r["id"],
                        "width": r["width"],
                        "width_class": r["width_class"],
                        "expr": r[ecol],
                        "strand": r["strand"],
                        "tss": r["tss"],
                        "chrom": r["seqnames"],
                    }
                )

            # species summary vs background
            bg_t = bg[bg[ecol] >= 1.0]
            species_sum.append(
                {
                    "species": sp,
                    "mass_g": meta["mass_g"],
                    "tissue": tissue,
                    "n_oxphos_peaks": len(sub),
                    "n_oxphos_genes": sub["symbol"].nunique(),
                    "median_width_oxphos": float(sub["width"].median()),
                    "mean_width_oxphos": float(sub["width"].mean()),
                    "frac_sharp_oxphos": float((sub["width"] <= 100).mean()),
                    "frac_broad_oxphos": float((sub["width"] >= 300).mean()),
                    "median_width_bg": float(bg_t["width"].median()) if len(bg_t) else np.nan,
                    "frac_sharp_bg": float((bg_t["width"] <= 100).mean()) if len(bg_t) else np.nan,
                    "frac_broad_bg": float((bg_t["width"] >= 300).mean()) if len(bg_t) else np.nan,
                }
            )

            # bidirectionality: for each gene dominant peak, opposite-strand peak within 500 bp
            for gene, gdf in sub.groupby("symbol"):
                gdf = gdf.sort_values(ecol, ascending=False)
                top = gdf.iloc[0]
                tss = float(top["tss"])
                chrom = str(top["seqnames"])
                st = str(top["strand"])
                # opposite strand peaks same chrom near tss (any promoter with expr, including other genes / unassigned)
                near = df[
                    (df["seqnames"].astype(str) == chrom)
                    & (df["is_promoter"])
                    & (df[ecol] >= 1.0)
                    & (df["strand"].astype(str) != st)
                    & (df["tss"].sub(tss).abs() <= 500)
                ]
                # same-gene opposite (true bidirectional for this locus)
                same_gene_opp = gdf[(gdf["strand"].astype(str) != st) & (gdf["tss"].sub(tss).abs() <= 500)]
                # antisense peak not assigned to same gene
                antisense = near[~near["id"].isin(gdf["id"])]
                bidir_rows.append(
                    {
                        "species": sp,
                        "tissue": tissue,
                        "gene": gene,
                        "mass_g": meta["mass_g"],
                        "dom_width": float(top["width"]),
                        "dom_expr": float(top[ecol]),
                        "has_same_gene_opp_strand": int(len(same_gene_opp) > 0),
                        "has_antisense_peak_500bp": int(len(antisense) > 0),
                        "n_antisense_peaks_500bp": int(len(antisense)),
                        "max_antisense_expr": float(antisense[ecol].max()) if len(antisense) else 0.0,
                        "antisense_symbols": ",".join(
                            sorted({s for s in antisense["symbol"].astype(str) if s and s != "nan"})[:5]
                        ),
                    }
                )

    peaks = pd.DataFrame(peak_rows)
    peaks.to_csv(RES / "oxphos_peak_widths_n7.tsv", sep="\t", index=False)
    bidir = pd.DataFrame(bidir_rows)
    bidir.to_csv(RES / "oxphos_bidirectionality_n7.tsv", sep="\t", index=False)
    ssum = pd.DataFrame(species_sum)
    ssum.to_csv(RES / "species_width_summary.tsv", sep="\t", index=False)

    # OLS width / bidir vs mass (fib)
    fit_rows = []
    for tissue in ("fib", "epi"):
        sub = ssum[ssum.tissue == tissue]
        for col in ["median_width_oxphos", "frac_sharp_oxphos", "frac_broad_oxphos"]:
            x = np.log10(sub.mass_g)
            y = sub[col].astype(float)
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() >= 5:
                sl, ic, r, p, se = stats.linregress(x[m], y[m])
                fit_rows.append({"tissue": tissue, "metric": col, "slope": sl, "r": r, "p": p})
        bsub = bidir[bidir.tissue == tissue]
        bmean = bsub.groupby(["species", "mass_g"], as_index=False).agg(
            frac_antisense=("has_antisense_peak_500bp", "mean"),
            frac_same_gene_opp=("has_same_gene_opp_strand", "mean"),
        )
        for col in ["frac_antisense", "frac_same_gene_opp"]:
            x = np.log10(bmean.mass_g)
            y = bmean[col].astype(float)
            if len(y) >= 5:
                sl, ic, r, p, se = stats.linregress(x, y)
                fit_rows.append({"tissue": tissue, "metric": col, "slope": sl, "r": r, "p": p})
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(RES / "ols_shape_vs_mass.tsv", sep="\t", index=False)

    # ---------- Chicken sequence windows ----------
    seq_rows = []
    if GENOME.exists():
        ch = peaks[(peaks.species == "Chicken")].copy()
        # dominant peak per gene×tissue
        for tissue in ("fib", "epi"):
            sub = ch[ch.tissue == tissue].sort_values("expr", ascending=False).groupby("gene").head(1)
            # background sample: random non-oxphos promoters from chicken table
            cdf = pd.read_csv(EXP / "Chicken_fib_epi_exp.txt", sep="\t")
            cdf = annotate(cdf, RUSLAN / SPECIES["Chicken"]["de"], ox)
            cdf["expr"] = lib_mean(cdf, tissue)
            cdf = cdf[cdf["class"].astype(str).str.lower().eq("promoter") & (cdf["expr"] >= 1) & ~cdf["symbol"].isin(ox)]
            bg_sample = cdf.sample(n=min(200, len(cdf)), random_state=0) if len(cdf) else cdf

            def process(row, group: str):
                chrom_raw = str(row.get("chrom", row.get("seqnames")))
                chrom = alias.get(chrom_raw, chrom_raw if str(chrom_raw).startswith("chr") else None)
                if not chrom:
                    return
                tss = int(row["tss"])
                strand = str(row["strand"])
                for wname, (a, b) in windows_for_tss(tss, strand).items():
                    lo, hi = sorted((a, b))
                    seq = fetch_seq(chrom, lo, hi)
                    if strand == "-" and seq:
                        seq = revcomp(seq)
                    st = seq_stats(seq)
                    seq_rows.append(
                        {
                            "tissue": tissue,
                            "group": group,
                            "gene": row.get("gene", row.get("symbol", "")),
                            "window": wname,
                            "chrom": chrom,
                            "tss": tss,
                            "strand": strand,
                            "width": float(row.get("width", np.nan)),
                            **st,
                        }
                    )

            for _, r in sub.iterrows():
                process(r, "OXPHOS")
            for _, r in bg_sample.iterrows():
                process(
                    {
                        "chrom": r["seqnames"],
                        "tss": r["tss"],
                        "strand": r["strand"],
                        "gene": r.get("symbol", ""),
                        "width": r["width"],
                    },
                    "background",
                )

    seq = pd.DataFrame(seq_rows)
    if len(seq):
        seq.to_csv(RES / "chicken_sequence_windows_oxphos_vs_bg.tsv", sep="\t", index=False)
        # summarize
        summ = (
            seq.groupby(["tissue", "group", "window"])
            .agg(
                n=("gc", "count"),
                mean_gc=("gc", "mean"),
                mean_cpg=("cpg", "mean"),
                frac_tata=("tata", "mean"),
                frac_inr=("inr", "mean"),
                frac_gcbox=("gcbox", "mean"),
                frac_caat=("caat", "mean"),
            )
            .reset_index()
        )
        summ.to_csv(RES / "chicken_sequence_window_summary.tsv", sep="\t", index=False)

        # MW tests OXPHOS vs bg for core window fib
        tests = []
        for tissue in ("fib", "epi"):
            for window in seq.window.unique():
                a = seq[(seq.tissue == tissue) & (seq.group == "OXPHOS") & (seq.window == window)]["gc"].dropna()
                b = seq[(seq.tissue == tissue) & (seq.group == "background") & (seq.window == window)]["gc"].dropna()
                if len(a) >= 10 and len(b) >= 10:
                    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                    tests.append(
                        {
                            "tissue": tissue,
                            "window": window,
                            "metric": "gc",
                            "oxphos_median": float(a.median()),
                            "bg_median": float(b.median()),
                            "p_mw": float(p),
                        }
                    )
                for metric in ("cpg",):
                    a = seq[(seq.tissue == tissue) & (seq.group == "OXPHOS") & (seq.window == window)][metric].dropna()
                    b = seq[(seq.tissue == tissue) & (seq.group == "background") & (seq.window == window)][metric].dropna()
                    if len(a) >= 10 and len(b) >= 10:
                        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                        tests.append(
                            {
                                "tissue": tissue,
                                "window": window,
                                "metric": metric,
                                "oxphos_median": float(a.median()),
                                "bg_median": float(b.median()),
                                "p_mw": float(p),
                            }
                        )
        pd.DataFrame(tests).to_csv(RES / "chicken_oxphos_vs_bg_mw_tests.tsv", sep="\t", index=False)

    # ---------- figures ----------
    # 1 width distributions
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for ax, tissue in zip(axes, ("fib", "epi")):
        sub = peaks[peaks.tissue == tissue]
        if sub.empty:
            continue
        parts = [sub.loc[sub.species == s, "width"].dropna().values for s in SPECIES]
        ax.boxplot(parts, tick_labels=list(SPECIES.keys()), showfliers=False)
        ax.axhline(100, color="#2a9d8f", ls=":", lw=1, label="sharp≤100")
        ax.axhline(300, color="#e76f51", ls=":", lw=1, label="broad≥300")
        ax.set_ylabel("peak width (bp)")
        ax.set_title(f"OXPHOS promoter widths — {tissue}")
        ax.tick_params(axis="x", rotation=35)
        ax.legend(fontsize=7)
    fig.suptitle("Are OXPHOS CAGE peaks sharp or broad?", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_peak_widths_n7.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_peak_widths_n7.png", dpi=160)
    plt.close(fig)

    # 2 frac sharp / bidir vs mass
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), constrained_layout=True)
    fib_s = ssum[ssum.tissue == "fib"].sort_values("mass_g")
    fib_b = (
        bidir[bidir.tissue == "fib"]
        .groupby(["species", "mass_g"], as_index=False)
        .agg(frac_antisense=("has_antisense_peak_500bp", "mean"))
        .sort_values("mass_g")
    )
    for ax, y, ylab in zip(
        axes,
        [fib_s.frac_sharp_oxphos, fib_s.frac_broad_oxphos, fib_b.frac_antisense],
        ["frac sharp (≤100bp)", "frac broad (≥300bp)", "frac genes with antisense ≤500bp"],
    ):
        x = np.log10(fib_s.mass_g if len(y) == len(fib_s) else fib_b.mass_g)
        yy = y.values if hasattr(y, "values") else y
        ax.scatter(x, yy, c="#1d3557", s=50)
        labs = fib_s.species if len(y) == len(fib_s) else fib_b.species
        for xi, yi, s in zip(x, yy, labs):
            ax.text(xi, yi, s[:3], fontsize=7)
        if len(x) >= 5:
            sl, ic, r, p, se = stats.linregress(x, yy)
            xx = np.linspace(x.min(), x.max(), 50)
            ax.plot(xx, sl * xx + ic, color="#e76f51", lw=1)
            ax.set_title(f"p={p:.3f} r={r:.2f}")
        ax.set_xlabel("log10 mass")
        ax.set_ylabel(ylab)
    fig.suptitle("OXPHOS promoter shape / bidirectionality vs mass (fib)", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_shape_bidir_vs_mass.png", dpi=160)
    fig.savefig(RES / "fig_oxphos_shape_bidir_vs_mass.png", dpi=160)
    plt.close(fig)

    # 3 chicken GC by window
    if len(seq):
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
        for ax, tissue in zip(axes, ("fib", "epi")):
            sub = seq[(seq.tissue == tissue) & (seq.window == "core_50_plus50")]
            data = [sub.loc[sub.group == g, "gc"].dropna() for g in ("OXPHOS", "background")]
            ax.boxplot(data, tick_labels=["OXPHOS", "background"], showfliers=False)
            ax.set_ylabel("GC fraction")
            ax.set_title(f"Chicken core ±50 bp — {tissue}")
        fig.suptitle("GC at OXPHOS TSS core vs background promoters", fontsize=11)
        fig.savefig(FIG / "fig_oxphos_chicken_gc_core.png", dpi=160)
        fig.savefig(RES / "fig_oxphos_chicken_gc_core.png", dpi=160)
        plt.close(fig)

        # motif rates
        fig, ax = plt.subplots(figsize=(7.5, 4.0), constrained_layout=True)
        core = seq[(seq.tissue == "fib") & (seq.window == "full_200_plus100")]
        rows = []
        for g in ("OXPHOS", "background"):
            s = core[core.group == g]
            rows.append(
                {
                    "group": g,
                    "TATA": s.tata.mean(),
                    "Inr": s.inr.mean(),
                    "GC-box": s.gcbox.mean(),
                    "CAAT": s.caat.mean(),
                }
            )
        mdf = pd.DataFrame(rows).set_index("group")
        mdf.T.plot(kind="bar", ax=ax, color=["#e76f51", "#1d3557"])
        ax.set_ylabel("fraction of promoters with motif")
        ax.set_title("Chicken fib: crude motifs in −200/+100 around dominant OXPHOS TSS")
        ax.tick_params(axis="x", rotation=0)
        fig.savefig(FIG / "fig_oxphos_chicken_motifs.png", dpi=160)
        fig.savefig(RES / "fig_oxphos_chicken_motifs.png", dpi=160)
        plt.close(fig)

    # summary
    summary = {
        "width_defs": {"sharp": "≤100 bp", "intermediate": "100–300", "broad": "≥300 bp (cluster width from Dropbox)"},
        "bidir_def": "opposite-strand promoter peak within 500 bp of dominant OXPHOS TSS, expr≥1",
        "species_fib_width": ssum[ssum.tissue == "fib"].to_dict(orient="records"),
        "bidir_fib_frac": bidir[bidir.tissue == "fib"]
        .groupby("species")
        .agg(frac_antisense=("has_antisense_peak_500bp", "mean"), n=("gene", "count"))
        .reset_index()
        .to_dict(orient="records"),
        "ols": fits.to_dict(orient="records"),
        "chicken_seq_note": "galGal6 windows; motifs are loose regex (not TF ChIP).",
    }
    if len(seq):
        summary["chicken_core_gc"] = (
            seq[(seq.window == "core_50_plus50")]
            .groupby(["tissue", "group"])["gc"]
            .median()
            .reset_index()
            .to_dict(orient="records")
        )
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps({k: summary[k] for k in ("ols", "bidir_fib_frac", "chicken_core_gc") if k in summary}, indent=2))
    print("Wrote", RES)


if __name__ == "__main__":
    main()
