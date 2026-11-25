#!/usr/bin/env python3
"""Honest promoter reanalysis (not a claim of 'simplification').

1) Threshold sensitivity for OXPHOS CAGE peak-share metrics vs mass (n=7)
2) Same-strand filter reminder + dominant-peak fraction
3) Chicken: nucleotide composition around OXPHOS TSS (galGal6)
4) Chicken: orthogonal annotations from ChickCAGE peak passport (if joinable)

Language: report what metrics do; do NOT claim promoter architecture biology
until orthogonal support exists.
"""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
DATA = ROOT / "data"
RES = ROOT / "results" / "promoter_honest_reanalysis"
FIG = ROOT / "figures"
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")
GENOME = Path("/path/to/user-home/chicken/CAGE_atlases/qc_reviewer/bpnet/galGal6.main.fa")
ALIAS = Path("/path/to/storage/kleiber_birds_gemcat/ref/galGal6.chromAlias.txt")
PASSPORT = Path("/path/to/user-home/chicken/CAGE_atlases/release/atlas_peak_passport.tsv.gz")

SPECIES = {
    "Chicken": {"cage": "GSE213253_Chicken_CAGE.xlsx", "de": "Chicken_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 912.5},
    "Duck": {"cage": "GSE213253_Duck_CAGE.xlsx", "de": "Duck_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 1121.0},
    "Quail": {"cage": "GSE213253_Quail_CAGE.xlsx", "de": "Quail_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 96.6},
    "ZebraFinch": {"cage": "GSE213253_ZebraFinch_CAGE.xlsx", "de": "Zebrafinch_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 11.93},
    "Turkey": {"cage": "GSE213253_Turkey_CAGE.xlsx", "de": "Turkey_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 5811.0},
    "Peacock": {"cage": "GSE213253_Peacock_CAGE.xlsx", "de": "Peacock_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 4093.75},
    "Emu": {"cage": "GSE213253_Emu_CAGE.xlsx", "de": "Emu_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 36200.0},
}


def load_oxphos() -> set[str]:
    return {ln.strip().upper() for ln in (DATA / "genesets/OXPHOS.txt").read_text().splitlines() if ln.strip()}


def peak_to_symbol(de_path: Path) -> pd.Series:
    de = pd.read_excel(de_path, usecols=["id", "genes"])
    de["genes"] = de["genes"].astype(str).str.upper()
    de = de[de["genes"].notna() & (de["genes"] != "NAN")]
    return de.drop_duplicates("id").set_index("id")["genes"]


def gene_peak_table(cage: pd.DataFrame, id2sym: pd.Series, genes: set[str]) -> pd.DataFrame:
    norms = [c for c in cage.columns if str(c).startswith("norm.")]
    df = cage.copy()
    df["peak_id"] = df["CAGE_cluster_id"].astype(str)
    df["symbol"] = df["peak_id"].map(id2sym)
    # chicken GEO often has geneId as symbol
    if df["symbol"].isna().all() and "geneId" in df.columns:
        df["symbol"] = df["geneId"].astype(str).str.upper()
    df = df[df["symbol"].isin(genes)].copy()
    if df.empty:
        return df
    df["_expr"] = df[norms].apply(pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
    df["strand"] = df["CAGE_cluster_strand"].astype(str)
    df["tss"] = pd.to_numeric(df["CAGE_TSS"], errors="coerce")
    df["chr"] = df["CAGE_cluster_chr"].astype(str)
    return df


def metrics_for_gene(gdf: pd.DataFrame, tpm_thr: float, same_strand: bool) -> dict:
    g = gdf.copy()
    g = g[g["_expr"] >= tpm_thr]
    if g.empty:
        return {}
    if same_strand and len(g):
        # major strand by total expr
        maj = g.groupby("strand")["_expr"].sum().idxmax()
        g = g[g["strand"] == maj]
    if g.empty:
        return {}
    vals = g["_expr"].astype(float).values
    p = vals / vals.sum()
    H = float(-(p * np.log2(p + 1e-12)).sum())
    return {
        "n_active": int(len(g)),
        "shannon": H,
        "dominant_frac": float(p.max()),
        "total_tpm": float(vals.sum()),
        "simplicity_S": float(1.0 - H / np.log2(max(len(g), 2))),
    }


def species_module_means(gene_metrics: pd.DataFrame) -> pd.Series:
    return gene_metrics[["n_active", "shannon", "dominant_frac", "total_tpm", "simplicity_S"]].mean()


def fit_mass(df: pd.DataFrame, col: str) -> dict:
    x = np.log10(df["mass_g"].values)
    y = df[col].values.astype(float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return {"n": int(m.sum()), "slope": np.nan, "r": np.nan, "p": np.nan}
    slope, icept, r, p, se = stats.linregress(x[m], y[m])
    return {"n": int(m.sum()), "slope": float(slope), "r": float(r), "p": float(p), "se": float(se)}


def fetch_seq(fa_path: Path, chrom: str, start: int, end: int) -> str:
    """Fetch subsequence using samtools faidx if available, else manual scan (slow)."""
    import subprocess

    # galGal6 fa often uses chr1; GEO may use NC_... — caller must normalize
    region = f"{chrom}:{start}-{end}"
    try:
        out = subprocess.check_output(["samtools", "faidx", str(fa_path), region], text=True, stderr=subprocess.DEVNULL)
        seq = "".join(out.splitlines()[1:]).upper()
        return seq
    except Exception:
        return ""


def load_alias() -> dict[str, str]:
    """Map RefSeq NC_* / Ensembl / GenBank → UCSC chr*."""
    out: dict[str, str] = {}
    if not ALIAS.exists():
        return out
    # file starts with '# ucsc\tensembl\tgenbank\trefseq'
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


def nc_to_chr(nc: str, alias: dict[str, str]) -> str | None:
    s = str(nc)
    if s.startswith("chr") and s in alias.values():
        return s
    return alias.get(s)


def composition(seq: str) -> dict:
    if not seq:
        return {"gc": np.nan, "cpg": np.nan, "n": 0, "tata_like": np.nan, "initiator_like": np.nan}
    s = seq.upper().replace("N", "")
    n = len(s)
    if n == 0:
        return {"gc": np.nan, "cpg": np.nan, "n": 0, "tata_like": np.nan, "initiator_like": np.nan}
    gc = (s.count("G") + s.count("C")) / n
    cpg = len(re.findall("CG", s)) / max(n - 1, 1)
    # crude motifs
    tata = 1.0 if re.search("TATA[AT]A", s) else 0.0
    inr = 1.0 if re.search("TCA[GTC]T", s) else 0.0  # loose Initiator-like
    return {"gc": gc, "cpg": cpg, "n": n, "tata_like": tata, "initiator_like": inr}


def main() -> None:
    RES.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    ox = load_oxphos()

    # ---------- 1) threshold sensitivity ----------
    thresholds = [0.0, 1.0, 5.0, 10.0, 20.0]
    sens_rows = []
    per_gene_rows = []
    chicken_peaks = None

    for sp, meta in SPECIES.items():
        id2sym = peak_to_symbol(RUSLAN / meta["de"])
        cage = pd.read_excel(DATA / meta["cage"])
        tab = gene_peak_table(cage, id2sym, ox)
        if tab.empty:
            continue
        if sp == "Chicken":
            chicken_peaks = tab.copy()
        for thr in thresholds:
            for same in (False, True):
                gmets = []
                for gene, gdf in tab.groupby("symbol"):
                    m = metrics_for_gene(gdf, thr, same_strand=same)
                    if not m:
                        continue
                    gmets.append({"gene": gene, **m})
                    if thr == 1.0 and same:
                        per_gene_rows.append({"species": sp, "mass_g": meta["mass_g"], "gene": gene, **m})
                if not gmets:
                    continue
                gm = pd.DataFrame(gmets)
                means = gm[["n_active", "shannon", "dominant_frac", "total_tpm", "simplicity_S"]].mean()
                sens_rows.append(
                    {
                        "species": sp,
                        "mass_g": meta["mass_g"],
                        "tpm_thr": thr,
                        "same_strand": same,
                        "n_genes": len(gm),
                        **means.to_dict(),
                    }
                )

    sens = pd.DataFrame(sens_rows)
    sens.to_csv(RES / "threshold_sensitivity_species_means.tsv", sep="\t", index=False)
    pd.DataFrame(per_gene_rows).to_csv(RES / "per_gene_metrics_thr1_samestrand.tsv", sep="\t", index=False)

    fit_rows = []
    for thr in thresholds:
        for same in (False, True):
            sub = sens[(sens.tpm_thr == thr) & (sens.same_strand == same)]
            if len(sub) < 5:
                continue
            for col in ["shannon", "simplicity_S", "dominant_frac", "total_tpm", "n_active"]:
                f = fit_mass(sub, col)
                fit_rows.append({"tpm_thr": thr, "same_strand": same, "metric": col, **f, "n_spp_rows": len(sub)})
    fits = pd.DataFrame(fit_rows)
    fits.to_csv(RES / "threshold_sensitivity_ols_vs_mass.tsv", sep="\t", index=False)

    # ---------- 2) figure: sensitivity of Shannon p-value ----------
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for ax, same, title in zip(axes, [False, True], ["all strands (old-style)", "same-strand only"]):
        sub = fits[(fits.same_strand == same) & (fits.metric == "shannon")].sort_values("tpm_thr")
        ax.plot(sub.tpm_thr, sub.p, "o-", color="#1d3557", label="p (Shannon~mass)")
        subS = fits[(fits.same_strand == same) & (fits.metric == "total_tpm")].sort_values("tpm_thr")
        ax.plot(subS.tpm_thr, subS.p, "s--", color="#e76f51", label="p (total TPM~mass)")
        ax.axhline(0.05, color="#999", ls=":", lw=1)
        ax.set_xlabel("TPM threshold for 'active' peak")
        ax.set_ylabel("OLS p-value vs log10 mass")
        ax.set_title(title)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=8)
    fig.suptitle("Honest check: does 'promoter Shannon~mass' survive thresholds?\n(also compare to total TPM — expression confound)", fontsize=10)
    fig.savefig(FIG / "fig_promoter_threshold_sensitivity.png", dpi=160)
    fig.savefig(RES / "fig_promoter_threshold_sensitivity.png", dpi=160)
    plt.close(fig)

    # ---------- 3) chicken nucleotide composition around dominant OXPHOS TSS ----------
    alias = load_alias()
    comp_rows = []
    if chicken_peaks is not None and GENOME.exists():
        import subprocess

        # ensure fai
        if not Path(str(GENOME) + ".fai").exists():
            subprocess.run(["samtools", "faidx", str(GENOME)], check=False)
        for gene, gdf in chicken_peaks.groupby("symbol"):
            gdf = gdf.sort_values("_expr", ascending=False)
            top = gdf.iloc[0]
            chrom = nc_to_chr(str(top["chr"]), alias)
            if not chrom:
                continue
            tss = int(top["tss"])
            strand = str(top["strand"])
            if strand == "+":
                a, b = tss - 100, tss + 50
            else:
                a, b = tss - 50, tss + 100
            a = max(1, a)
            seq = fetch_seq(GENOME, chrom, a, b)
            if strand == "-" and seq:
                tr = str.maketrans("ACGT", "TGCA")
                seq = seq.translate(tr)[::-1]
            c = composition(seq)
            comp_rows.append(
                {
                    "gene": gene,
                    "chr_ucsc": chrom,
                    "chr_raw": str(top["chr"]),
                    "tss": tss,
                    "strand": strand,
                    "peak_expr": float(top["_expr"]),
                    "n_peaks_gene": int(len(gdf)),
                    **c,
                }
            )
    comp = pd.DataFrame(comp_rows)
    if len(comp):
        comp.to_csv(RES / "chicken_oxphos_tss_sequence_composition.tsv", sep="\t", index=False)
        summary_comp = {
            "n_genes_with_seq": int(len(comp)),
            "mean_gc": float(comp.gc.mean()),
            "mean_cpg": float(comp.cpg.mean()),
            "frac_tata_like": float(comp.tata_like.mean()),
            "frac_initiator_like": float(comp.initiator_like.mean()),
            "note": "Crude motifs; chicken fib dominant OXPHOS TSS ± window. Not a cross-species test.",
        }
    else:
        summary_comp = {"n_genes_with_seq": 0, "note": "No sequences fetched (chrom name mismatch or no chicken peaks)."}

    # ---------- 4) orthogonal passport join (chicken) ----------
    ortho_rows = []
    ortho_summary = {"passport_join": "skipped"}
    if chicken_peaks is not None and PASSPORT.exists():
        head = pd.read_csv(PASSPORT, sep="\t", nrows=2)
        cols = head.columns.tolist()
        (RES / "passport_columns.json").write_text(json.dumps(cols, indent=2))
        want_keys = [
            "chrom", "chr", "start", "end", "strand", "atlas_class", "create_original_class",
            "tissue_class", "expression_tier", "atac", "h3k4me3", "h3k27ac", "rampage",
            "bpnet", "tssclass", "is_likely_noise", "call_confidence", "peak_id", "cre_id",
        ]
        keep = []
        for key in want_keys:
            keep += [c for c in cols if key == c.lower() or key in c.lower()]
        keep = list(dict.fromkeys(keep))[:50]
        usecols = [c for c in cols if c in keep or c.lower() in {"chrom", "chr", "start", "end", "strand"}]
        pp = pd.read_csv(PASSPORT, sep="\t", usecols=lambda c: c in set(usecols))
        chrcol = "chrom" if "chrom" in pp.columns else ("chr" if "chr" in pp.columns else None)
        if chrcol and "start" in pp.columns:
            chicken_peaks = chicken_peaks.copy()
            chicken_peaks["_chr"] = chicken_peaks["chr"].map(lambda x: nc_to_chr(str(x), alias))
            pp["_chr"] = pp[chrcol].astype(str)
            strand_col = "strand" if "strand" in pp.columns else None
            # index by chrom for speed
            for chrom, gpeaks in chicken_peaks.dropna(subset=["_chr"]).groupby("_chr"):
                subpp = pp[pp["_chr"] == chrom]
                if subpp.empty:
                    continue
                starts = pd.to_numeric(subpp["start"], errors="coerce")
                for _, r in gpeaks.iterrows():
                    mask = pd.Series(True, index=subpp.index)
                    if strand_col:
                        mask &= subpp[strand_col].astype(str) == str(r["strand"])
                    cand = subpp[mask]
                    if cand.empty:
                        continue
                    dist = (starts.loc[cand.index] - r["tss"]).abs()
                    jidx = dist.idxmin()
                    if dist.loc[jidx] > 50:
                        continue
                    j = cand.loc[jidx]
                    ortho_rows.append(
                        {
                            "gene": r["symbol"],
                            "tss": int(r["tss"]),
                            "chr_ucsc": chrom,
                            "peak_expr": float(r["_expr"]),
                            "dist_to_passport": float(dist.loc[jidx]),
                            **{k: j[k] for k in keep if k in j.index},
                        }
                    )
            if ortho_rows:
                od = pd.DataFrame(ortho_rows).drop_duplicates(subset=["gene", "tss"])
                od.to_csv(RES / "chicken_oxphos_peaks_passport_join.tsv", sep="\t", index=False)
                ortho_summary = {"passport_join": "coordinate_nearest_alias", "n_joined": int(len(od)), "columns_kept": keep[:25]}
            else:
                ortho_summary = {"passport_join": "no_hits_within_50bp", "n_oxphos_peaks": int(len(chicken_peaks))}

    # ---------- 5) partial correlation: Shannon ~ mass | TPM ----------
    partial = {}
    sub = sens[(sens.tpm_thr == 1.0) & (sens.same_strand == True)].copy()
    if len(sub) >= 5:
        # residualize shannon and logM on log(total_tpm)
        y = sub["shannon"].values.astype(float)
        x = np.log10(sub["mass_g"].values.astype(float))
        z = np.log10(sub["total_tpm"].values.astype(float) + 1e-6)
        def resid(a, b):
            slope, icept, *_ = stats.linregress(b, a)
            return a - (icept + slope * b)
        y_r = resid(y, z)
        x_r = resid(x, z)
        slope, icept, r, p, se = stats.linregress(x_r, y_r)
        partial = {
            "method": "OLS residualize both Shannon and log10(mass) on log10(total_tpm)",
            "slope": float(slope),
            "r": float(r),
            "p": float(p),
            "interpretation": "If p stays small after removing TPM, mass trend is not only expression level.",
        }
        (RES / "partial_shannon_mass_given_tpm.json").write_text(json.dumps(partial, indent=2))

    # ---------- summary ----------
    # highlight: at thr=1 same-strand, is shannon still sig? is TPM also?
    def pick(thr, same, metric):
        hit = fits[(fits.tpm_thr == thr) & (fits.same_strand == same) & (fits.metric == metric)]
        return hit.iloc[0].to_dict() if len(hit) else {}

    summary = {
        "language": "Do NOT claim promoter simplification; report threshold sensitivity and confounds.",
        "shannon_thr1_allstrand": pick(1.0, False, "shannon"),
        "shannon_thr1_samestrand": pick(1.0, True, "shannon"),
        "tpm_thr1_samestrand": pick(1.0, True, "total_tpm"),
        "shannon_thr10_samestrand": pick(10.0, True, "shannon"),
        "sequence_composition_chicken": summary_comp,
        "orthogonal": ortho_summary,
        "partial_shannon_mass_given_tpm": partial,
        "next_steps": [
            "Build NC_→chr map for galGal6 and re-join passport/ATAC/H3K4me3",
            "Browser QC of multi-peak OXPHOS genes (same-strand only)",
            "Partial correlation: Shannon ~ mass | total_tpm",
            "Await salmon all7 for denser OXPHOS gene set in expression space",
        ],
    }
    (RES / "SUMMARY.json").write_text(json.dumps(summary, indent=2))
    (RES / "README_honest.md").write_text(
        """# Promoter honest reanalysis

This folder does **not** claim that OXPHOS promoters 'simplify' with body mass.

Shannon / S here are statistics of **how CAGE TPM is split among already-called peaks**.
They are entangled with expression level (total TPM).

See:
- `threshold_sensitivity_ols_vs_mass.tsv`
- `fig_promoter_threshold_sensitivity.png`
- `chicken_oxphos_tss_sequence_composition.tsv` (if sequences fetched)
- `SUMMARY.json`

Human-readable project verdict:
`corpus_kleiber_2026-08-06/03_документный_слой/Честный_вердикт_что_сделали_и_чего_не_видим_v1.0.md`
"""
    )
    print(json.dumps(summary, indent=2))
    print("Wrote", RES)


if __name__ == "__main__":
    main()
