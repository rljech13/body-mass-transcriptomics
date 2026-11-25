#!/usr/bin/env python3
"""Shared helpers for OXPHOS species-thread pipeline."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
DATA = ROOT / "data"
EXP = Path("/path/to/storage/kleiber_birds_gemcat/dropbox_20260808/unpacked/exp_tab")
RUSLAN = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes")
RES = ROOT / "results" / "oxphos_species_thread"
FIG = ROOT / "figures"
GENOME_ROOT = Path("/path/to/storage/kleiber_birds_gemcat/ref/genomes_with_mt")
FIMO = Path("/path/to/user-home/miniconda3/envs/meme/bin/fimo")
SALMON = ROOT / "results" / "salmon_all7_chicken" / "gene_tpm_by_species.tsv"

MTDNA = {
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
    "CYTB", "COX1", "COX2", "COX3", "ATP6", "ATP8",
}

PREFERRED_ALIAS = {
    "ATP5A1": "ATP5F1A", "ATP5B": "ATP5F1B", "ATP5C1": "ATP5F1C",
    "ATP5D": "ATP5F1D", "ATP5E": "ATP5F1E", "ATP5F1": "ATP5PB",
    "ATP5G1": "ATP5MC1", "ATP5G3": "ATP5MC3", "ATP5H": "ATP5PD",
    "ATP5J": "ATP5PF", "ATP5J2": "ATP5MF", "ATP5L": "ATP5MG",
    "ATP5O": "ATP5PO", "COX4": "COX4I1", "PGC1A": "PPARGC1A",
    "NRF2": "NFE2L2", "ERRA": "ESRRA",
}

SPECIES = {
    "Zebra": {"file": "Zebra_fib_epi_exp.txt", "de": "Zebrafinch_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 11.93, "salmon": "ZebraFinch", "geo_rna": "GSE213253_Zebrafinch_RNAseq.xlsx"},
    "Quail": {"file": "Quail_fib_epi_exp.txt", "de": "Quail_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 96.6, "salmon": "Quail", "geo_rna": "GSE213253_Quail_RNAseq.xlsx"},
    "Chicken": {"file": "Chicken_fib_epi_exp.txt", "de": "Chicken_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 912.5, "salmon": "Chicken", "geo_rna": "GSE213253_Chicken_RNAseq.xlsx"},
    "Duck": {"file": "Duck_fib_epi_exp.txt", "de": "Duck_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 1121.0, "salmon": "Duck", "geo_rna": "GSE213253_Duck_RNAseq.xlsx"},
    "Peacock": {"file": "Peacock_fib_epi_exp.txt", "de": "Peacock_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 4093.75, "salmon": "Peacock", "geo_rna": "GSE213253_Peacock_RNAseq.xlsx"},
    "Turkey": {"file": "Turkey_fib_epi_exp.txt", "de": "Turkey_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 5811.0, "salmon": "Turkey", "geo_rna": "GSE213253_Turkey_RNAseq.xlsx"},
    "Emu": {"file": "Emu_fib_epi_exp.txt", "de": "Emu_DE_epi(up)_fib(down)_promoters.xlsx", "mass_g": 36200.0, "salmon": "Emu", "geo_rna": "GSE213253_Emu_RNAseq.xlsx"},
}
SPECIES_ORDER = list(SPECIES.keys())

# High DMEM 1F1M sample codes that overlap RNA-seq subset (from GSE213253)
MATCHED_HIGH_CODES = {
    "Chicken": ["C10", "C12"],
    "Quail": ["C26", "C28"],
    "Duck": ["C32", "C34"],  # may appear as C32.D.F style in RNA; CAGE cols use Cxx
    "Peacock": ["C06", "C04"],
    "Turkey": ["C16", "C18"],
    "Emu": ["C22", "C24"],
    "Zebra": ["C39", "C40"],
}

DRIVERS = ["NRF1", "GABPA", "ESRRA", "YY1"]
DRIVERS_EXTRA = ["TFAM", "POLRMT", "SIRT1", "PPRC1"]

# galGal6 RefSeq → chr (common)
CHICKEN_NC_TO_CHR = {
    "NC_006088.5": "chr1", "NC_006089.5": "chr2", "NC_006090.5": "chr3",
    "NC_006091.5": "chr4", "NC_006092.5": "chr5", "NC_006093.5": "chr6",
    "NC_006094.5": "chr7", "NC_006095.5": "chr8", "NC_006096.5": "chr9",
    "NC_006097.5": "chr10", "NC_006098.5": "chr11", "NC_006099.5": "chr12",
    "NC_006100.5": "chr13", "NC_006101.5": "chr14", "NC_006102.5": "chr15",
    "NC_006103.5": "chr16", "NC_006104.5": "chr17", "NC_006105.5": "chr18",
    "NC_006106.5": "chr19", "NC_006107.5": "chr20", "NC_006108.5": "chr21",
    "NC_006109.5": "chr22", "NC_006110.5": "chr23", "NC_006111.5": "chr24",
    "NC_006112.5": "chr25", "NC_006113.5": "chr26", "NC_006114.5": "chr27",
    "NC_006115.5": "chr28", "NC_006119.5": "chr32", "NC_006126.5": "chr33",
    "NC_006127.5": "chrW", "NC_006088": "chr1",
    "NC_001323.1": "chrM",
}


def load_geneset(name: str) -> set[str]:
    path = DATA / "genesets" / f"{name}.txt"
    return {ln.strip().upper() for ln in path.read_text().splitlines() if ln.strip()}


def nuclear_oxphos() -> set[str]:
    return {g for g in load_geneset("OXPHOS") if g not in MTDNA}


def is_augustus_id(s: str) -> bool:
    u = str(s).upper().strip()
    return bool(re.match(r"^G\d+(\.T\d+)?$", u)) or u in {"NAN", "NA", "", "NONE", "."}


def normalize_symbol(sym: str) -> str:
    s = str(sym).upper().strip()
    if not s or s in {"NAN", "NA", "NONE", "."} or is_augustus_id(s):
        return ""
    return PREFERRED_ALIAS.get(s, s)


def peak_to_symbol_de(de_path: Path) -> pd.Series:
    de = pd.read_excel(de_path, usecols=["id", "genes", "geneId"])
    for col in ("genes", "geneId"):
        de[col] = de[col].astype(str).str.upper()
    de["sym"] = de["genes"].where(~de["genes"].isin(["NAN", "NA", ""]), de["geneId"])
    de["sym"] = de["sym"].map(normalize_symbol)
    de = de[de["sym"] != ""].drop_duplicates("id")
    return de.set_index("id")["sym"]


def annotate_table(df: pd.DataFrame, id2sym: pd.Series) -> pd.DataFrame:
    out = df.copy()
    gid = out["geneId"].astype(str)
    cg = out["chicken_genes"].astype(str) if "chicken_genes" in out.columns else pd.Series([""] * len(out))
    sym = []
    for g, c in zip(gid, cg):
        cu = str(c).upper().strip()
        gu = str(g).upper().strip()
        pick = ""
        if cu and not is_augustus_id(cu):
            pick = cu
        elif not is_augustus_id(gu):
            pick = gu
        sym.append(normalize_symbol(pick))
    out["symbol"] = sym
    fill = out["symbol"].eq("")
    out.loc[fill, "symbol"] = out.loc[fill, "id"].map(id2sym).fillna("").map(normalize_symbol)
    return out


def lib_cols(df: pd.DataFrame, tissue: str) -> list[str]:
    pref = "norm.C" if tissue == "fib" else "norm.E"
    return [c for c in df.columns if str(c).startswith(pref)]


def filter_matched_cols(cols: list[str], codes: list[str]) -> list[str]:
    keep = []
    for c in cols:
        for code in codes:
            if re.search(rf"(?:^|[._]){code}(?:[._]|$)", str(c), flags=re.I) or code in str(c):
                keep.append(c)
                break
    return keep


def gene_by_lib_matrix(
    df: pd.DataFrame,
    genes: set[str],
    tissue: str,
    thr: float = 1.0,
    cols: list[str] | None = None,
) -> pd.DataFrame:
    """Rows=genes, cols=libs; total same-strand promoter TPM per lib."""
    cols = cols or lib_cols(df, tissue)
    if not cols:
        return pd.DataFrame()
    sub = df[df["class"].astype(str).str.lower().eq("promoter") & df["symbol"].isin(genes)]
    if sub.empty:
        return pd.DataFrame(columns=cols)
    sub = sub.copy()
    mat = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    sub = sub.assign(**{c: mat[c] for c in cols})
    sub["_tot"] = mat.sum(axis=1)
    # major strand per gene by total TPM
    strand_tot = sub.groupby(["symbol", "strand"], observed=True)["_tot"].sum().reset_index()
    maj = strand_tot.sort_values("_tot", ascending=False).drop_duplicates("symbol")
    maj_map = maj.set_index("symbol")["strand"]
    sub = sub[sub["strand"].values == sub["symbol"].map(maj_map).values]
    # threshold per peak per lib, then sum
    for c in cols:
        sub[c] = sub[c].where(sub[c] >= thr, 0.0)
    out = sub.groupby("symbol", observed=True)[cols].sum()
    return out


def median_pairwise_spearman(mat: pd.DataFrame) -> float:
    """mat: genes x samples. Return median pairwise Spearman among genes with variance."""
    m = mat.loc[mat.std(axis=1) > 1e-12]
    if len(m) < 3:
        return float("nan")
    # rank along samples, then pearson = spearman
    ranks = m.rank(axis=1, method="average").to_numpy(dtype=float)
    # correlate genes (rows)
    c = np.corrcoef(ranks)
    if c.ndim < 2:
        return float("nan")
    iu = np.triu_indices(c.shape[0], k=1)
    vals = c[iu]
    vals = vals[np.isfinite(vals)]
    return float(np.median(vals)) if len(vals) else float("nan")


def ols_logmass(y: pd.Series, masses: pd.Series) -> dict:
    x = np.log10(masses.astype(float).values)
    yy = y.astype(float).values
    mask = np.isfinite(x) & np.isfinite(yy)
    if mask.sum() < 4:
        return {"n": int(mask.sum()), "slope": np.nan, "r": np.nan, "p": np.nan, "se": np.nan, "intercept": np.nan}
    sl, ic, r, p, se = stats.linregress(x[mask], yy[mask])
    return {"n": int(mask.sum()), "slope": float(sl), "r": float(r), "p": float(p), "se": float(se), "intercept": float(ic)}


def residual_after_covariate(y: pd.Series, cov: pd.Series) -> pd.Series:
    x = cov.astype(float).values
    yy = y.astype(float).values
    mask = np.isfinite(x) & np.isfinite(yy)
    out = pd.Series(np.nan, index=y.index)
    if mask.sum() < 3:
        return out
    sl, ic, *_ = stats.linregress(x[mask], yy[mask])
    pred = ic + sl * cov.astype(float)
    out = y.astype(float) - pred
    return out


def spearman_logmass(y: pd.Series, masses: pd.Series) -> dict:
    x = np.log10(masses.astype(float).values)
    yy = y.astype(float).values
    mask = np.isfinite(x) & np.isfinite(yy)
    if mask.sum() < 4:
        return {"n": int(mask.sum()), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[mask], yy[mask])
    return {"n": int(mask.sum()), "rho": float(rho), "p": float(p)}


def de_augustus_to_symbol(de_path: Path) -> dict[str, str]:
    """Map Augustus transcript/gene IDs (g9190.t1, g9190) → symbol via Ruslan DE."""
    de = pd.read_excel(de_path, usecols=["geneId", "transcriptId", "genes"])
    out: dict[str, str] = {}
    for _, r in de.iterrows():
        sym = normalize_symbol(r["genes"]) or normalize_symbol(r["geneId"])
        if not sym:
            continue
        tid = str(r["transcriptId"]).strip()
        gid = str(r["geneId"]).strip()
        keys = {tid, gid}
        if "." in tid:
            keys.add(tid.split(".", 1)[0])
        for key in keys:
            if key and key.lower() not in {"nan", "none", ""}:
                out[key.upper()] = sym
    return out


def load_geo_rnaseq(species: str) -> dict:
    """Official GSE213253 RNA-seq: symbol-indexed gene×lib norm matrix + library totals.

    Peacock/Turkey/Emu use Augustus IDs; mapped through DE transcriptId.
    GEO 'norm' totals are NOT a shared 1e6 TPM scale across species — use
    within-table fractions, not raw interspecific means, as the primary metric.
    """
    meta = SPECIES[species]
    df = pd.read_excel(DATA / meta["geo_rna"])
    norms = [c for c in df.columns if str(c).startswith("norm.")]
    raw_id = df["gene"].astype(str)
    mat = df[norms].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib_total = mat.sum(axis=0)
    if raw_id.map(is_augustus_id).mean() < 0.5:
        sym = raw_id.map(normalize_symbol)
        mapped_how = "native_symbol"
    else:
        id2 = de_augustus_to_symbol(RUSLAN / meta["de"])
        keys = raw_id.str.upper()
        keys2 = keys.str.replace(r"\.T\d+$", "", regex=True)
        mapped = keys.map(id2)
        miss = mapped.isna() | (mapped == "")
        mapped = mapped.where(~miss, keys2.map(id2))
        sym = mapped.fillna("")
        mapped_how = "de_transcriptId"
    keep = sym != ""
    tmp = mat.loc[keep].copy()
    tmp["_sym"] = sym[keep].values
    gene_mat = tmp.groupby("_sym", observed=True)[norms].mean()
    return {
        "gene_mat": gene_mat,
        "lib_total": lib_total,
        "mapped_how": mapped_how,
        "n_raw": int(len(df)),
        "n_mapped_rows": int(keep.sum()),
        "n_symbols": int(gene_mat.shape[0]),
        "norm_cols": norms,
    }
