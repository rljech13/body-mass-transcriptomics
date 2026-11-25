#!/usr/bin/env python3
"""FANTOM5 cultured mesenchymal CAGE: OXPHOS mean TPM vs mass.

Closest public analog to Ruslan's fibroblast CAGE that actually adds species.
Primary cell type: bone-marrow MSC. Sensitivity: unmatched aortic SMC
(official tables in all five species).

Do not pool with GSE213253 embryonic fibroblasts. Chicken MSC is missing
from the official peak-TPM matrix; it is recovered from the datahub TPM
bigWigs of the single donor library CNhs11292.
"""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd
import pyBigWig

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    EXP,
    PREFERRED_ALIAS,
    RES,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    nuclear_oxphos,
    normalize_symbol,
    ols_logmass,
    peak_to_symbol_de,
    spearman_logmass,
)
from oxphos_thread_lib import RUSLAN  # noqa: E402

FANTOM = Path("/path/to/storage/kleiber_birds_gemcat/fantom5_msc")
OUT = RES / "fantom_msc_oxphos"

# Adult masses (g). Dog: FANTOM donors are beagle (~10 kg), not AnAge 40 kg.
MASS = {
    "mouse": 20.5,
    "rat": 300.0,
    "chicken": 912.5,
    "dog": 10000.0,
    "human": 62035.0,
}
MASS_DOG_ANAGE = 40000.0

HUMAN_MSC = ["CNhs11344", "CNhs12100", "CNhs12126", "CNhs11316"]
MOUSE_MSC = ["CNhs12628", "CNhs12633", "CNhs12616"]
HUMAN_AORTIC = ["CNhs11085", "CNhs11305", "CNhs11309"]
MOUSE_AORTIC = ["CNhs11297"]
HUMAN_DERMAL = ["CNhs12499", "CNhs11379", "CNhs12028", "CNhs12052", "CNhs12055", "CNhs12059"]
MOUSE_MEF = ["CNhs12130"]

JUNK_NAME = re.compile(r"^(g\d+|aabr|ens[a-z]{3}g|loc\d+|chr\d)", re.I)


def oxphos_targets() -> set[str]:
    genes = nuclear_oxphos()
    extra = set()
    for g in list(genes):
        extra.add(PREFERRED_ALIAS.get(g, g))
    return {normalize_symbol(g) for g in genes | extra if normalize_symbol(g)}


def parse_osc_header(path: Path) -> tuple[list[str], int]:
    """Return (column names of first data row, n comment lines)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        n = 0
        for line in fh:
            if line.startswith("#"):
                n += 1
                continue
            cols = line.rstrip("\n").split("\t")
            return cols, n
    raise RuntimeError(f"no header in {path}")


def load_osc_tpm(path: Path) -> pd.DataFrame:
    header, skip = parse_osc_header(path)
    df = pd.read_csv(path, sep="\t", skiprows=skip + 1, header=None, names=header)
    df = df.rename(columns={header[0]: "peak_id"})
    return df


def classify_col(name: str) -> str | None:
    n = unquote(name).replace("%20", " ").lower()
    if "universal rna" in n or "whole body" in n or "extraembryonic" in n:
        return None
    if "leg bud" in n or "wing bud" in n:
        return None
    if "st2" in n or "response to" in n or "fraction" in n:
        return None
    if "mesenchymal stem" in n:
        return "msc"
    if "aortic" in n and "smooth muscle" in n and "differentiated" in n:
        return "smc_diff"
    if "aortic" in n and "smooth muscle" in n:
        return "smc"
    if "hepatocyte" in n:
        return "hep"
    return None


def extract_cnhs(name: str) -> str:
    m = re.search(r"(CNhs\d+)", name)
    return m.group(1) if m else ""


def map_rat_dog(ann: pd.DataFrame) -> pd.Series:
    name = ann["Ensembl_gene_name"].astype(str)
    ref = ann.get("RefSeq_transcript_id", pd.Series([""] * len(ann))).astype(str)
    pid = ann["CAGE_peaks_ID"].astype(str)
    out = []
    for n, r in zip(name, ref):
        s = normalize_symbol(n)
        if s and not JUNK_NAME.match(s) and not s.startswith("ENS"):
            out.append(s)
            continue
        out.append("")
    return pd.Series(out, index=pid, name="symbol")


def chicken_tx2sym(gtf: Path) -> dict[str, str]:
    tx = {}
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if "\ttranscript\t" not in line and '\tgene\t' not in line:
                continue
            m_tid = re.search(r'transcript_id "([^"]+)"', line)
            m_gn = re.search(r'gene_name "([^"]+)"', line)
            if not m_gn:
                continue
            sym = normalize_symbol(m_gn.group(1))
            if not sym:
                continue
            if m_tid:
                tid = m_tid.group(1)
                tx[tid] = sym
                tx[tid.split(".")[0]] = sym
    return tx


def map_chicken(ann: pd.DataFrame, tx2sym: dict[str, str]) -> pd.Series:
    pid = ann.iloc[:, 0].astype(str)
    raw = ann.iloc[:, 1].astype(str)
    out = []
    for cell in raw:
        hit = ""
        for tok in str(cell).split(","):
            tok = tok.strip()
            hit = tx2sym.get(tok) or tx2sym.get(tok.split(".")[0], "")
            if hit:
                break
        out.append(hit)
    return pd.Series(out, index=pid, name="symbol")


def score_peak_table(tpm: pd.DataFrame, id2sym: pd.Series, genes: set[str], cols: list[str]) -> dict[str, dict]:
    sub = tpm.set_index("peak_id")
    keep = [c for c in cols if c in sub.columns]
    if not keep:
        return {}
    num = sub[keep].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    num["symbol"] = id2sym.reindex(num.index).fillna("")
    num = num[num["symbol"].isin(genes)]
    rows = {}
    for c in keep:
        gmax = num.groupby("symbol")[c].max()
        gsum = num.groupby("symbol")[c].sum()
        present = gmax[gmax > 0]
        rows[c] = {
            "n_genes_any": int((gmax > 0).sum()),
            "mean_max": float(present.mean()) if len(present) else float("nan"),
            "mean_sum": float(gsum[gsum > 0].mean()) if (gsum > 0).any() else float("nan"),
            "cnhs": extract_cnhs(c),
        }
    return rows


def chicken_msc_from_bw(id2sym: pd.Series, genes: set[str], peak_bed: Path) -> dict:
    fwd = pyBigWig.open(str(FANTOM / "chicken_msc.fwd.bw"))
    rev = pyBigWig.open(str(FANTOM / "chicken_msc.rev.bw"))
    chroms = set(fwd.chroms()) | set(rev.chroms())
    recs = []
    with gzip.open(peak_bed, "rt") as fh:
        for line in fh:
            p = line.split("\t")
            chrom, start, end, name, strand = p[0], int(p[1]), int(p[2]), p[3], p[5]
            if chrom not in chroms or end <= start:
                continue
            bw = fwd if strand == "+" else rev
            try:
                val = bw.stats(chrom, start, end, type="sum")[0]
            except RuntimeError:
                val = 0.0
            recs.append((name, 0.0 if val is None else float(val)))
    fwd.close()
    rev.close()
    s = pd.Series(dict(recs), dtype=float)
    s.index = s.index.map(str)
    sym = id2sym.reindex(s.index).fillna("")
    df = pd.DataFrame({"tpm": s.values, "symbol": sym.values}, index=s.index)
    df = df[df["symbol"].isin(genes)]
    gmax = df.groupby("symbol")["tpm"].max()
    gsum = df.groupby("symbol")["tpm"].sum()
    present = gmax[gmax > 0]
    return {
        "n_genes_any": int((gmax > 0).sum()),
        "mean_max": float(present.mean()) if len(present) else float("nan"),
        "mean_sum": float(gsum[gsum > 0].mean()) if (gsum > 0).any() else float("nan"),
        "cnhs": "CNhs11292",
        "source": "datahub_tpm_bigwig",
    }


def gene_col_matches(colnames: list[str], ids: list[str]) -> list[str]:
    hit = []
    for c in colnames:
        for i in ids:
            if i in c:
                hit.append(c)
                break
    return hit


def fantom_gene_symbol(raw: str) -> str:
    """Parse FANTOM gene-level 00Annotation / short_description."""
    s = str(raw)
    # common: entrez:123,hgnc:ATP5F1A,ensg:...
    m = re.search(r"(?:hgnc:|HGNC:|symbol:)([A-Za-z0-9-]+)", s, re.I)
    if m:
        return normalize_symbol(m.group(1))
    # p1@SYMBOL or p2@Symbol
    m = re.search(r"p\d+@([A-Za-z0-9-]+)", s)
    if m:
        return normalize_symbol(m.group(1))
    # last pipe / comma token that looks like a symbol
    for tok in re.split(r"[,;| ]+", s):
        t = normalize_symbol(tok)
        if t and not t.startswith("ENS") and not t.isdigit() and len(t) >= 3:
            if not JUNK_NAME.match(t):
                return t
    return normalize_symbol(s)


def score_gene_level(path: Path, lib_ids: list[str], genes: set[str]) -> list[dict]:
    header, skip = parse_osc_header(path)
    cols = gene_col_matches(header, lib_ids)
    if not cols:
        return []
    usecols = [header[0]] + cols
    df = pd.read_csv(path, sep="\t", skiprows=skip + 1, header=None, names=header, usecols=usecols)
    ann = df[header[0]].map(fantom_gene_symbol)
    df = df.assign(symbol=ann)
    df = df[df["symbol"].isin(genes)]
    rows = []
    for c in cols:
        v = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
        by = pd.DataFrame({"symbol": df["symbol"].values, "tpm": v.values}).groupby("symbol")["tpm"].sum()
        present = by[by > 0]
        rows.append({
            "col": c,
            "cnhs": extract_cnhs(c),
            "n_genes_any": int((by > 0).sum()),
            "mean_sum": float(present.mean()) if len(present) else float("nan"),
            "mean_max": float("nan"),
        })
    return rows


def bird_fib_mean_max(genes: set[str]) -> pd.DataFrame:
    recs = []
    for sp in SPECIES_ORDER:
        meta = SPECIES[sp]
        exp = pd.read_csv(EXP / meta["file"], sep="\t")
        id2 = peak_to_symbol_de(RUSLAN / meta["de"])
        exp = annotate_table(exp, id2)
        cols = lib_cols(exp, "fib")
        sub = exp[exp["symbol"].isin(genes)].copy()
        if sub.empty or not cols:
            continue
        mat = sub[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        sub = sub.assign(**{c: mat[c] for c in cols})
        per_lib = []
        n_genes = []
        for c in cols:
            gmax = sub.groupby("symbol")[c].max()
            present = gmax[gmax > 0]
            n_genes.append(int(len(present)))
            per_lib.append(float(present.mean()) if len(present) else np.nan)
        recs.append({
            "species": sp.lower(),
            "clade": "bird",
            "cell": "embryonic_fibroblast",
            "source": "GSE213253_exp_tab",
            "mass_g": meta["mass_g"],
            "n_libs": len(cols),
            "n_genes_med": int(np.median(n_genes)) if n_genes else 0,
            "mean_tpm": float(np.nanmean(per_lib)),
            "agg": "max_peak",
        })
    return pd.DataFrame(recs)


def fit_block(df: pd.DataFrame, label: str) -> dict:
    y = np.log10(df["mean_tpm"].astype(float))
    fit = ols_logmass(pd.Series(y.values, index=df.index), df["mass_g"])
    sp = spearman_logmass(pd.Series(y.values, index=df.index), df["mass_g"])
    return {
        "panel": label,
        "n": fit["n"],
        "r": fit["r"],
        "p": fit["p"],
        "slope": fit["slope"],
        "rho": sp["rho"],
        "p_spearman": sp["p"],
        "species": ",".join(df["species"].tolist()),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = oxphos_targets()
    print(f"nuclear OXPHOS after alias collapse: {len(genes)}", flush=True)

    lib_rows = []

    # --- rat / dog official peak TPM ---
    for species, pref in [("rat", "rn6"), ("dog", "canFam3")]:
        tpm = load_osc_tpm(FANTOM / f"{pref}.cage_peak_tpm.osc.txt.gz")
        ann = pd.read_csv(FANTOM / f"{pref}.cage_peak_ann.txt.gz", sep="\t")
        id2 = map_rat_dog(ann)
        groups: dict[str, list[str]] = {}
        for c in tpm.columns:
            if c == "peak_id":
                continue
            g = classify_col(c)
            if g:
                groups.setdefault(g, []).append(c)
        for cell, cols in groups.items():
            scored = score_peak_table(tpm, id2, genes, cols)
            for col, rec in scored.items():
                lib_rows.append({
                    "species": species,
                    "cell": cell,
                    "source": "official_peak_tpm",
                    "col": col,
                    "agg_note": "peak_table",
                    **rec,
                })
        print(f"{species} groups: { {k: len(v) for k, v in groups.items()} }", flush=True)

    # --- chicken official peak TPM (SMC, hep; no MSC column) ---
    tx2 = chicken_tx2sym(FANTOM / "Gallus_gallus.Gallus_gallus-5.0.92.gtf.gz")
    print(f"chicken GTF symbols: {len(tx2)}", flush=True)
    tpm_gg = load_osc_tpm(FANTOM / "galGal5.cage_peak_tpm.osc.txt.gz")
    ann_gg = pd.read_csv(FANTOM / "galGal5.cage_peak_ann.txt.gz", sep="\t")
    id2_gg = map_chicken(ann_gg, tx2)
    groups_gg: dict[str, list[str]] = {}
    for c in tpm_gg.columns:
        if c == "peak_id":
            continue
        g = classify_col(c)
        if g:
            groups_gg.setdefault(g, []).append(c)
    for cell, cols in groups_gg.items():
        scored = score_peak_table(tpm_gg, id2_gg, genes, cols)
        for col, rec in scored.items():
            lib_rows.append({
                "species": "chicken",
                "cell": cell,
                "source": "official_peak_tpm",
                "col": col,
                "agg_note": "peak_table",
                **rec,
            })
    print(f"chicken groups: { {k: len(v) for k, v in groups_gg.items()} }", flush=True)

    msc_bw = FANTOM / "chicken_msc.fwd.bw"
    if msc_bw.exists():
        rec = chicken_msc_from_bw(id2_gg, genes, FANTOM / "galGal5.cage_peak_coord.bed.gz")
        lib_rows.append({
            "species": "chicken",
            "cell": "msc",
            "source": rec.pop("source"),
            "col": "CNhs11292_bigwig",
            "agg_note": "peak_sum_from_bw",
            **rec,
        })
        print(f"chicken MSC bigwig genes={rec['n_genes_any']} mean_max={rec['mean_max']:.3f}", flush=True)

    # --- human / mouse gene-level ---
    hg = FANTOM / "hg19.gene_phase1and2combined_tpm.osc.txt.gz"
    mm = FANTOM / "mm9.gene_phase1and2combined_tpm.osc.txt.gz"
    if hg.exists() and hg.stat().st_size > 200_000_000:
        for cell, ids in [("msc", HUMAN_MSC), ("smc", HUMAN_AORTIC), ("dermal_fib", HUMAN_DERMAL)]:
            for rec in score_gene_level(hg, ids, genes):
                lib_rows.append({
                    "species": "human",
                    "cell": cell,
                    "source": "official_gene_tpm",
                    "agg_note": "gene_sum",
                    **rec,
                })
        print("human gene-level done", flush=True)
    else:
        print(f"human gene TPM not ready ({hg} size={hg.stat().st_size if hg.exists() else 0})", flush=True)

    if mm.exists() and mm.stat().st_size > 100_000_000:
        for cell, ids in [("msc", MOUSE_MSC), ("smc", MOUSE_AORTIC), ("mef", MOUSE_MEF)]:
            for rec in score_gene_level(mm, ids, genes):
                lib_rows.append({
                    "species": "mouse",
                    "cell": cell,
                    "source": "official_gene_tpm",
                    "agg_note": "gene_sum",
                    **rec,
                })
        print("mouse gene-level done", flush=True)

    libs = pd.DataFrame(lib_rows)
    # Official chicken MSC peak TPM already contains CNhs11292; keep bigwig
    # only as a QC row in libraries.tsv, not in the species mean.
    libs.to_csv(OUT / "libraries.tsv", sep="\t", index=False)
    libs_fit = libs[~((libs["species"] == "chicken") & (libs["cell"] == "msc") & (libs["source"] == "datahub_tpm_bigwig"))]

    # species × cell: mean of library mean_sum (gene-level) or mean_max (peaks)
    # Headline for mixed sources: mean_sum when gene-level, mean_max when peak-level
    # so that each library uses the native FANTOM gene aggregation where possible.
    def pick_metric(r: pd.Series) -> float:
        # Prefer sum-of-peaks: that is FANTOM's gene-level aggregation
        # and is available for every library here.
        if np.isfinite(r.get("mean_sum", np.nan)):
            return float(r["mean_sum"])
        return float(r["mean_max"])

    libs_fit["metric"] = libs_fit.apply(pick_metric, axis=1)

    sp_rows = []
    for (sp, cell), g in libs_fit.groupby(["species", "cell"], sort=False):
        g = g[np.isfinite(g["metric"])]
        if g.empty:
            continue
        sp_rows.append({
            "species": sp,
            "cell": cell,
            "n_libs": len(g),
            "n_genes_med": int(g["n_genes_any"].median()),
            "mean_tpm": float(g["metric"].mean()),
            "source": ";".join(sorted(set(g["source"]))),
            "mass_g": MASS.get(sp, np.nan),
        })
    species = pd.DataFrame(sp_rows)
    species.to_csv(OUT / "species_cell.tsv", sep="\t", index=False)

    birds = bird_fib_mean_max(genes)
    birds.to_csv(OUT / "bird_fib_same_geneset.tsv", sep="\t", index=False)

    fits = []
    msc = species[species["cell"].eq("msc")].copy()
    if len(msc) >= 4:
        fits.append(fit_block(msc, "fantom_msc"))
        if "dog" in set(msc["species"]):
            alt = msc.copy()
            alt.loc[alt["species"].eq("dog"), "mass_g"] = MASS_DOG_ANAGE
            fits.append(fit_block(alt, "fantom_msc_dog40kg"))
    smc = species[species["cell"].eq("smc")].copy()
    if len(smc) >= 4:
        fits.append(fit_block(smc, "fantom_aortic_smc"))
    if len(birds) >= 4:
        fits.append(fit_block(birds.assign(mean_tpm=birds["mean_tpm"]), "gse213253_ef_maxpeak"))

    fit_df = pd.DataFrame(fits)
    fit_df.to_csv(OUT / "fits.tsv", sep="\t", index=False)
    print(species.to_string(index=False))
    print(fit_df.to_string(index=False))
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
