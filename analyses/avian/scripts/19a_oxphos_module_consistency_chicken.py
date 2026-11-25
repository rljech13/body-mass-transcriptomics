#!/usr/bin/env python3
"""Phase 1: CEF coexpression module + regulon-lite (prior TFs × FIMO × coupling)."""
from __future__ import annotations

import json
import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from oxphos_thread_lib import (
    CHICKEN_NC_TO_CHR,
    DATA,
    DRIVERS,
    DRIVERS_EXTRA,
    EXP,
    FIMO,
    GENOME_ROOT,
    RES,
    RUSLAN,
    annotate_table,
    gene_by_lib_matrix,
    lib_cols,
    load_geneset,
    median_pairwise_spearman,
    nuclear_oxphos,
    peak_to_symbol_de,
)

warnings.filterwarnings("ignore", category=UserWarning)

OUT = RES / "module_consistency"
DEF = RES / "definitions"
N_RANDOM = 200
RNG = np.random.default_rng(42)
FIMO_THRESH = "1e-4"
JASPAR_MEME_URL = (
    "https://jaspar.elixir.no/download/data/2024/CORE/"
    "JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt"
)
# JASPAR matrix IDs (vertebrates CORE)
JASPAR_IDS = {
    "NRF1": "MA0506.3",
    "GABPA": "MA0062.4",
    "ESRRA": "MA0592.4",
    "YY1": "MA0095.4",
}
JASPAR_NAME_ALIASES = {
    "NRF1": {"NRF1", "NRF-1", "NRF1"},
    "GABPA": {"GABPA"},
    "ESRRA": {"ESRRA", "ERR1", "NR3B1"},
    "YY1": {"YY1", "YY-1"},
}


def load_chicken() -> pd.DataFrame:
    path = EXP / "Chicken_fib_epi_exp.txt"
    df = pd.read_csv(path, sep="\t")
    id2 = peak_to_symbol_de(RUSLAN / "Chicken_DE_epi(up)_fib(down)_promoters.xlsx")
    return annotate_table(df, id2)


def background_genes(df: pd.DataFrame, tissue: str, thr: float = 1.0) -> pd.DataFrame:
    """Mean TPM across libs for promoter-assigned symbols (matched random pool)."""
    genes = set(
        df.loc[df["class"].astype(str).str.lower().eq("promoter") & df["symbol"].ne(""), "symbol"]
    )
    mat = gene_by_lib_matrix(df, genes, tissue, thr=thr)
    if mat.empty:
        return pd.DataFrame(columns=["gene", "mean_tpm"])
    mean_tpm = mat.mean(axis=1)
    mean_tpm = mean_tpm[mean_tpm > 0]
    return mean_tpm.rename("mean_tpm").rename_axis("gene").reset_index()


def coexpr_vs_random(mat: pd.DataFrame, bg: pd.DataFrame, label: str) -> dict:
    genes = [g for g in mat.index if g in set(bg.gene)]
    mat = mat.loc[genes]
    obs = median_pairwise_spearman(mat)
    # matched random: sample genes with similar mean TPM ranks
    bg = bg.set_index("gene")
    means = bg.loc[genes, "mean_tpm"].values
    lo, hi = np.percentile(means, [10, 90]) if len(means) else (0, 1)
    pool = bg[(bg.mean_tpm >= lo * 0.5) & (bg.mean_tpm <= hi * 2)].index.tolist()
    if len(pool) < len(genes) + 5:
        pool = bg.index.tolist()
    rand_medians = []
    n = len(genes)
    # need lib matrix for random genes — use mean_tpm ranks only via resampling from mat's parent?
    # We'll pass full gene matrix separately; here only bootstrap gene labels from mat if available
    return {
        "label": label,
        "n_genes": n,
        "median_rho": obs,
        "n_libs": int(mat.shape[1]),
        "pool_note": "random filled in caller",
    }


def random_coexpr(full_mat: pd.DataFrame, bg: pd.DataFrame, n_genes: int, n_rand: int = N_RANDOM) -> np.ndarray:
    """Return array of median pairwise Spearman for random gene sets."""
    pool = [g for g in bg.gene.tolist() if g in full_mat.index]
    if len(pool) < n_genes + 2:
        return np.array([])
    meds = []
    for _ in range(n_rand):
        pick = RNG.choice(pool, size=n_genes, replace=False)
        meds.append(median_pairwise_spearman(full_mat.loc[list(pick)]))
    return np.array([m for m in meds if np.isfinite(m)])


def build_full_gene_mat(df: pd.DataFrame, tissue: str, genes: list[str] | None = None) -> pd.DataFrame:
    if genes is None:
        genes = background_genes(df, tissue).gene.tolist()
    return gene_by_lib_matrix(df, set(genes), tissue)


def dominant_tss_windows(df: pd.DataFrame, genes: set[str], tissue: str = "fib") -> pd.DataFrame:
    cols = lib_cols(df, tissue)
    sub = df[df["class"].astype(str).str.lower().eq("promoter") & df["symbol"].isin(genes)].copy()
    for c in cols:
        sub[c] = pd.to_numeric(sub[c], errors="coerce").fillna(0.0)
    sub["_e"] = sub[cols].mean(axis=1)
    rows = []
    for gene, gdf in sub.groupby("symbol"):
        maj = gdf.groupby("strand")["_e"].sum().idxmax()
        gdf = gdf[gdf.strand == maj]
        gdf = gdf[gdf["_e"] >= 1.0]
        if gdf.empty:
            continue
        peak = gdf.loc[gdf["_e"].idxmax()]
        chrom = str(peak["seqnames"])
        strand = str(peak["strand"])
        start = int(peak["start"])
        end = int(peak["end"])
        # TSS edge: 5' of peak
        if strand == "+":
            tss = start
            w0, w1 = tss - 300, tss + 100
        else:
            tss = end
            w0, w1 = tss - 100, tss + 300
        rows.append(
            dict(
                gene=gene,
                chrom_raw=chrom,
                chrom=CHICKEN_NC_TO_CHR.get(chrom, chrom),
                strand=strand,
                tss=tss,
                win_start=max(0, w0),
                win_end=w1,
                peak_tpm=float(peak["_e"]),
            )
        )
    return pd.DataFrame(rows)


def fetch_jaspar_meme(out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 1000:
        return out_path
    import urllib.request

    urllib.request.urlretrieve(JASPAR_MEME_URL, out_path)
    return out_path


def extract_selected_meme(full_meme: Path, out_meme: Path) -> Path:
    text = full_meme.read_text(errors="replace")
    blocks = text.split("MOTIF ")
    header = blocks[0]
    keep = []
    wanted_ids = set(JASPAR_IDS.values())
    wanted_names = {n.upper() for names in JASPAR_NAME_ALIASES.values() for n in names}
    # also accept Nrf1 / Yy1 casing from JASPAR
    wanted_names |= {"NRF1", "GABPA", "ESRRA", "YY1"}
    for b in blocks[1:]:
        first = b.splitlines()[0].strip()
        parts = first.split()
        mid = parts[0]
        name = parts[1].upper() if len(parts) > 1 else ""
        if mid in wanted_ids or name in wanted_names:
            keep.append("MOTIF " + b)
    if not keep:
        raise RuntimeError("No JASPAR motifs extracted for OXPHOS drivers")
    out_meme.write_text(header + "".join(keep))
    return out_meme


def run_fimo_on_fasta(meme: Path, fasta: Path, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / "fimo.tsv"
    cmd = [
        str(FIMO),
        "--thresh",
        FIMO_THRESH,
        "--oc",
        str(out_dir),
        str(meme),
        str(fasta),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    if not tsv.exists():
        return pd.DataFrame()
    # skip comment lines
    rows = []
    with tsv.open() as fh:
        header = None
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            if header is None:
                header = ln.strip().split("\t")
                continue
            rows.append(ln.strip().split("\t"))
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows, columns=header)


def fasta_from_windows(wins: pd.DataFrame, genome_fa: Path, out_fa: Path) -> int:
    # index genome with pyfaidx if available else simple
    try:
        from pyfaidx import Fasta
    except ImportError:
        subprocess.run(
            ["pip", "install", "pyfaidx", "-q"],
            check=False,
            capture_output=True,
        )
        from pyfaidx import Fasta

    fa = Fasta(str(genome_fa))
    n = 0
    with out_fa.open("w") as out:
        for r in wins.itertuples(index=False):
            chrom = r.chrom
            if chrom not in fa:
                # try raw
                chrom = r.chrom_raw
            if chrom not in fa:
                continue
            s = int(r.win_start)
            e = int(r.win_end)
            seq = fa[chrom][s:e]
            if r.strand == "-":
                seq = seq.complement.reverse
            seq = str(seq).upper()
            if len(seq) < 50:
                continue
            out.write(f">{r.gene}\n{seq}\n")
            n += 1
    return n


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    DEF.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)

    df = load_chicken()
    ox_set = nuclear_oxphos()
    glyc = load_geneset("GLYCOLYSIS")

    # --- set freeze ---
    set_df = pd.DataFrame({"gene": sorted(ox_set), "definition": "set", "source": "OXPHOS.txt_nuclear"})
    set_df.to_csv(DEF / "oxphos_set.tsv", sep="\t", index=False)

    # --- CEF fib coexpr ---
    bg_fib = background_genes(df, "fib")
    ox_mat = gene_by_lib_matrix(df, ox_set, "fib")
    ox_mat = ox_mat.loc[(ox_mat > 0).any(axis=1)]
    glyc_mat = gene_by_lib_matrix(df, glyc, "fib")
    glyc_mat = glyc_mat.loc[(glyc_mat > 0).any(axis=1)]

    # full mat for random: top 800 expressed genes (speed)
    bg_top = bg_fib.nlargest(min(800, len(bg_fib)), "mean_tpm")
    print(f"[19a] building full mat n_genes={len(bg_top)} ...", flush=True)
    full_mat = gene_by_lib_matrix(df, set(bg_top.gene), "fib")
    print(f"[19a] full_mat {full_mat.shape}", flush=True)

    ox_rho = median_pairwise_spearman(ox_mat)
    glyc_rho = median_pairwise_spearman(glyc_mat)
    rand_meds = random_coexpr(full_mat, bg_top, n_genes=len(ox_mat), n_rand=N_RANDOM)
    p_ox = float(np.mean(rand_meds >= ox_rho)) if len(rand_meds) and np.isfinite(ox_rho) else np.nan

    coexpr_rows = [
        dict(group="OXPHOS_set", n_genes=len(ox_mat), n_libs=ox_mat.shape[1], median_rho=ox_rho, random_median_rho=float(np.median(rand_meds)) if len(rand_meds) else np.nan, empir_p_ge_obs=p_ox),
        dict(group="GLYCOLYSIS", n_genes=len(glyc_mat), n_libs=glyc_mat.shape[1], median_rho=glyc_rho, random_median_rho=np.nan, empir_p_ge_obs=np.nan),
    ]
    pd.DataFrame(coexpr_rows).to_csv(OUT / "cef_coexpr_vs_random.tsv", sep="\t", index=False)

    # Module: genes with mean pairwise rho to other OXPHOS genes above random median (or keep all if module ok)
    # Per-gene: median Spearman to other OXPHOS genes
    # per-gene median corr to others via one Spearman corr matrix
    m = ox_mat.loc[ox_mat.std(axis=1) > 1e-12]
    ranks = m.rank(axis=1, method="average")
    cmat = np.corrcoef(ranks.to_numpy(dtype=float))
    np.fill_diagonal(cmat, np.nan)
    gene_scores = []
    for i, g in enumerate(m.index):
        vals = cmat[i]
        vals = vals[np.isfinite(vals)]
        gene_scores.append(
            dict(
                gene=g,
                median_rho_to_others=float(np.median(vals)) if len(vals) else np.nan,
                mean_tpm=float(ox_mat.loc[g].mean()),
            )
        )
    gsc = pd.DataFrame(gene_scores).sort_values("median_rho_to_others", ascending=False)
    gsc.to_csv(OUT / "cef_oxphos_per_gene_coexpr.tsv", sep="\t", index=False)

    rand_med = float(np.median(rand_meds)) if len(rand_meds) else 0.0
    # keep genes with median_rho_to_others >= max(0, rand_med) OR if overall module fails, keep top 50% by score
    module_ok = bool(np.isfinite(ox_rho) and np.isfinite(p_ox) and ox_rho > rand_med and p_ox < 0.05)
    if module_ok:
        mod_genes = gsc.loc[gsc.median_rho_to_others >= rand_med, "gene"].tolist()
        if len(mod_genes) < 10:
            mod_genes = gsc.head(max(15, len(gsc) // 2)).gene.tolist()
    else:
        mod_genes = gsc.head(max(15, len(gsc) // 2)).gene.tolist()

    mod_df = pd.DataFrame({"gene": sorted(mod_genes), "definition": "module", "source": "CEF_coexpr_filter"})
    mod_df.to_csv(DEF / "oxphos_module.tsv", sep="\t", index=False)

    # recompute module coexpr
    mod_mat = ox_mat.loc[[g for g in mod_genes if g in ox_mat.index]]
    mod_rho = median_pairwise_spearman(mod_mat)
    rand_mod = random_coexpr(full_mat, bg_top, n_genes=len(mod_mat), n_rand=N_RANDOM)
    p_mod = float(np.mean(rand_mod >= mod_rho)) if len(rand_mod) and np.isfinite(mod_rho) else np.nan
    pd.DataFrame(
        [
            dict(group="OXPHOS_module", n_genes=len(mod_mat), median_rho=mod_rho, random_median_rho=float(np.median(rand_mod)) if len(rand_mod) else np.nan, empir_p_ge_obs=p_mod, module_ok_set_level=module_ok),
        ]
    ).to_csv(OUT / "cef_module_coexpr.tsv", sep="\t", index=False)

    # --- drivers coupling ---
    drivers = [d for d in DRIVERS + DRIVERS_EXTRA]
    dmat = gene_by_lib_matrix(df, set(drivers), "fib")
    coup_rows = []
    module_mean = mod_mat.mean(axis=0)
    for d in drivers:
        if d not in dmat.index:
            coup_rows.append(dict(driver=d, spearman_vs_module_mean=np.nan, n_libs=mod_mat.shape[1], detected=False))
            continue
        r, p = stats.spearmanr(dmat.loc[d], module_mean)
        coup_rows.append(dict(driver=d, spearman_vs_module_mean=float(r), p=float(p), n_libs=mod_mat.shape[1], detected=True))
    coup = pd.DataFrame(coup_rows)
    coup.to_csv(OUT / "cef_driver_module_coupling.tsv", sep="\t", index=False)

    # per-gene TF coupling for prior drivers
    per_tf = []
    for d in DRIVERS:
        if d not in dmat.index:
            continue
        for g in mod_mat.index:
            r, p = stats.spearmanr(dmat.loc[d], mod_mat.loc[g])
            per_tf.append(dict(driver=d, gene=g, spearman=float(r) if np.isfinite(r) else np.nan, p=float(p) if np.isfinite(p) else np.nan))
    pd.DataFrame(per_tf).to_csv(OUT / "cef_driver_per_gene_coupling.tsv", sep="\t", index=False)

    # --- FIMO on chicken module promoters ---
    wins = dominant_tss_windows(df, set(mod_genes), "fib")
    wins.to_csv(OUT / "chicken_module_tss_windows.tsv", sep="\t", index=False)
    genome = GENOME_ROOT / "Chicken" / "genome.fa"
    motif_dir = OUT / "fimo_chicken"
    motif_dir.mkdir(parents=True, exist_ok=True)
    full_meme = motif_dir / "JASPAR2024_vertebrates.meme"
    sel_meme = motif_dir / "oxphos_drivers.meme"
    fasta = motif_dir / "module_promoters.fa"
    fimo_hits = pd.DataFrame()
    fimo_ok = False
    try:
        fetch_jaspar_meme(full_meme)
        extract_selected_meme(full_meme, sel_meme)
        n_fa = fasta_from_windows(wins, genome, fasta)
        if n_fa > 0 and FIMO.exists():
            fimo_hits = run_fimo_on_fasta(sel_meme, fasta, motif_dir / "run")
            fimo_ok = not fimo_hits.empty
            fimo_hits.to_csv(OUT / "fimo_hits_raw.tsv", sep="\t", index=False)
    except Exception as e:
        (OUT / "fimo_error.txt").write_text(str(e))

    # map motif id → TF
    id2tf = {v: k for k, v in JASPAR_IDS.items()}
    genes_with_motif = set()
    if fimo_ok:
        # columns: motif_id, motif_alt_id, sequence_name, ...
        col_seq = "sequence_name" if "sequence_name" in fimo_hits.columns else fimo_hits.columns[2]
        col_motif = "motif_id" if "motif_id" in fimo_hits.columns else fimo_hits.columns[0]
        for _, r in fimo_hits.iterrows():
            genes_with_motif.add(str(r[col_seq]))
        # also per TF
        hit_summary = []
        for mid, tf in id2tf.items():
            sub = fimo_hits[fimo_hits[col_motif].astype(str).str.contains(mid.split(".")[0]) | fimo_hits[col_motif].eq(mid)]
            # JASPAR often uses MA0506.1 as motif_id
            sub = fimo_hits[fimo_hits[col_motif].astype(str) == mid]
            hit_summary.append(dict(tf=tf, motif_id=mid, n_genes_with_hit=sub[col_seq].nunique() if len(sub) else 0))
        pd.DataFrame(hit_summary).to_csv(OUT / "fimo_hit_summary.tsv", sep="\t", index=False)

    # regulon-lite: module genes with motif OR positive corr to any prior driver
    pos_corr_genes = set()
    if per_tf:
        ptf = pd.DataFrame(per_tf)
        pos_corr_genes = set(ptf.loc[ptf.spearman > 0, "gene"])

    if genes_with_motif:
        regulon_genes = sorted(set(mod_genes) & (genes_with_motif | pos_corr_genes))
        method = "lite_fimo_or_poscorr"
    else:
        regulon_genes = sorted(set(mod_genes) & pos_corr_genes) or sorted(mod_genes)
        method = "lite_poscorr_only_fimo_failed" if not fimo_ok else "lite_poscorr_only"

    reg_df = pd.DataFrame({"gene": regulon_genes, "definition": "regulon", "source": method})
    reg_df.to_csv(DEF / "oxphos_regulon.tsv", sep="\t", index=False)
    (DEF / "oxphos_driver_tfs.txt").write_text("\n".join(DRIVERS) + "\n")

    # chicken epi coherence check for three defs
    epi_rows = []
    for name, genes in [("set", ox_set), ("module", set(mod_genes)), ("regulon", set(regulon_genes))]:
        mat = gene_by_lib_matrix(df, genes, "epi")
        mat = mat.loc[(mat > 0).any(axis=1)]
        rho = median_pairwise_spearman(mat)
        epi_rows.append(dict(definition=name, n_genes_detected=len(mat), n_libs=mat.shape[1], median_rho=rho, mean_tpm=float(mat.mean().mean()) if len(mat) else np.nan))
    pd.DataFrame(epi_rows).to_csv(OUT / "chicken_epi_coherence_three_defs.tsv", sep="\t", index=False)

    # also fib coherence three defs
    fib_rows = []
    for name, genes in [("set", ox_set), ("module", set(mod_genes)), ("regulon", set(regulon_genes))]:
        mat = gene_by_lib_matrix(df, genes, "fib")
        mat = mat.loc[(mat > 0).any(axis=1)]
        rho = median_pairwise_spearman(mat)
        fib_rows.append(dict(definition=name, n_genes_detected=len(mat), n_libs=mat.shape[1], median_rho=rho, mean_tpm=float(mat.mean().mean()) if len(mat) else np.nan))
    pd.DataFrame(fib_rows).to_csv(OUT / "chicken_fib_coherence_three_defs.tsv", sep="\t", index=False)

    summary = {
        "module_ok_set_level": module_ok,
        "oxphos_set_median_rho": ox_rho,
        "random_median_rho": float(np.median(rand_meds)) if len(rand_meds) else None,
        "empir_p": p_ox,
        "n_set_detected_fib": len(ox_mat),
        "n_module": len(mod_genes),
        "module_median_rho": mod_rho,
        "module_empir_p": p_mod,
        "n_regulon": len(regulon_genes),
        "regulon_method": method,
        "fimo_ok": fimo_ok,
        "n_genes_with_motif": len(genes_with_motif),
        "driver_coupling": coup.to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
