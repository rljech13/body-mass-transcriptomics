#!/usr/bin/env python3
"""Phase 2 cis: GC / CpG o/e / prior-motif density in OXPHOS TSS windows across n=7.

Not motif discovery — scalar X vs log10(mass) for frozen module/regulon (+ set, GLYCOLYSIS).

Sequence source:
  - Chicken: local genomes_with_mt + galGal6 NC→chr alias
  - Other spp with NCBI-style chroms (NC_/CM_/HG_/OZ_): NCBI sequence viewer
  - Peacock (local scaffolds only): skipped (no fetchable accession)
"""
from __future__ import annotations

import json
import re
import subprocess
import time
import urllib.error
import urllib.request
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from oxphos_thread_lib import (
    CHICKEN_NC_TO_CHR,
    EXP,
    FIG,
    FIMO,
    GENOME_ROOT,
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    load_geneset,
    ols_logmass,
    peak_to_symbol_de,
)

warnings.filterwarnings("ignore", category=UserWarning)

OUT = RES / "cis"
DEF = RES / "definitions"
MEME_SHUFFLE = Path("/path/to/user-home/miniconda3/envs/meme/bin/fasta-dinucleotide-shuffle")
FIMO_THRESH = "1e-4"
JASPAR_IDS = {
    "NRF1": "MA0506.3",
    "GABPA": "MA0062.4",
    "ESRRA": "MA0592.4",
    "YY1": "MA0095.4",
}
JASPAR_MEME = (
    RES / "module_consistency" / "fimo_chicken" / "oxphos_drivers.meme"
)
WIN_UP = 300
WIN_DOWN = 100
NCBI_SLEEP = 0.34


def load_def(name: str) -> set[str]:
    return set(pd.read_csv(DEF / f"oxphos_{name}.tsv", sep="\t")["gene"].astype(str).str.upper())


def load_fai(path: Path) -> dict[str, int]:
    out = {}
    for ln in path.read_text().splitlines():
        c, L, *_ = ln.split("\t")
        out[c] = int(L)
    return out


def load_galgal_alias() -> dict[str, str]:
    alias_path = Path("/path/to/storage/kleiber_birds_gemcat/ref/galGal6.chromAlias.txt")
    m = dict(CHICKEN_NC_TO_CHR)
    if not alias_path.exists():
        return m
    for ln in alias_path.read_text().splitlines():
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) < 4:
            continue
        ucsc, _ens, _gb, refseq = parts[0], parts[1], parts[2], parts[3]
        m[refseq] = ucsc
        m[refseq.split(".")[0]] = ucsc
    return m


def dominant_tss_windows(df: pd.DataFrame, genes: set[str], tissue: str = "fib") -> pd.DataFrame:
    cols = [c for c in df.columns if str(c).startswith("norm.C" if tissue == "fib" else "norm.E")]
    sub = df[df["class"].astype(str).str.lower().eq("promoter") & df["symbol"].isin(genes)].copy()
    if sub.empty or not cols:
        return pd.DataFrame()
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
        if strand == "+":
            tss = start
            w0, w1 = tss - WIN_UP, tss + WIN_DOWN
        else:
            tss = end
            w0, w1 = tss - WIN_DOWN, tss + WIN_UP
        rows.append(
            dict(
                gene=gene,
                chrom_raw=chrom,
                strand=strand,
                tss=tss,
                win_start=max(0, w0),
                win_end=w1,
                peak_tpm=float(peak["_e"]),
            )
        )
    return pd.DataFrame(rows)


def resolve_chrom(species: str, chrom_raw: str, alias: dict[str, str], fai: dict[str, int]) -> str | None:
    if chrom_raw in fai:
        return chrom_raw
    if species == "Chicken":
        mapped = alias.get(chrom_raw) or alias.get(chrom_raw.split(".")[0])
        if mapped and mapped in fai:
            return mapped
    return None


def fetch_local(genome: Path, chrom: str, start: int, end: int, strand: str) -> str:
    region = f"{chrom}:{start + 1}-{end}"  # faidx is 1-based inclusive
    r = subprocess.run(
        ["samtools", "faidx", str(genome), region],
        capture_output=True,
        text=True,
        check=False,
    )
    if r.returncode != 0 or not r.stdout.strip():
        return ""
    seq = "".join(ln.strip() for ln in r.stdout.splitlines() if not ln.startswith(">"))
    seq = seq.upper()
    if strand == "-":
        comp = str.maketrans("ACGTN", "TGCAN")
        seq = seq.translate(comp)[::-1]
    return seq


def fetch_ncbi(chrom: str, start: int, end: int, strand: str) -> str:
    """1-based inclusive window via NCBI sequence viewer; strand applied locally."""
    # GenBank chroms are CM######## / HG###### (no underscore); RefSeq NC_/NW_/NZ_
    if not re.match(r"^(NC_|NW_|NZ_|CM|HG|OZ|CP|JA)", chrom):
        return ""
    # viewer from/to are 1-based
    fr, to = start + 1, end
    if fr < 1:
        fr = 1
    if to < fr:
        return ""
    url = (
        "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?"
        f"id={chrom}&db=nuccore&report=fasta&from={fr}&to={to}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "oxphos-cis/1.0"})
    for attempt in range(4):
        try:
            time.sleep(NCBI_SLEEP + 0.15 * attempt)
            with urllib.request.urlopen(req, timeout=60) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            if "<html" in text[:200].lower():
                time.sleep(2.0 * (attempt + 1))
                continue
            seq = "".join(ln.strip() for ln in text.splitlines() if not ln.startswith(">") and ln.strip())
            seq = re.sub(r"[^ACGTacgtNn]", "", seq).upper()
            if len(seq) < 50:
                time.sleep(1.0 * (attempt + 1))
                continue
            if strand == "-":
                comp = str.maketrans("ACGTN", "TGCAN")
                seq = seq.translate(comp)[::-1]
            return seq
        except (urllib.error.URLError, TimeoutError, ValueError):
            time.sleep(1.5 * (attempt + 1))
    return ""


def gc_frac(seq: str) -> float:
    s = re.sub(r"[^ACGT]", "", seq.upper())
    if not s:
        return float("nan")
    return (s.count("G") + s.count("C")) / len(s)


def cpg_oe(seq: str) -> float:
    s = re.sub(r"[^ACGT]", "", seq.upper())
    L = len(s)
    if L < 2:
        return float("nan")
    n_c = s.count("C")
    n_g = s.count("G")
    n_cg = s.count("CG")
    denom = (n_c / L) * (n_g / L)
    if denom <= 0:
        return float("nan")
    return (n_cg / L) / denom


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def dinuc_shuffle_fasta(infile: Path, outfile: Path, seed: int = 42) -> None:
    cmd = [str(MEME_SHUFFLE), "-f", str(infile), "-s", str(seed)]
    r = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if r.returncode != 0:
        raise RuntimeError(r.stderr or r.stdout)
    outfile.write_text(r.stdout)


def run_fimo(meme: Path, fasta: Path, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(FIMO), "--thresh", FIMO_THRESH, "--oc", str(out_dir), str(meme), str(fasta)]
    r = subprocess.run(cmd, capture_output=True, text=True, check=False)
    tsv = out_dir / "fimo.tsv"
    if not tsv.exists():
        (out_dir / "fimo_stderr.txt").write_text(r.stderr or r.stdout or "")
        return pd.DataFrame()
    rows, header = [], None
    with tsv.open() as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            if header is None:
                header = ln.rstrip("\n").split("\t")
                continue
            rows.append(ln.rstrip("\n").split("\t"))
    if not rows:
        return pd.DataFrame(columns=header or [])
    return pd.DataFrame(rows, columns=header)


def motif_density_from_fimo(hits: pd.DataFrame, seq_len: int = 400) -> pd.DataFrame:
    """Per sequence_name × TF: hits per 100 bp."""
    if hits.empty:
        return pd.DataFrame(columns=["seq_id", "tf", "n_hits", "density_per_100bp"])
    col_seq = "sequence_name" if "sequence_name" in hits.columns else hits.columns[2]
    col_motif = "motif_id" if "motif_id" in hits.columns else hits.columns[0]
    id2tf = {v: k for k, v in JASPAR_IDS.items()}
    # also match without version
    for k, v in list(id2tf.items()):
        id2tf[k.split(".")[0]] = v
    rows = []
    for seq_id, g in hits.groupby(col_seq):
        for mid, gg in g.groupby(col_motif):
            tf = id2tf.get(str(mid), id2tf.get(str(mid).split(".")[0], str(mid)))
            n = len(gg)
            rows.append(
                dict(
                    seq_id=str(seq_id),
                    tf=tf,
                    n_hits=n,
                    density_per_100bp=100.0 * n / seq_len,
                )
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    alias = load_galgal_alias()
    defs = {
        "module": load_def("module"),
        "regulon": load_def("regulon"),
        "set": load_def("set"),
        "GLYCOLYSIS": load_geneset("GLYCOLYSIS"),
    }
    masses = pd.Series({sp: SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER})

    if not JASPAR_MEME.exists():
        raise SystemExit(f"Missing meme file {JASPAR_MEME}; run 19a first")

    win_rows = []
    seq_records: list[tuple[str, str]] = []
    fetch_log = []

    for sp in SPECIES_ORDER:
        meta = SPECIES[sp]
        df = pd.read_csv(EXP / meta["file"], sep="\t")
        id2 = peak_to_symbol_de(RUSLAN / meta["de"])
        df = annotate_table(df, id2)
        genome = GENOME_ROOT / sp / "genome.fa"
        fai = load_fai(genome.with_suffix(".fa.fai")) if genome.exists() else {}
        for dname, genes in defs.items():
            wins = dominant_tss_windows(df, genes, tissue="fib")
            if wins.empty:
                fetch_log.append(dict(species=sp, definition=dname, status="no_windows", n_ok=0, n_fail=0))
                continue
            n_ok = n_fail = 0
            for r in wins.itertuples(index=False):
                chrom_local = resolve_chrom(sp, r.chrom_raw, alias, fai)
                seq = ""
                source = ""
                if chrom_local and genome.exists():
                    seq = fetch_local(genome, chrom_local, int(r.win_start), int(r.win_end), r.strand)
                    source = "local"
                if not seq:
                    seq = fetch_ncbi(r.chrom_raw, int(r.win_start), int(r.win_end), r.strand)
                    source = "ncbi" if seq else "fail"
                status = "ok" if seq else "fail"
                if seq:
                    n_ok += 1
                    sid = f"{sp}|{dname}|{r.gene}"
                    seq_records.append((sid, seq))
                    win_rows.append(
                        dict(
                            species=sp,
                            definition=dname,
                            gene=r.gene,
                            chrom_raw=r.chrom_raw,
                            chrom_resolved=chrom_local or "",
                            strand=r.strand,
                            tss=r.tss,
                            win_start=r.win_start,
                            win_end=r.win_end,
                            peak_tpm=r.peak_tpm,
                            seq_len=len(seq),
                            GC=gc_frac(seq),
                            CpG_oe=cpg_oe(seq),
                            source=source,
                            status=status,
                        )
                    )
                else:
                    n_fail += 1
                    win_rows.append(
                        dict(
                            species=sp,
                            definition=dname,
                            gene=r.gene,
                            chrom_raw=r.chrom_raw,
                            chrom_resolved=chrom_local or "",
                            strand=r.strand,
                            tss=r.tss,
                            win_start=r.win_start,
                            win_end=r.win_end,
                            peak_tpm=r.peak_tpm,
                            seq_len=0,
                            GC=np.nan,
                            CpG_oe=np.nan,
                            source=source,
                            status=status,
                        )
                    )
            fetch_log.append(dict(species=sp, definition=dname, status="done", n_ok=n_ok, n_fail=n_fail))

    wins_df = pd.DataFrame(win_rows)
    wins_df.to_csv(OUT / "promoter_windows_features.tsv", sep="\t", index=False)
    pd.DataFrame(fetch_log).to_csv(OUT / "fetch_status.tsv", sep="\t", index=False)

    fasta = OUT / "promoters_real.fa"
    write_fasta(fasta, seq_records)
    scramble_fa = OUT / "promoters_scramble.fa"
    if seq_records:
        dinuc_shuffle_fasta(fasta, scramble_fa)
    else:
        scramble_fa.write_text("")

    hits_real = run_fimo(JASPAR_MEME, fasta, OUT / "fimo_real") if seq_records else pd.DataFrame()
    hits_scr = run_fimo(JASPAR_MEME, scramble_fa, OUT / "fimo_scramble") if seq_records else pd.DataFrame()
    hits_real.to_csv(OUT / "fimo_hits_real.tsv", sep="\t", index=False)
    hits_scr.to_csv(OUT / "fimo_hits_scramble.tsv", sep="\t", index=False)

    dens_real = motif_density_from_fimo(hits_real)
    dens_scr = motif_density_from_fimo(hits_scr)
    dens_real["kind"] = "real"
    dens_scr["kind"] = "scramble"
    dens = pd.concat([dens_real, dens_scr], ignore_index=True)
    dens.to_csv(OUT / "motif_density_long.tsv", sep="\t", index=False)

    # attach total motif density to window table
    tot_real = (
        dens_real.groupby("seq_id")["density_per_100bp"].sum()
        if not dens_real.empty
        else pd.Series(dtype=float)
    )
    tot_scr = (
        dens_scr.groupby("seq_id")["density_per_100bp"].sum()
        if not dens_scr.empty
        else pd.Series(dtype=float)
    )
    wins_df["seq_id"] = wins_df["species"] + "|" + wins_df["definition"] + "|" + wins_df["gene"]
    wins_df["motif_density"] = wins_df["seq_id"].map(tot_real)
    wins_df["motif_density_scramble"] = wins_df["seq_id"].map(tot_scr)
    wins_df["motif_density_delta"] = wins_df["motif_density"].fillna(0) - wins_df["motif_density_scramble"].fillna(0)
    # genes with no hits → 0 density if sequence ok
    ok = wins_df.status.eq("ok")
    wins_df.loc[ok & wins_df.motif_density.isna(), "motif_density"] = 0.0
    wins_df.loc[ok & wins_df.motif_density_scramble.isna(), "motif_density_scramble"] = 0.0
    wins_df.loc[ok, "motif_density_delta"] = (
        wins_df.loc[ok, "motif_density"] - wins_df.loc[ok, "motif_density_scramble"]
    )
    wins_df.to_csv(OUT / "promoter_windows_features.tsv", sep="\t", index=False)

    # per-TF densities wide
    if not dens_real.empty:
        wide = dens_real.pivot_table(index="seq_id", columns="tf", values="density_per_100bp", aggfunc="sum").fillna(0)
        wide.columns = [f"dens_{c}" for c in wide.columns]
        wins_df = wins_df.merge(wide, left_on="seq_id", right_index=True, how="left")
        for c in wide.columns:
            wins_df.loc[ok & wins_df[c].isna(), c] = 0.0
        wins_df.to_csv(OUT / "promoter_windows_features.tsv", sep="\t", index=False)

    # species medians
    metrics = ["GC", "CpG_oe", "motif_density", "motif_density_delta"]
    for tf in JASPAR_IDS:
        col = f"dens_{tf}"
        if col in wins_df.columns:
            metrics.append(col)

    agg_rows = []
    ols_rows = []
    for dname in defs:
        for metric in metrics:
            if metric not in wins_df.columns:
                continue
            med = {}
            for sp in SPECIES_ORDER:
                sub = wins_df[(wins_df.species == sp) & (wins_df.definition == dname) & (wins_df.status == "ok")]
                med[sp] = float(sub[metric].median()) if len(sub) else np.nan
                agg_rows.append(
                    dict(
                        definition=dname,
                        species=sp,
                        mass_g=masses[sp],
                        metric=metric,
                        median=med[sp],
                        n_genes=int(len(sub)),
                    )
                )
            fit = ols_logmass(pd.Series(med).reindex(SPECIES_ORDER), masses)
            ols_rows.append(dict(definition=dname, metric=metric, **fit))

    agg = pd.DataFrame(agg_rows)
    ols = pd.DataFrame(ols_rows)
    agg.to_csv(OUT / "species_median_scalars.tsv", sep="\t", index=False)
    ols.to_csv(OUT / "ols_cis_vs_log10mass.tsv", sep="\t", index=False)

    # figure: module GC / CpG / motif_density vs mass
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), constrained_layout=True)
    x = np.log10(masses.reindex(SPECIES_ORDER).values)
    for ax, metric, ylab in zip(
        axes,
        ["GC", "CpG_oe", "motif_density"],
        ["median GC", "median CpG o/e", "median motif dens /100bp"],
    ):
        sub = agg[(agg.definition == "module") & (agg.metric == metric)].set_index("species").reindex(SPECIES_ORDER)
        y = sub["median"]
        ax.scatter(x, y.values, c="#1d3557", s=40)
        fit = ols_logmass(y, masses.reindex(SPECIES_ORDER))
        if np.isfinite(fit["slope"]):
            xx = np.linspace(x.min(), x.max(), 50)
            ax.plot(xx, fit["intercept"] + fit["slope"] * xx, color="#e76f51", lw=1)
        for sp, xi, yi in zip(SPECIES_ORDER, x, y.values):
            if np.isfinite(yi):
                ax.text(xi, yi, sp[:3], fontsize=7)
        ax.set_title(f"module {metric}\np={fit['p']:.3g} n={fit['n']}")
        ax.set_xlabel("log10 mass")
        ax.set_ylabel(ylab)
    fig.suptitle("OXPHOS module cis scalars vs mass (fib TSS −300/+100)", fontsize=11)
    fig.savefig(FIG / "fig_oxphos_thread_cis_module_vs_mass.png", dpi=160)
    fig.savefig(OUT / "fig_oxphos_thread_cis_module_vs_mass.png", dpi=160)
    plt.close()

    summary = {
        "window": "[-300,+100] from dominant fib TSS",
        "fimo_thresh": FIMO_THRESH,
        "jaspar_ids": JASPAR_IDS,
        "n_sequences": len(seq_records),
        "fetch_log": fetch_log,
        "ols_module": ols[ols.definition == "module"].to_dict(orient="records"),
        "note": (
            "Local genomes_with_mt mostly differ from Dropbox CAGE assemblies; "
            "Chicken uses local+alias; others use NCBI fetch by CAGE chrom accession; "
            "Peacock scaffolds are not fetchable."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(ols[ols.definition.isin(["module", "regulon"])].to_string(index=False))
    print("\nFetch:")
    print(pd.DataFrame(fetch_log).to_string(index=False))


if __name__ == "__main__":
    main()
