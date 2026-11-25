#!/usr/bin/env python3
"""Motif enrichment for the nine locked bird sets.

Does not rewrite story_cre_motifs/*_fimo_proximal.tsv, the HOCOMOCO
hit tables, or method_compare/. Those scans already cover every
proximal peak. This script only changes which genes are the foreground.

JASPAR FIMO and HOCOMOCO FIMO: two-sided Fisher, BH within set x species.
SEA and monaLisa: same settings as script 82, new output directory.
A motif is reported when BH < 0.05 in at least three species.
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
SCRIPTS = ROOT / "scripts"
STORY = ROOT / "results" / "oxphos_species_thread" / "story_cre_motifs"
OUT = STORY / "nine_sets"
KB = ROOT / "corpus_kleiber_2026-08-06" / "02_база_знаний"
HOCO = STORY / "method_compare" / "hocomoco"
JASPAR_MEME = STORY / "method_compare" / "JASPAR2024_one_per_symbol.meme"
SEA_BIN = STORY / "method_compare" / "meme55" / "bin" / "sea"
R_SCRIPT = SCRIPTS / "94_monalisa_nine.R"

SETS = [
    ("oxphos_16", "OXPHOS 16"),
    ("folate_12", "фолат 12"),
    ("ribo_5", "ribo 5"),
    ("green_hubs_7", "green hubs 7"),
    ("green_51_rest", "green, остальные 41"),
    ("carbon_7", "углерод 7"),
    ("amino_5", "аминокислоты 5"),
    ("nafld_4", "NAFLD 4"),
    ("epi_lipid_3", "липиды 3"),
]


def load_s82():
    spec = importlib.util.spec_from_file_location("motif82", SCRIPTS / "82_motif_method_compare.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_sets() -> dict[str, set[str]]:
    dis = pd.read_csv(KB / "working_sets_disjoint_2026-09-22.tsv", sep="\t")
    ov = pd.read_csv(KB / "working_sets_overlap_ok_2026-09-22.tsv", sep="\t")
    genes = {stem: set(dis.loc[dis.set_id == stem, "gene"].astype(str)) for stem, _ in SETS}
    for stem, _ in SETS:
        extra = set(ov.loc[ov.set_id == stem, "gene"].astype(str))
        if extra:
            genes[stem] = extra
    gene_dir = OUT / "genes"
    gene_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for stem, title in SETS:
        if not genes[stem]:
            raise SystemExit(f"пустой набор {stem}")
        path = gene_dir / f"{stem}.txt"
        path.write_text("".join(g + "\n" for g in sorted(genes[stem])))
        rows.append(dict(set=stem, title=title, n=len(genes[stem])))
    pd.DataFrame(rows).to_csv(OUT / "sets.tsv", sep="\t", index=False)
    return genes


def species_order(s82) -> list[str]:
    s81 = s82.load_s81()
    return list(s81.SPECIES_ORDER)


def sequences_by_gene(s82, sp: str) -> tuple[dict[str, set[str]], set[str]]:
    recs = s82.fasta_records(STORY / f"{sp}_proximal.fa")
    by_gene: dict[str, set[str]] = {}
    for sid, gene, _raw in recs:
        by_gene.setdefault(gene, set()).add(sid)
    return by_gene, {sid for sid, _, _ in recs}


def fisher_catalog(s82, genes: dict[str, set[str]], catalog: str, hit_path, tf_of) -> None:
    """tf_of(dataframe) -> Series of TF symbols aligned to hits."""
    fisher_block = s82.load_s81().fisher_block
    sps = species_order(s82)
    rows = []
    for sp in sps:
        by_gene, all_sid = sequences_by_gene(s82, sp)
        hits = pd.read_csv(hit_path(sp), sep="\t", comment="#")
        hits = hits[hits["sequence_name"].notna()]
        tf = tf_of(hits)
        by_tf: dict[str, set[str]] = {}
        for name, seq in zip(tf, hits["sequence_name"].astype(str)):
            if not isinstance(name, str) or not name:
                continue
            by_tf.setdefault(name, set()).add(seq)
        print(f"  {catalog} {sp} seq {len(all_sid)} motifs {len(by_tf)}", flush=True)
        for stem, title in SETS:
            pos = set().union(*(by_gene.get(g, set()) for g in genes[stem]))
            neg = all_sid - pos
            if len(pos) < 2 or len(neg) < 20:
                print(f"    {stem} пропуск pos={len(pos)}", flush=True)
                continue
            for name, hit_seq in by_tf.items():
                rec = fisher_block(hit_seq, pos, neg)
                rec.update(set=stem, title=title, species=sp, tf=name,
                           n_pos=len(pos), n_neg=len(neg))
                rows.append(rec)
    fish = s82.bh_within(pd.DataFrame(rows))
    dest = OUT / catalog
    dest.mkdir(parents=True, exist_ok=True)
    fish.to_csv(dest / "motif_fisher_proximal.tsv", sep="\t", index=False)
    ge3 = s82.ge3_table(fish)
    ge3.to_csv(dest / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
    n_bh = int((fish.bh < 0.05).sum())
    print(f"{catalog} BH<0.05 {n_bh}  ge3 {len(ge3)}", flush=True)
    if len(ge3):
        print(ge3.to_string(index=False), flush=True)


def jaspar_tf(hits: pd.DataFrame) -> pd.Series:
    return hits["motif_alt_id"].astype(str).str.upper()


def hoco_tf(hits: pd.DataFrame) -> pd.Series:
    symbol = pd.read_csv(HOCO / "models_kept.tsv", sep="\t").set_index("motif_id")["tf"]
    return hits["motif_id"].map(symbol)


def run_fisher(s82, genes) -> None:
    print("JASPAR FIMO, девять наборов", flush=True)
    fisher_catalog(
        s82, genes, "fimo_jaspar",
        lambda sp: STORY / f"{sp}_fimo_proximal.tsv",
        jaspar_tf,
    )
    print("HOCOMOCO FIMO, девять наборов", flush=True)
    fisher_catalog(
        s82, genes, "hocomoco",
        lambda sp: HOCO / f"{sp}_fimo.tsv",
        hoco_tf,
    )


def run_sea(s82, genes) -> None:
    s81 = s82.load_s81()
    sea_dir = OUT / "sea"
    jobs = []
    for sp in s81.SPECIES_ORDER:
        recs = s82.fasta_records(STORY / f"{sp}_proximal.fa")
        for stem, title in SETS:
            pos = sea_dir / "fasta" / f"{sp}_{stem}_pos.fa"
            neg = sea_dir / "fasta" / f"{sp}_{stem}_neg.fa"
            n_pos, n_neg = s82.split_set_fasta(recs, genes[stem], pos, neg)
            h, note = s82.sea_hofract(n_pos)
            jobs.append(dict(
                sp=sp, stem=stem, title=title, pos=pos, neg=neg,
                n_pos=n_pos, n_neg=n_neg, hofract=h, note=note,
                outdir=sea_dir / "runs" / sp / stem,
            ))

    def one(job):
        tsv = job["outdir"] / "sea.tsv"
        err_path = job["outdir"] / "sea.err"
        if tsv.exists() and tsv.stat().st_size > 50 and (job["outdir"] / "sea.ok").exists():
            return job["sp"], job["stem"], "cached"
        if job["n_pos"] < 2 or job["n_neg"] < 20:
            return job["sp"], job["stem"], "skip"
        job["outdir"].mkdir(parents=True, exist_ok=True)
        cmd = [
            str(SEA_BIN), "--p", str(job["pos"]), "--n", str(job["neg"]),
            "--m", str(JASPAR_MEME), "--text", "--noseqs", "--notrim",
            "--hofract", f"{job['hofract']:.6f}", "--seed", "1",
            "--pvalue", "--thresh", "1", "--verbosity", "1",
        ]
        with tsv.open("w") as out, err_path.open("w") as err:
            proc = subprocess.run(cmd, stdout=out, stderr=err)
        if proc.returncode != 0 or tsv.stat().st_size < 200:
            raise SystemExit(f"SEA {job['sp']} {job['stem']} rc={proc.returncode}")
        (job["outdir"] / "sea.ok").write_text(job["note"] + "\n")
        return job["sp"], job["stem"], "ran"

    print(f"SEA, заданий {len(jobs)}", flush=True)
    with ThreadPoolExecutor(max_workers=12) as pool:
        futs = [pool.submit(one, job) for job in jobs]
        for fut in as_completed(futs):
            sp, stem, how = fut.result()
            print(f"  {sp} {stem} {how}", flush=True)
    # collect like script 82: SEA table columns vary; reuse its parser by reading sea.tsv
    fish_rows = []
    for job in jobs:
        tsv = job["outdir"] / "sea.tsv"
        if not tsv.exists() or tsv.stat().st_size < 50:
            continue
        tab = pd.read_csv(tsv, sep="\t")
        if "ALT_ID" not in tab.columns or "PVALUE" not in tab.columns:
            raise SystemExit(f"непонятный sea.tsv {tsv}")
        tab = tab.dropna(subset=["ALT_ID", "ID", "PVALUE"])
        tab["tf"] = tab["ALT_ID"].astype(str).str.upper()
        if tab["tf"].duplicated().any():
            raise SystemExit(f"SEA {job['sp']} {job['stem']}: повтор символа")
        for rec in tab.itertuples(index=False):
            fish_rows.append(dict(
                set=job["stem"], title=job["title"], species=job["sp"],
                tf=rec.tf, p=float(rec.PVALUE),
                n_pos=job["n_pos"], n_neg=job["n_neg"],
            ))
    if not fish_rows:
        print("SEA: таблица p не собрана, сырые sea.tsv лежат в runs/", flush=True)
        return
    fish = s82.bh_within(pd.DataFrame(fish_rows))
    fish.to_csv(sea_dir / "motif_fisher_proximal.tsv", sep="\t", index=False)
    ge3 = s82.ge3_table(fish)
    ge3.to_csv(sea_dir / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
    print(f"SEA BH<0.05 {int((fish.bh < 0.05).sum())}  ge3 {len(ge3)}", flush=True)


def run_monalisa(s82) -> None:
    s81 = s82.load_s81()
    mona = OUT / "monalisa"
    mona.mkdir(parents=True, exist_ok=True)

    def one(sp):
        out = mona / f"{sp}_binned.tsv"
        log = mona / f"{sp}.log"
        if (mona / f"{sp}.ok").exists() and out.exists() and out.stat().st_size > 50:
            return sp, "cached"
        cmd = [
            "Rscript", "--vanilla", str(R_SCRIPT),
            "--meme", str(JASPAR_MEME),
            "--fasta", str(STORY / f"{sp}_proximal.fa"),
            "--species", sp,
            "--geneset-dir", str(OUT / "genes"),
            "--sets", str(OUT / "sets.tsv"),
            "--out", str(out),
            "--cores", "16",
        ]
        with log.open("w") as fh:
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            raise SystemExit(f"monaLisa {sp} rc={proc.returncode} см. {log}")
        (mona / f"{sp}.ok").write_text("ok\n")
        return sp, "ran"

    print("monaLisa, семь видов", flush=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        futs = [pool.submit(one, sp) for sp in s81.SPECIES_ORDER]
        for fut in as_completed(futs):
            sp, how = fut.result()
            print(f"  {sp} {how}", flush=True)
    parts = [pd.read_csv(mona / f"{sp}_binned.tsv", sep="\t") for sp in s81.SPECIES_ORDER]
    fish = s82.bh_within(pd.concat(parts, ignore_index=True))
    fish.to_csv(mona / "motif_fisher_proximal.tsv", sep="\t", index=False)
    ge3 = s82.ge3_table(fish)
    ge3.to_csv(mona / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
    print(f"monaLisa BH<0.05 {int((fish.bh < 0.05).sum())}  ge3 {len(ge3)}", flush=True)


def main() -> None:
    stage = sys.argv[1] if len(sys.argv) > 1 else "fisher"
    s82 = load_s82()
    genes = load_sets()
    print("наборы", {k: len(v) for k, v in genes.items()}, flush=True)
    if stage in ("fisher", "all"):
        run_fisher(s82, genes)
    if stage in ("sea", "all"):
        run_sea(s82, genes)
    if stage in ("monalisa", "all"):
        run_monalisa(s82)


if __name__ == "__main__":
    main()
