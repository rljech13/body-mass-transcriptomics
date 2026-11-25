#!/usr/bin/env python3
"""Three lightweight motif tests next to the frozen JASPAR FIMO.

Does not rewrite story_cre_motifs/*_fimo_proximal.tsv or the frozen tables.

  hocomoco  FIMO 5.3.0, threshold 1e-4, two-sided Fisher, BH within
            set x species. Catalog is HOCOMOCO v12 H12CORE, one model
            per factor (index 0). Names are HGNC symbols (ERR1 -> ESRRA).
  sea       SEA from MEME 5.5 against the other proximal promoters.
            Score threshold is chosen on a hold-out of at least 6
            sequences when the set is large enough for that.
  monalisa  monaLisa::calcBinnedMotifEnrR, bins set vs rest, weights
            by GC and 3-mers. Calls the R script.
  compare   Gene-symbol lists. No averaged score.
"""
from __future__ import annotations

import importlib.util
import math
import re
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
CMP = STORY / "method_compare"
FROZEN = CMP / "fimo_jaspar_frozen"
HOCO = CMP / "hocomoco"
SEA_DIR = CMP / "sea"
MONA = CMP / "monalisa"
REPORT = CMP / "compare"
FIMO = Path("/path/to/user-home/miniconda3/envs/meme/bin/fimo")
SEA_BIN = CMP / "meme55" / "bin" / "sea"
R_SCRIPT = SCRIPTS / "82_monalisa_binned.R"
FIMO_THRESH = "1e-4"
MOTIF_RE = re.compile(r"^(.+)\.H12CORE\.(\d+)\.([^.]+)\.([A-D])$")
JASPAR_MOTIF_RE = re.compile(r"^MOTIF\s+(\S+)(?:\s+(\S+))?")


def load_s81():
    sys.path.insert(0, str(SCRIPTS))
    spec = importlib.util.spec_from_file_location(
        "story81", SCRIPTS / "81_story_cre_motifs.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def split_meme(text: str) -> tuple[list[str], list[tuple[str, list[str]]]]:
    header: list[str] = []
    blocks: list[tuple[str, list[str]]] = []
    cur_name = None
    cur: list[str] = []
    for line in text.splitlines(keepends=True):
        if line.startswith("MOTIF "):
            if cur_name is not None:
                blocks.append((cur_name, cur))
            cur_name = line.split()[1]
            cur = [line]
        elif cur_name is None:
            header.append(line)
        else:
            cur.append(line)
    if cur_name is not None:
        blocks.append((cur_name, cur))
    return header, blocks


def write_blocks(path: Path, header: list[str], blocks: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.writelines(header)
        for block in blocks:
            fh.writelines(block)
            if not block[-1].endswith("\n"):
                fh.write("\n")


def hocomoco_symbol_table() -> pd.DataFrame:
    master = pd.read_csv(HOCO / "hocomoco_v12_masterlist.tsv", sep="\t")
    mapping = pd.read_csv(HOCO / "hocomoco_benchmarking_ID_mapping.tsv", sep="\t")
    human = master[master["curated:uniprot_id"].astype(str).str.endswith("_HUMAN")]
    u2g = dict(zip(human["curated:uniprot_id"], human["auto:gene_symbol"]))
    mapping = mapping.copy()
    mapping["prefix"] = mapping["HOCOMOCO_ID"].str.split(".").str[0]
    rows = []
    seen = set()
    for rec in mapping.itertuples(index=False):
        if rec.prefix in seen:
            continue
        if ".H12CORE.0." not in rec.HOCOMOCO_ID:
            continue
        seen.add(rec.prefix)
        gene = u2g.get(rec.Uniprot_HUMAN)
        if gene is None or (isinstance(gene, float) and math.isnan(gene)):
            raise SystemExit(f"no HGNC for {rec.prefix} / {rec.Uniprot_HUMAN}")
        rows.append(dict(prefix=rec.prefix, uniprot=rec.Uniprot_HUMAN, tf=str(gene).upper()))
    tab = pd.DataFrame(rows)
    if tab["tf"].duplicated().any():
        bad = tab[tab["tf"].duplicated(keep=False)]
        raise SystemExit(f"two HOCOMOCO factors map to one symbol\n{bad}")
    return tab


def build_hocomoco_meme() -> Path:
    src = HOCO / "H12CORE_meme_format.meme"
    dst = HOCO / "H12CORE_one_per_tf.meme"
    symbols = hocomoco_symbol_table().set_index("prefix")
    header, blocks = split_meme(src.read_text())
    kept = []
    meta = []
    for name, block in blocks:
        m = MOTIF_RE.match(name)
        if m is None or int(m.group(2)) != 0:
            continue
        prefix = m.group(1)
        if prefix not in symbols.index:
            raise SystemExit(f"motif {name} has no HGNC map")
        kept.append(block)
        meta.append(dict(
            motif_id=name,
            prefix=prefix,
            origin=m.group(3),
            grade=m.group(4),
            tf=symbols.loc[prefix, "tf"],
            uniprot=symbols.loc[prefix, "uniprot"],
        ))
    write_blocks(dst, header, kept)
    pd.DataFrame(meta).to_csv(HOCO / "models_kept.tsv", sep="\t", index=False)
    print(f"HOCOMOCO одна модель на фактор: {len(kept)}  -> {dst}", flush=True)
    return dst


def jaspar_nsites(block: list[str]) -> int:
    for line in block:
        m = re.search(r"nsites=\s*([0-9]+)", line)
        if m:
            return int(m.group(1))
    return 0


def jaspar_version(motif_id: str) -> int:
    m = re.search(r"\.(\d+)$", motif_id)
    return int(m.group(1)) if m else 0


def build_jaspar_one_per_symbol(meme: Path) -> Path:
    """One JASPAR matrix per uppercase alt id.

    The frozen FIMO ORs every version that shares an alt id. Versions
    that share a name are 53 alt ids; the OXPHOS six and the watch list
    each have a single matrix. Here the kept matrix is the one with the
    most sites, then the higher JASPAR version.
    """
    dst = CMP / "JASPAR2024_one_per_symbol.meme"
    header, blocks = split_meme(meme.read_text())
    best: dict[str, tuple[int, int, str, list[str]]] = {}
    for name, block in blocks:
        m = JASPAR_MOTIF_RE.match(block[0].rstrip("\n"))
        if m is None:
            continue
        motif_id, alt = m.group(1), (m.group(2) or m.group(1))
        symbol = alt.upper()
        key = (jaspar_nsites(block), jaspar_version(motif_id))
        prev = best.get(symbol)
        if prev is None or key > (prev[0], prev[1]):
            best[symbol] = (key[0], key[1], motif_id, block)
    chosen = list(best.values())
    write_blocks(dst, header, [c[3] for c in chosen])
    pd.DataFrame(
        [dict(tf=sym, motif_id=rec[2], nsites=rec[0]) for sym, rec in best.items()]
    ).to_csv(CMP / "jaspar_one_per_symbol.tsv", sep="\t", index=False)
    print(f"JASPAR одна матрица на имя: {len(chosen)} из {len(blocks)}", flush=True)
    return dst


def fasta_records(path: Path) -> list[tuple[str, str, str]]:
    """Return (header_without_gt, gene, raw_record_including_header)."""
    out = []
    header = None
    buf: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    out.append((header, header.split("|", 1)[0], "".join(buf)))
                header = line[1:].strip()
                buf = [line]
            else:
                buf.append(line)
    if header is not None:
        out.append((header, header.split("|", 1)[0], "".join(buf)))
    return out


def bh_within(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["bh"] = np.nan
    for (_, _), idx in df.groupby(["set", "species"]).groups.items():
        p = df.loc[idx, "p"].to_numpy(float)
        ok = np.isfinite(p)
        if ok.sum() == 0:
            continue
        adj = np.full(len(p), np.nan)
        adj[ok] = stats.false_discovery_control(p[ok], method="bh")
        df.loc[idx, "bh"] = adj
    return df


def ge3_table(df: pd.DataFrame) -> pd.DataFrame:
    sig = df[df["bh"] < 0.05]
    rows = []
    if sig.empty:
        return pd.DataFrame(columns=["set", "tf", "n_species"])
    for (stem, tf), g in sig.groupby(["set", "tf"]):
        k = int(g["species"].nunique())
        if k >= 3:
            rows.append(dict(set=stem, tf=tf, n_species=k))
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["set", "tf", "n_species"])
    return out.sort_values(["set", "n_species", "tf"], ascending=[True, False, True])


def watch_vs_mass(watch: pd.DataFrame, species_order: list[str], logm: np.ndarray) -> pd.DataFrame:
    rows = []
    for (stem, tf), g in watch.groupby(["set", "tf"]):
        g = g.set_index("species").reindex(species_order)
        p = g["p"].to_numpy(float)
        y = -np.log10(np.clip(p, 1e-300, 1))
        ok = np.isfinite(y)
        rec = dict(set=stem, tf=tf, n_points=int(ok.sum()), r=np.nan, p_mass=np.nan)
        if ok.sum() >= 5:
            fit = stats.linregress(logm[ok], y[ok])
            rec["r"] = float(fit.rvalue)
            rec["p_mass"] = float(fit.pvalue)
        bh = g["bh"] if "bh" in g.columns else pd.Series(np.nan, index=g.index)
        rec["n_species_bh"] = int((bh < 0.05).sum())
        for sp in species_order:
            rec[f"p_{sp.lower()}"] = float(g.loc[sp, "p"]) if sp in g.index and np.isfinite(g.loc[sp, "p"]) else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def attach_watch_bh(watch: pd.DataFrame, fisher: pd.DataFrame) -> pd.DataFrame:
    key = fisher.set_index(["set", "species", "tf"])["bh"]
    watch = watch.copy()
    watch["bh"] = [
        key.get((r.set, r.species, r.tf), np.nan) for r in watch.itertuples(index=False)
    ]
    return watch


def run_hocomoco(s81) -> None:
    meme = build_hocomoco_meme()
    symbol = pd.read_csv(HOCO / "models_kept.tsv", sep="\t").set_index("motif_id")["tf"].to_dict()
    HOCO.mkdir(parents=True, exist_ok=True)
    jobs = []
    for sp in s81.SPECIES_ORDER:
        fasta = STORY / f"{sp}_proximal.fa"
        if not fasta.exists():
            raise SystemExit(f"нет fasta {fasta}")
        out = HOCO / f"{sp}_fimo.tsv"
        log = HOCO / f"{sp}_fimo.log"
        jobs.append((sp, fasta, meme, out, log))

    def one(job):
        sp, fasta, meme_path, out, log = job
        if out.exists() and out.stat().st_size > 1000 and log.exists() and "fimo failed" not in log.read_text(errors="replace"):
            # a finished fimo --text file starts with a header; require the log to be closed by us
            if (HOCO / f"{sp}_fimo.ok").exists():
                return sp, "cached", out.stat().st_size
        cmd = [
            str(FIMO), "--text", "--no-qvalue", "--verbosity", "1",
            "--thresh", FIMO_THRESH, str(meme_path), str(fasta),
        ]
        with out.open("w") as fh, log.open("w") as err:
            proc = subprocess.run(cmd, stdout=fh, stderr=err)
        if proc.returncode != 0:
            log.write_text(log.read_text() + f"\nfimo failed {proc.returncode}\n")
            raise SystemExit(f"FIMO {sp} rc={proc.returncode}")
        (HOCO / f"{sp}_fimo.ok").write_text("ok\n")
        return sp, "ran", out.stat().st_size

    print("FIMO HOCOMOCO, семь видов", flush=True)
    with ThreadPoolExecutor(max_workers=7) as pool:
        futs = [pool.submit(one, job) for job in jobs]
        for fut in as_completed(futs):
            sp, how, size = fut.result()
            print(f"  {sp} {how} {size}", flush=True)

    story = {stem: set(s81.genes_of(stem)) for stem, _ in s81.SETS}
    fish_rows = []
    watch_rows = []
    for sp in s81.SPECIES_ORDER:
        recs = fasta_records(STORY / f"{sp}_proximal.fa")
        by_gene: dict[str, set[str]] = {}
        for sid, gene, _ in recs:
            by_gene.setdefault(gene, set()).add(sid)
        all_sid = {sid for sid, _, _ in recs}
        hits = pd.read_csv(
            HOCO / f"{sp}_fimo.tsv", sep="\t", comment="#",
            usecols=lambda c: c in {"motif_id", "motif_alt_id", "sequence_name"},
        )
        if hits.empty:
            print(f"  {sp} нет хитов", flush=True)
            continue
        id_col = "motif_id"
        hits["tf"] = hits[id_col].map(symbol)
        missing = hits["tf"].isna().sum()
        if missing:
            raise SystemExit(f"{sp}: {missing} хитов без символа")
        by_tf: dict[str, set[str]] = {}
        for tf, seq in zip(hits["tf"], hits["sequence_name"].astype(str)):
            by_tf.setdefault(str(tf), set()).add(seq)
        print(f"  {sp} seq {len(all_sid)} motifs {len(by_tf)}", flush=True)
        for stem, title in s81.SETS:
            pos = set().union(*(by_gene.get(g, set()) for g in story[stem])) if story[stem] else set()
            neg = all_sid - pos
            if len(pos) < 2 or len(neg) < 20:
                print(f"    {stem} пропуск pos={len(pos)} neg={len(neg)}", flush=True)
                continue
            for tf, hit_seq in by_tf.items():
                rec = s81.fisher_block(hit_seq, pos, neg)
                rec.update(set=stem, title=title, species=sp, tf=tf, n_pos=len(pos), n_neg=len(neg))
                fish_rows.append(rec)
            for w in s81.WATCH:
                rec = s81.fisher_block(by_tf.get(w, set()), pos, neg)
                rec.update(set=stem, title=title, species=sp, tf=w, n_pos=len(pos), n_neg=len(neg))
                watch_rows.append(rec)
    fish = bh_within(pd.DataFrame(fish_rows))
    fish.to_csv(HOCO / "motif_fisher_proximal.tsv", sep="\t", index=False)
    ge3 = ge3_table(fish)
    ge3.to_csv(HOCO / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
    watch = attach_watch_bh(pd.DataFrame(watch_rows), fish)
    watch.to_csv(HOCO / "motif_watch_fisher.tsv", sep="\t", index=False)
    logm = np.log10([s81.SPECIES[sp]["mass_g"] for sp in s81.SPECIES_ORDER])
    watch_vs_mass(watch, list(s81.SPECIES_ORDER), logm).to_csv(
        HOCO / "watch_tf_enrichment_vs_mass.tsv", sep="\t", index=False
    )
    print("HOCOMOCO BH<0.05", int((fish.bh < 0.05).sum()), "ge3", ge3.to_dict("records"), flush=True)


def sea_hofract(n_pos: int) -> tuple[float, str]:
    """SEA drops the hold-out when it would contain fewer than 5 sequences.

    Default 0.1 is kept when it already yields at least 6. Otherwise the
    fraction is raised so the hold-out has 6, and only if at least 5
    sequences remain to be scored. Smaller sets cannot host that hold-out.
    """
    if n_pos < 11:
        return 0.1, "too_few_for_holdout"
    if 0.1 * n_pos >= 6:
        return 0.1, "default_0.1"
    return 6 / n_pos, "raised_holdout_6"


def split_set_fasta(recs, genes: set[str], pos_path: Path, neg_path: Path) -> tuple[int, int]:
    pos_path.parent.mkdir(parents=True, exist_ok=True)
    n_pos = n_neg = 0
    with pos_path.open("w") as fp, neg_path.open("w") as fn:
        for _sid, gene, raw in recs:
            if gene in genes:
                fp.write(raw)
                n_pos += 1
            else:
                fn.write(raw)
                n_neg += 1
    return n_pos, n_neg


def run_sea(s81) -> None:
    meme = build_jaspar_one_per_symbol(s81.MEME)
    SEA_DIR.mkdir(parents=True, exist_ok=True)
    story = {stem: set(s81.genes_of(stem)) for stem, _ in s81.SETS}
    jobs = []
    for sp in s81.SPECIES_ORDER:
        recs = fasta_records(STORY / f"{sp}_proximal.fa")
        for stem, title in s81.SETS:
            pos = SEA_DIR / "fasta" / f"{sp}_{stem}_pos.fa"
            neg = SEA_DIR / "fasta" / f"{sp}_{stem}_neg.fa"
            n_pos, n_neg = split_set_fasta(recs, story[stem], pos, neg)
            h, note = sea_hofract(n_pos)
            outdir = SEA_DIR / "runs" / sp / stem
            jobs.append(dict(
                sp=sp, stem=stem, title=title, pos=pos, neg=neg,
                n_pos=n_pos, n_neg=n_neg, hofract=h, note=note, outdir=outdir,
            ))

    def one(job):
        tsv = job["outdir"] / "sea.tsv"
        err_path = job["outdir"] / "sea.err"
        if tsv.exists() and tsv.stat().st_size > 50 and (job["outdir"] / "sea.ok").exists():
            return job["sp"], job["stem"], "cached", err_path.read_text(errors="replace") if err_path.exists() else ""
        if job["n_pos"] < 2 or job["n_neg"] < 20:
            return job["sp"], job["stem"], "skip", ""
        job["outdir"].mkdir(parents=True, exist_ok=True)
        # --text skips the HTML template, which aborts with
        # "Template does not contain data section" on this MEME build.
        cmd = [
            str(SEA_BIN),
            "--p", str(job["pos"]),
            "--n", str(job["neg"]),
            "--m", str(meme),
            "--text",
            "--noseqs", "--notrim",
            "--hofract", f"{job['hofract']:.6f}",
            "--seed", "1",
            "--pvalue", "--thresh", "1",
            "--verbosity", "1",
        ]
        with tsv.open("w") as out, err_path.open("w") as err:
            proc = subprocess.run(cmd, stdout=out, stderr=err)
        if proc.returncode != 0 or tsv.stat().st_size < 200:
            raise SystemExit(f"SEA {job['sp']} {job['stem']} rc={proc.returncode}")
        (job["outdir"] / "sea.ok").write_text(job["note"] + "\n")
        return job["sp"], job["stem"], "ran", err_path.read_text(errors="replace")

    print(f"SEA, заданий {len(jobs)}", flush=True)
    notes = []
    with ThreadPoolExecutor(max_workers=12) as pool:
        futs = [pool.submit(one, job) for job in jobs]
        done = { (j["sp"], j["stem"]): j for j in jobs }
        for fut in as_completed(futs):
            sp, stem, how, err = fut.result()
            job = done[(sp, stem)]
            holdout_used = "No hold-out set was created" not in err
            notes.append(dict(
                species=sp, set=stem, n_pos=job["n_pos"], n_neg=job["n_neg"],
                hofract=job["hofract"], plan=job["note"], status=how,
                holdout_used=holdout_used and how != "skip",
                length_warning="different length distributions" in err,
            ))
            print(f"  {sp} {stem} {how} holdout={holdout_used}", flush=True)
    pd.DataFrame(notes).to_csv(SEA_DIR / "run_notes.tsv", sep="\t", index=False)

    fish_rows = []
    for job in jobs:
        tsv = job["outdir"] / "sea.tsv"
        if not tsv.exists():
            continue
        tab = pd.read_csv(tsv, sep="\t")
        if "ALT_ID" not in tab.columns or "PVALUE" not in tab.columns:
            raise SystemExit(f"непонятный sea.tsv {tsv}: {tab.columns.tolist()}")
        tab = tab.dropna(subset=["ALT_ID", "ID", "PVALUE"])
        tab["tf"] = tab["ALT_ID"].astype(str).str.upper()
        # The reduced MEME has one matrix per symbol. A duplicate name would
        # make the smallest p a hidden winner; stop instead of collapsing.
        if tab["tf"].duplicated().any():
            raise SystemExit(f"SEA {job['sp']} {job['stem']}: повтор символа")
        for rec in tab.itertuples(index=False):
            fish_rows.append(dict(
                set=job["stem"], title=job["title"], species=job["sp"], tf=rec.tf,
                motif_id=getattr(rec, "ID"),
                p=float(rec.PVALUE),
                enr_ratio=float(rec.ENR_RATIO) if np.isfinite(rec.ENR_RATIO) else np.nan,
                score_thr=float(rec.SCORE_THR) if np.isfinite(getattr(rec, "SCORE_THR", np.nan)) else np.nan,
                tp=float(rec.TP), fp=float(rec.FP),
                n_pos=job["n_pos"], n_neg=job["n_neg"],
                sea_q=float(rec.QVALUE) if np.isfinite(rec.QVALUE) else np.nan,
            ))
    fish = bh_within(pd.DataFrame(fish_rows))
    fish.to_csv(SEA_DIR / "motif_fisher_proximal.tsv", sep="\t", index=False)
    ge3 = ge3_table(fish)
    ge3.to_csv(SEA_DIR / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
    # watch: p=1 when the name was not reported after a finished run
    watch_rows = []
    have = fish.set_index(["set", "species", "tf"])
    for job in jobs:
        if not (job["outdir"] / "sea.tsv").exists():
            continue
        for w in s81.WATCH:
            key = (job["stem"], job["sp"], w)
            if key in have.index:
                row = have.loc[key]
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                p, bh = float(row["p"]), float(row["bh"])
            else:
                p, bh = 1.0, np.nan
            watch_rows.append(dict(set=job["stem"], species=job["sp"], tf=w, p=p, bh=bh))
    watch = pd.DataFrame(watch_rows)
    watch.to_csv(SEA_DIR / "motif_watch_fisher.tsv", sep="\t", index=False)
    logm = np.log10([s81.SPECIES[sp]["mass_g"] for sp in s81.SPECIES_ORDER])
    watch_vs_mass(watch, list(s81.SPECIES_ORDER), logm).to_csv(
        SEA_DIR / "watch_tf_enrichment_vs_mass.tsv", sep="\t", index=False
    )
    print("SEA BH<0.05", int((fish.bh < 0.05).sum()) if len(fish) else 0, "ge3", ge3.to_dict("records"), flush=True)


def run_monalisa(s81) -> None:
    meme = build_jaspar_one_per_symbol(s81.MEME)
    MONA.mkdir(parents=True, exist_ok=True)
    geneset_dir = s81.DATA / "genesets"

    def one(sp):
        out = MONA / f"{sp}_binned.tsv"
        log = MONA / f"{sp}.log"
        if (MONA / f"{sp}.ok").exists() and out.exists() and out.stat().st_size > 50:
            return sp, "cached"
        cmd = [
            "Rscript", "--vanilla", str(R_SCRIPT),
            "--meme", str(meme),
            "--fasta", str(STORY / f"{sp}_proximal.fa"),
            "--species", sp,
            "--geneset-dir", str(geneset_dir),
            "--out", str(out),
            "--cores", "12",
        ]
        with log.open("w") as fh:
            proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            raise SystemExit(f"monaLisa {sp} rc={proc.returncode} см. {log}")
        (MONA / f"{sp}.ok").write_text("ok\n")
        return sp, "ran"

    print("monaLisa, семь видов", flush=True)
    with ThreadPoolExecutor(max_workers=7) as pool:
        futs = [pool.submit(one, sp) for sp in s81.SPECIES_ORDER]
        for fut in as_completed(futs):
            sp, how = fut.result()
            print(f"  {sp} {how}", flush=True)

    parts = [pd.read_csv(MONA / f"{sp}_binned.tsv", sep="\t") for sp in s81.SPECIES_ORDER]
    fish = bh_within(pd.concat(parts, ignore_index=True))
    fish.to_csv(MONA / "motif_fisher_proximal.tsv", sep="\t", index=False)
    ge3 = ge3_table(fish)
    ge3.to_csv(MONA / "motifs_bh_in_ge3_species.tsv", sep="\t", index=False)
    have = fish.set_index(["set", "species", "tf"])
    watch_rows = []
    for sp in s81.SPECIES_ORDER:
        for stem, _title in s81.SETS:
            for w in s81.WATCH:
                key = (stem, sp, w)
                if key in have.index:
                    row = have.loc[key]
                    if isinstance(row, pd.DataFrame):
                        row = row.iloc[0]
                    p, bh = float(row["p"]), float(row["bh"])
                else:
                    p, bh = 1.0, np.nan
                watch_rows.append(dict(set=stem, species=sp, tf=w, p=p, bh=bh))
    watch = pd.DataFrame(watch_rows)
    watch.to_csv(MONA / "motif_watch_fisher.tsv", sep="\t", index=False)
    logm = np.log10([s81.SPECIES[sp]["mass_g"] for sp in s81.SPECIES_ORDER])
    watch_vs_mass(watch, list(s81.SPECIES_ORDER), logm).to_csv(
        MONA / "watch_tf_enrichment_vs_mass.tsv", sep="\t", index=False
    )
    print("monaLisa BH<0.05", int((fish.bh < 0.05).sum()), "ge3", ge3.to_dict("records"), flush=True)


def _ge3_sets(path: Path) -> dict[str, set[str]]:
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t")
    if df.empty or "tf" not in df.columns:
        return {}
    out: dict[str, set[str]] = {}
    for stem, g in df.groupby("set"):
        out[str(stem)] = set(g["tf"].astype(str))
    return out


def uniform_watch(method: str, path: Path) -> pd.DataFrame:
    w = pd.read_csv(path, sep="\t")
    out = w[["set", "tf", "n_points", "n_species_bh", "r", "p_mass"]].copy()
    out.insert(0, "method", method)
    return out


def fimo_watch_uniform() -> pd.DataFrame:
    """Frozen watch file stores the mass-correlation p in column `p` and has no BH count."""
    w = pd.read_csv(FROZEN / "watch_tf_enrichment_vs_mass.tsv", sep="\t")
    fish = pd.read_csv(FROZEN / "motif_fisher_proximal.tsv", sep="\t")
    counts = (
        fish.loc[fish["bh"] < 0.05]
        .groupby(["set", "tf"])["species"]
        .nunique()
    )
    rows = []
    for rec in w.itertuples(index=False):
        rows.append(dict(
            method="fimo_jaspar",
            set=rec.set,
            tf=rec.tf,
            n_points=7,
            n_species_bh=int(counts.get((rec.set, rec.tf), 0)),
            r=float(rec.r),
            p_mass=float(rec.p),
        ))
    return pd.DataFrame(rows)


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def run_compare(s81) -> None:
    methods = {
        "fimo_jaspar": FROZEN / "motifs_bh_in_ge3_species.tsv",
        "hocomoco": HOCO / "motifs_bh_in_ge3_species.tsv",
        "monalisa": MONA / "motifs_bh_in_ge3_species.tsv",
        "sea": SEA_DIR / "motifs_bh_in_ge3_species.tsv",
    }
    for name, path in methods.items():
        if not path.exists():
            raise SystemExit(f"нет таблицы {name}: {path}")
    # frozen reference must still match the live FIMO table
    live = (STORY / "motifs_bh_in_ge3_species.tsv").read_bytes()
    frozen = (FROZEN / "motifs_bh_in_ge3_species.tsv").read_bytes()
    if live != frozen:
        raise SystemExit("живая таблица FIMO разъехалась с заморозкой")
    lists = {name: _ge3_sets(path) for name, path in methods.items()}
    stems = [stem for stem, _ in s81.SETS]
    jac_rows = []
    ret_rows = []
    long_rows = []
    for stem in stems:
        base = lists["fimo_jaspar"].get(stem, set())
        for name, grouped in lists.items():
            got = grouped.get(stem, set())
            for tf in sorted(got):
                src = pd.read_csv(methods[name], sep="\t")
                n = int(src[(src["set"] == stem) & (src["tf"] == tf)]["n_species"].iloc[0])
                long_rows.append(dict(set=stem, method=name, tf=tf, n_species=n))
            if name == "fimo_jaspar":
                continue
            jac_rows.append(dict(
                set=stem, method=name,
                n_fimo=len(base), n_method=len(got), n_both=len(base & got),
                jaccard=jaccard(base, got),
                only_fimo=";".join(sorted(base - got)),
                only_method=";".join(sorted(got - base)),
            ))
        for tf in sorted(base):
            rec = dict(set=stem, tf=tf)
            for name in ("hocomoco", "monalisa", "sea"):
                src = pd.read_csv(methods[name], sep="\t")
                hit = src[(src["set"] == stem) & (src["tf"] == tf)]
                rec[f"n_{name}"] = int(hit["n_species"].iloc[0]) if len(hit) else 0
                rec[f"kept_{name}"] = bool(len(hit))
            ret_rows.append(rec)
    REPORT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(long_rows).to_csv(REPORT / "ge3_by_method.tsv", sep="\t", index=False)
    pd.DataFrame(jac_rows).to_csv(REPORT / "jaccard_vs_fimo.tsv", sep="\t", index=False)
    pd.DataFrame(ret_rows).to_csv(REPORT / "fimo_ge3_retention.tsv", sep="\t", index=False)

    watch = pd.concat(
        [
            fimo_watch_uniform(),
            uniform_watch("hocomoco", HOCO / "watch_tf_enrichment_vs_mass.tsv"),
            uniform_watch("monalisa", MONA / "watch_tf_enrichment_vs_mass.tsv"),
            uniform_watch("sea", SEA_DIR / "watch_tf_enrichment_vs_mass.tsv"),
        ],
        ignore_index=True,
    )
    watch.to_csv(REPORT / "watch_vs_mass.tsv", sep="\t", index=False)
    pd.DataFrame([
        dict(method="fimo_jaspar", catalog="JASPAR2024 CORE vertebrates",
             scan="FIMO 5.3.0, thresh 1e-4", test="Fisher two-sided",
             background="other proximal peaks", composition="none",
             matrices="879, versions of one name OR-ed"),
        dict(method="hocomoco", catalog="HOCOMOCO v12 H12CORE, index 0, HGNC",
             scan="FIMO 5.3.0, thresh 1e-4", test="Fisher two-sided",
             background="other proximal peaks", composition="none",
             matrices="949, one model per symbol"),
        dict(method="monalisa", catalog="JASPAR2024, one matrix per symbol",
             scan="monaLisa 1.18.0, min.score 10", test="Fisher one-sided enrichment",
             background="other bin", composition="GC bins and k-mers to 3",
             matrices="816"),
        dict(method="sea", catalog="JASPAR2024, one matrix per symbol",
             scan="SEA 5.5.9, --pvalue --thresh 1, --notrim",
             test="SEA score, BH recomputed inside set x species",
             background="other proximal peaks",
             composition="none; hold-out raises to 6 seq when n>=11, else refused",
             matrices="816"),
    ]).to_csv(REPORT / "method_notes.tsv", sep="\t", index=False)
    print(pd.DataFrame(jac_rows).to_string(index=False), flush=True)
    print(pd.DataFrame(ret_rows).to_string(index=False), flush=True)


def main() -> None:
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    s81 = load_s81()
    if stage in ("hocomoco", "all"):
        run_hocomoco(s81)
    if stage in ("sea", "all"):
        run_sea(s81)
    if stage in ("monalisa", "all"):
        run_monalisa(s81)
    if stage in ("compare", "all"):
        run_compare(s81)


if __name__ == "__main__":
    main()
