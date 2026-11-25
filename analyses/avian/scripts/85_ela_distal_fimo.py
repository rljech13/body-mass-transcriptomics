#!/usr/bin/env python3
"""Motifs on distal ELA of the story sets, the layer HEY2 came from.

Ruslan's enhancer PDF (TFBS enh to sizes.pdf) is FIMO on bidirectional
enhancers linked to his 133, Fisher against other enhancers, then the
top 10 names per species, then -log10(p) against log10 mass. HEY2 and
MEF2B are that second step.

Here the positive sequences are CREate ELA whose nearest transcript TSS
is a story gene and the core sits 2–100 kb away. The window is expanded
to 400 bp, as on the proximal FIMO. Background is a fixed sample of
distal ELA assigned to other genes. BH is inside set x species.

The top-10-then-correlate step is also computed and written apart, so it
can be compared with his picture without being treated as a finding.
"""

from __future__ import annotations

import gzip
import importlib.util
import random
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import fisher_exact, pearsonr

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
CREATE = ROOT / "results/oxphos_species_thread/create_fib"
PEAKS = ROOT / "results/oxphos_species_thread/create_to_genes/peaks_nearest_tss.tsv"
OUT = ROOT / "results/oxphos_species_thread/create_to_genes/ela_fimo"
GENOME = Path("/path/to/storage/kleiber_birds_gemcat/ref/genomes_with_mt")
MEME = (
    ROOT
    / "results/oxphos_species_thread/module_consistency/fimo_chicken"
    / "JASPAR2024_vertebrates.meme"
)
FIMO = Path("/path/to/user-home/miniconda3/envs/meme/bin/fimo")
BEDTOOLS = Path("/path/to/user-home/miniconda3/envs/cage/bin/bedtools")

SPECIES = ["Zebra", "Quail", "Chicken", "Duck", "Peacock", "Turkey", "Emu"]
MASS = {
    "Zebra": 11.93,
    "Quail": 96.6,
    "Chicken": 912.5,
    "Duck": 1121.0,
    "Peacock": 4093.75,
    "Turkey": 5811.0,
    "Emu": 36200.0,
}
LOGM = np.log10([MASS[s] for s in SPECIES])
SETS = ["oxphos_16", "folate_3", "ribo_5", "green_hubs_7", "green_51"]
WATCH = ["HEY2", "MEF2B", "SNAI3", "PHOX2B", "ZEB1", "SOX14", "ESRRA"]
WIN = 400
BG_N = 4000
THRESH = "1e-4"


def load_s84():
    path = ROOT / "scripts/84_create_peaks_to_story_genes.py"
    spec = importlib.util.spec_from_file_location("s84", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def chrom_lengths(species: str) -> dict[str, int]:
    lengths = {}
    with (GENOME / species / "genome.fa.fai").open() as handle:
        for line in handle:
            chrom, length, *_ = line.split("\t")
            lengths[chrom] = int(length)
    return lengths


def window(chrom: str, core_start: int, core_end: int, lengths: dict[str, int]):
    length = lengths.get(chrom)
    if length is None:
        return None
    mid = (core_start + core_end) // 2
    half = WIN // 2
    start = max(0, mid - half)
    end = min(length, mid + half)
    if end - start < 50:
        return None
    return chrom, start, end


def positives() -> dict[str, list[tuple]]:
    """Distal ELA already assigned to a story gene."""
    out: dict[str, list[tuple]] = {sp: [] for sp in SPECIES}
    with PEAKS.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            rec = dict(zip(header, line.rstrip("\n").split("\t")))
            if rec["cre_class"] != "ELA" or rec["shelf"] != "distal":
                continue
            out[rec["species"]].append(
                (
                    rec["chrom"],
                    int(rec["core_start"]),
                    int(rec["core_end"]),
                    rec["gene"],
                    rec["sets"],
                )
            )
    return out


def background(s84, story: set[str]) -> dict[str, list[tuple]]:
    rng = random.Random(22)
    out: dict[str, list[tuple]] = {}
    for species in SPECIES:
        print(f"background {species}", flush=True)
        tss = s84.load_tss(species)
        reservoir: list[tuple] = []
        seen = 0
        bed = CREATE / species / f"{species}_region_filter_classify.bed.gz"
        with gzip.open(bed, "rt") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8 or "class:ELA" not in parts[3]:
                    continue
                rows = tss.get(parts[0])
                if not rows:
                    continue
                core_start, core_end = int(parts[6]), int(parts[7])
                if core_end <= core_start:
                    core_start, core_end = int(parts[1]), int(parts[2])
                hit = s84.nearest(rows, core_start, core_end)
                if hit is None:
                    continue
                gene, _tss, dist = hit
                if not (s84.PROX < dist <= s84.DIST) or gene in story:
                    continue
                seen += 1
                item = (parts[0], core_start, core_end, gene)
                if len(reservoir) < BG_N:
                    reservoir.append(item)
                else:
                    j = rng.randrange(seen)
                    if j < BG_N:
                        reservoir[j] = item
        out[species] = reservoir
        print(f"  distal ELA outside the sets {seen}, kept {len(reservoir)}", flush=True)
    return out


def write_fasta(species: str, pos: list[tuple], bg: list[tuple]) -> Path:
    lengths = chrom_lengths(species)
    bed_path = OUT / f"{species}.windows.bed"
    fa_path = OUT / f"{species}.fa"
    with bed_path.open("w") as bed:
        n = 0
        for chrom, a, b, gene, sets in pos:
            win = window(chrom, a, b, lengths)
            if win is None:
                continue
            name = f"pos|{n}|{gene}|{sets}"
            bed.write(f"{win[0]}\t{win[1]}\t{win[2]}\t{name}\n")
            n += 1
        for chrom, a, b, gene in bg:
            win = window(chrom, a, b, lengths)
            if win is None:
                continue
            name = f"bg|{n}|{gene}|."
            bed.write(f"{win[0]}\t{win[1]}\t{win[2]}\t{name}\n")
            n += 1
    genome = GENOME / species / "genome.fa"
    subprocess.check_call(
        [str(BEDTOOLS), "getfasta", "-fi", str(genome), "-bed", str(bed_path), "-name", "-fo", str(fa_path)],
        stderr=subprocess.DEVNULL,
    )
    return fa_path


def run_fimo(species: str, fasta: Path) -> Path:
    tsv = OUT / f"{species}.fimo.tsv"
    if tsv.exists() and tsv.stat().st_size > 100:
        print(f"fimo {species} exists", flush=True)
        return tsv
    print(f"fimo {species}", flush=True)
    cmd = [str(FIMO), "--text", "--no-qvalue", "--verbosity", "1", "--thresh", THRESH, str(MEME), str(fasta)]
    with tsv.open("w") as handle:
        subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.DEVNULL)
    return tsv


def sequence_sets(fasta: Path) -> tuple[dict[str, set[str]], set[str]]:
    """Name -> set membership. Background names have an empty set."""
    membership: dict[str, set[str]] = {}
    with fasta.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            name = line[1:].strip().split("::")[0]
            kind, _i, _gene, sets = name.split("|", 3)
            membership[name] = set(sets.split(",")) if kind == "pos" and sets != "." else set()
    return membership


def bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [1.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        val = min(running, pvals[i] * m / rank)
        running = val
        adj[i] = min(val, 1.0)
    return adj


def fisher_species(species: str, membership: dict[str, set[str]], hits: dict[str, set[str]]) -> list[dict]:
    rows = []
    motifs = sorted(hits)
    for set_name in SETS:
        pos = {name for name, sets in membership.items() if set_name in sets}
        neg = {name for name, sets in membership.items() if set_name not in sets and not sets}
        if len(pos) < 3 or len(neg) < 20:
            continue
        raw = []
        packed = []
        for motif in motifs:
            hit = hits[motif]
            pp = len(pos & hit)
            pn = len(neg & hit)
            table = np.array([[pp, len(pos) - pp], [pn, len(neg) - pn]], int)
            odds, p = fisher_exact(table, alternative="two-sided")
            raw.append(float(p))
            packed.append((motif, pp, pn, float(odds), float(p)))
        adjusted = bh(raw)
        for (motif, pp, pn, odds, p), q in zip(packed, adjusted):
            rows.append(dict(
                species=species,
                set=set_name,
                motif=motif,
                pos_n=len(pos),
                neg_n=len(neg),
                pos_hit=pp,
                neg_hit=pn,
                odds=odds,
                p=p,
                bh=q,
                logp=-np.log10(p) if p > 0 else 300.0,
            ))
    return rows


def load_hits(tsv: Path) -> dict[str, set[str]]:
    hits: dict[str, set[str]] = defaultdict(set)
    with tsv.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            rec = dict(zip(header, line.rstrip("\n").split("\t")))
            name = rec["sequence_name"].split("::")[0]
            alt = rec.get("motif_alt_id") or rec["motif_id"]
            hits[alt].add(name)
    return hits


def mass_rows(records: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str], dict[str, float]] = defaultdict(dict)
    bh_n: dict[tuple[str, str], int] = defaultdict(int)
    for rec in records:
        grouped[(rec["set"], rec["motif"])][rec["species"]] = rec["logp"]
        if rec["bh"] < 0.05:
            bh_n[(rec["set"], rec["motif"])] += 1
    out = []
    for (set_name, motif), by_sp in grouped.items():
        xs = []
        ys = []
        for sp, mass in zip(SPECIES, LOGM):
            if sp not in by_sp:
                continue
            xs.append(mass)
            ys.append(by_sp[sp])
        r = p = float("nan")
        if len(xs) >= 5 and np.std(ys) > 0:
            r, p = pearsonr(xs, ys)
        out.append(dict(
            set=set_name,
            motif=motif,
            n_species=len(xs),
            n_bh=bh_n[(set_name, motif)],
            r=float(r),
            p=float(p),
            **{sp: by_sp.get(sp, float("nan")) for sp in SPECIES},
        ))
    return out


def top10_step(records: list[dict]) -> list[dict]:
    """His enhancer picture: union of top 10 raw p per species, then correlate."""
    out = []
    for set_name in SETS:
        chosen = set()
        by_sp: dict[str, list[dict]] = defaultdict(list)
        for rec in records:
            if rec["set"] == set_name:
                by_sp[rec["species"]].append(rec)
        for sp, rows in by_sp.items():
            rows.sort(key=lambda rec: rec["p"])
            for rec in rows[:10]:
                chosen.add(rec["motif"])
        lookup = {(rec["species"], rec["motif"]): rec["logp"] for rec in records if rec["set"] == set_name}
        for motif in sorted(chosen):
            xs = []
            ys = []
            for sp, mass in zip(SPECIES, LOGM):
                if (sp, motif) not in lookup:
                    continue
                xs.append(mass)
                ys.append(lookup[(sp, motif)])
            r = p = float("nan")
            if len(xs) >= 5 and np.std(ys) > 0:
                r, p = pearsonr(xs, ys)
            out.append(dict(set=set_name, motif=motif, r=float(r), p=float(p), n=len(xs)))
    return out


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys = list(rows[0])
    with path.open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for rec in rows:
            handle.write("\t".join(str(rec[k]) for k in keys) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    s84 = load_s84()
    membership_genes, _order = s84.load_sets()
    story = set(membership_genes)
    pos = positives()
    bg = background(s84, story)
    records: list[dict] = []
    for species in SPECIES:
        fasta = write_fasta(species, pos[species], bg[species])
        tsv = run_fimo(species, fasta)
        membership = sequence_sets(fasta)
        hits = load_hits(tsv)
        print(f"  {species} sequences {len(membership)} motifs {len(hits)}", flush=True)
        records.extend(fisher_species(species, membership, hits))
    write_tsv(OUT / "fisher.tsv", records)
    slopes = mass_rows(records)
    write_tsv(OUT / "logp_vs_mass.tsv", slopes)
    tops = top10_step(records)
    write_tsv(OUT / "top10_then_mass.tsv", tops)
    print("BH < 0.05 in >= 3 species:")
    for rec in sorted(slopes, key=lambda r: (-r["n_bh"], r["p"] if np.isfinite(r["p"]) else 1)):
        if rec["n_bh"] >= 3:
            print(f"  {rec['set']} {rec['motif']} species_bh {rec['n_bh']} r {rec['r']:.3f} p {rec['p']:.3g}")
    print("watch:")
    for rec in slopes:
        if rec["motif"] in WATCH:
            print(f"  {rec['set']} {rec['motif']} bh_species {rec['n_bh']} r {rec['r']:.3f} p {rec['p']:.3g}")
    print("top10 step, watch or p<0.05:")
    for rec in tops:
        if rec["motif"] in WATCH or (np.isfinite(rec["p"]) and rec["p"] < 0.05):
            print(f"  {rec['set']} {rec['motif']} r {rec['r']:.3f} p {rec['p']:.3g}")


if __name__ == "__main__":
    main()
