#!/usr/bin/env python3
"""Assign fibroblast CREate regions to the story gene sets.

Anchor is the transcript TSS already written for the remap genome
(bidir_enh/{species}/tss.bed, column 5). A region goes to the nearest
TSS of any gene. Distance is from that TSS to the CREate core
(thickStart–thickEnd). Two different genes at the same distance are
left unassigned.

Shelves, matching the CAGE sum and the distal look already used:
  proximal  distance <= 2 kb
  distal    2 kb < distance <= 100 kb

PLA/ELA is the CREate class (TPM >= 2 vs below). This does not test
correlation across libraries and does not test body mass.
"""

from __future__ import annotations

import gzip
from bisect import bisect_left
from pathlib import Path

CREATE = Path(
    "/path/to/user-home/chicken/kleiber_birds_gemcat/results/"
    "oxphos_species_thread/create_fib"
)
TSS_ROOT = Path("/path/to/storage/kleiber_birds_gemcat/bidir_enh")
GENESETS = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/genesets")
OUT = Path(
    "/path/to/user-home/chicken/kleiber_birds_gemcat/results/"
    "oxphos_species_thread/create_to_genes"
)

SPECIES = ["Zebra", "Quail", "Chicken", "Duck", "Peacock", "Turkey", "Emu"]
PROX = 2000
DIST = 100_000
ALIAS = {
    "ATP5E": "ATP5F1E",
    "ATP5A1": "ATP5F1A",
    "ATP5B": "ATP5F1B",
    "ATP5C1": "ATP5F1C",
    "ATP5D": "ATP5F1D",
}
SETS = [
    ("oxphos_16", "story_oxphos_16.txt"),
    ("folate_3", "story_folate_3.txt"),
    ("ribo_5", "story_ribo_5.txt"),
    ("green_hubs_7", "story_green_hubs_7.txt"),
    ("green_51", "story_green_sig_51.txt"),
]


def norm(symbol: str) -> str:
    s = symbol.strip().upper()
    return ALIAS.get(s, s)


def load_sets() -> tuple[dict[str, list[str]], list[str]]:
    membership: dict[str, list[str]] = {}
    order: list[str] = []
    for name, filename in SETS:
        for line in (GENESETS / filename).read_text().splitlines():
            gene = norm(line)
            if not gene:
                continue
            membership.setdefault(gene, [])
            if name not in membership[gene]:
                membership[gene].append(name)
            if gene not in order:
                order.append(gene)
    return membership, order


def load_tss(species: str) -> dict[str, list[tuple[int, str]]]:
    by_chrom: dict[str, list[tuple[int, str]]] = {}
    path = TSS_ROOT / species / "tss.bed"
    seen: set[tuple[str, int, str]] = set()
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, gene, tss_s = parts[0], norm(parts[3]), parts[4]
            if not gene:
                continue
            tss = int(tss_s)
            key = (chrom, tss, gene)
            if key in seen:
                continue
            seen.add(key)
            by_chrom.setdefault(chrom, []).append((tss, gene))
    for chrom, rows in by_chrom.items():
        rows.sort()
    return by_chrom


def nearest(tss_rows: list[tuple[int, str]], start: int, end: int):
    """Nearest TSS to the core interval. None if two genes tie."""
    if not tss_rows or end <= start:
        return None
    positions = [row[0] for row in tss_rows]
    left = bisect_left(positions, start) - 1
    right = bisect_left(positions, end)
    inside = tss_rows[left + 1 : right]
    if inside:
        genes = {gene for _tss, gene in inside}
        if len(genes) == 1:
            tss, gene = min(inside, key=lambda row: abs(row[0] - (start + end) / 2))
            return gene, tss, 0
        return None
    cands = []
    if left >= 0:
        tss, gene = tss_rows[left]
        cands.append((start - tss, tss, gene))
    if right < len(tss_rows):
        tss, gene = tss_rows[right]
        cands.append((tss - end, tss, gene))
    if not cands:
        return None
    cands.sort()
    if len(cands) == 2 and cands[0][2] != cands[1][2] and cands[0][0] == cands[1][0]:
        return None
    dist, tss, gene = cands[0]
    return gene, tss, dist


def class_of(name: str) -> str:
    if "class:PLA" in name:
        return "PLA"
    if "class:ELA" in name:
        return "ELA"
    return ""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    membership, gene_order = load_sets()
    story = set(membership)
    peak_path = OUT / "peaks_nearest_tss.tsv"
    locus_path = OUT / "locus_present.tsv"
    summary_path = OUT / "summary_by_gene.tsv"

    with peak_path.open("w") as peaks, locus_path.open("w") as loci, summary_path.open("w") as summary:
        peaks.write(
            "species\tgene\tsets\tcre_class\tshelf\tdist_bp\ttpm\t"
            "chrom\tstart\tend\tcore_start\tcore_end\ttss\n"
        )
        loci.write("species\tgene\tsets\tn_tss\n")
        summary.write(
            "species\tgene\tsets\tn_tss\t"
            "prox_pla\tprox_ela\tprox_pla_tpm\tprox_ela_tpm\t"
            "dist_pla\tdist_ela\tdist_pla_tpm\tdist_ela_tpm\n"
        )
        for species in SPECIES:
            print(species, flush=True)
            tss = load_tss(species)
            n_tss = {gene: 0 for gene in story}
            for rows in tss.values():
                for _pos, gene in rows:
                    if gene in n_tss:
                        n_tss[gene] += 1
            for gene in gene_order:
                loci.write(
                    f"{species}\t{gene}\t{','.join(membership[gene])}\t{n_tss[gene]}\n"
                )

            counts = {
                gene: {"prox": {"PLA": [0, 0.0], "ELA": [0, 0.0]},
                       "dist": {"PLA": [0, 0.0], "ELA": [0, 0.0]}}
                for gene in story
            }
            n_peaks = n_hit = n_story = n_tie = n_far = 0
            bed = CREATE / species / f"{species}_region_filter_classify.bed.gz"
            with gzip.open(bed, "rt") as handle:
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 8:
                        continue
                    n_peaks += 1
                    chrom = parts[0]
                    rows = tss.get(chrom)
                    if not rows:
                        continue
                    start, end = int(parts[1]), int(parts[2])
                    core_start, core_end = int(parts[6]), int(parts[7])
                    if core_end <= core_start:
                        core_start, core_end = start, end
                    hit = nearest(rows, core_start, core_end)
                    if hit is None:
                        n_tie += 1
                        continue
                    gene, tss_pos, dist = hit
                    if dist > DIST:
                        n_far += 1
                        continue
                    n_hit += 1
                    if gene not in story:
                        continue
                    n_story += 1
                    klass = class_of(parts[3])
                    if klass not in ("PLA", "ELA"):
                        continue
                    try:
                        tpm = float(parts[4])
                    except ValueError:
                        tpm = 0.0
                    shelf = "proximal" if dist <= PROX else "distal"
                    key = "prox" if shelf == "proximal" else "dist"
                    counts[gene][key][klass][0] += 1
                    counts[gene][key][klass][1] += tpm
                    peaks.write(
                        f"{species}\t{gene}\t{','.join(membership[gene])}\t{klass}\t"
                        f"{shelf}\t{dist}\t{tpm:.6g}\t{chrom}\t{start}\t{end}\t"
                        f"{core_start}\t{core_end}\t{tss_pos}\n"
                    )
            for gene in gene_order:
                c = counts[gene]
                summary.write(
                    f"{species}\t{gene}\t{','.join(membership[gene])}\t{n_tss[gene]}\t"
                    f"{c['prox']['PLA'][0]}\t{c['prox']['ELA'][0]}\t"
                    f"{c['prox']['PLA'][1]:.6g}\t{c['prox']['ELA'][1]:.6g}\t"
                    f"{c['dist']['PLA'][0]}\t{c['dist']['ELA'][0]}\t"
                    f"{c['dist']['PLA'][1]:.6g}\t{c['dist']['ELA'][1]:.6g}\n"
                )
            print(
                f"  peaks {n_peaks} assigned<={DIST} {n_hit} "
                f"story {n_story} ties {n_tie} beyond {n_far}",
                flush=True,
            )
    print("wrote", peak_path)


if __name__ == "__main__":
    main()
