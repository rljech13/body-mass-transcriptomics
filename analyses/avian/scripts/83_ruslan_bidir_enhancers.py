#!/usr/bin/env python3
"""Transcribed enhancers the way Ruslan called them.

bidir_enhancers (Andersson et al. 2014), defaults: gap 400, window 200,
|directionality| < 0.8. Mask is his TxDb mask: transcript TSS ±3 kb and
exons ±200 bp. The ROC promoter classifier is not in birds - ESC.R, so it
is not reapplied.

CTSS are the fibroblast libraries already remapped to genomes_with_mt.
Epiblast CTSS is not on this host. After the call: keep enhancers with
signal in >1 library and sum of TPM > 1, then link each to a promoter tag
cluster whose center is within ±250 kb, Pearson r > 0 and p < 0.05.
"""

import argparse
import gzip
import os
import subprocess
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr

BEDTOOLS = "/path/to/user-home/miniconda3/envs/cage/bin/bedtools"
BIDIR = "/path/to/storage/kleiber_birds_gemcat/tools/enhancers/scripts/bidir_enhancers"
CTSS_ROOT = Path("/path/to/storage/kleiber_birds_gemcat/cage_remap_mt/ctss")
GENOME_ROOT = Path("/path/to/storage/kleiber_birds_gemcat/ref/genomes_with_mt")
GFF_DIR = Path("/path/to/storage/kleiber_birds_gemcat/bidir_enh/gff")
ALIAS = Path("/path/to/storage/kleiber_birds_gemcat/ref/galGal6.chromAlias.txt")
OXPHOS = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/genesets/story_oxphos_16.txt")
CHICKEN_GFF = Path(
    "/path/to/storage/kleiber_birds_gemcat/ref/_fetch_tmp/Chicken/"
    "ncbi_dataset/data/GCF_000002315.6/genomic.gff"
)

TRANSCRIPT = {
    "mRNA", "transcript", "lnc_RNA", "ncRNA", "miRNA", "tRNA", "snoRNA",
    "snRNA", "rRNA", "primary_transcript", "pseudogenic_transcript",
    "V_gene_segment", "C_gene_segment",
}
EXON = {"exon", "CDS"}


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def attrs_of(field):
    out = {}
    for part in field.split(";"):
        if "=" in part:
            key, val = part.split("=", 1)
            out[key] = val
    return out


def gene_symbol(attrs):
    target = attrs.get("Target", "")
    if "##" in target:
        sym = target.split("##", 1)[1].split()[0]
        if sym:
            return sym
    for key in ("gene", "Name", "ID"):
        if attrs.get(key):
            return attrs[key]
    return ""


def chicken_chrom_map():
    mapping = {}
    with ALIAS.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            ucsc, _ens, genbank, refseq = line.rstrip("\n").split("\t")[:4]
            mapping[genbank] = ucsc
            mapping[refseq] = ucsc
            mapping[ucsc] = ucsc
    return mapping


def fai_lengths(species):
    lengths = {}
    with (GENOME_ROOT / species / "genome.fa.fai").open() as handle:
        for line in handle:
            chrom, length = line.split("\t")[:2]
            lengths[chrom] = int(length)
    return lengths


def gff_for(species):
    if species == "Chicken":
        return CHICKEN_GFF, chicken_chrom_map()
    lifted = GFF_DIR / f"{species}.liftoff.gff"
    native = GFF_DIR / f"{species}.genomic.gff.gz"
    if lifted.exists() and lifted.stat().st_size > 0:
        return lifted, {}
    if native.exists():
        return native, {}
    raise FileNotFoundError(f"no GFF for {species}")


def build_annotation(species, dest):
    gff, chrom_map = gff_for(species)
    lengths = fai_lengths(species)
    mask_raw = dest / "mask.raw.bed"
    tss_raw = dest / "tss.raw.bed"
    seen_tss = set()
    n_tx = n_exon = 0
    with open_text(gff) as handle, mask_raw.open("w") as mask, tss_raw.open("w") as tss:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, _source, feature, start_s, end_s, _score, strand = parts[:7]
            chrom = chrom_map.get(chrom, chrom)
            if chrom not in lengths:
                continue
            if strand not in "+-":
                continue
            start = int(start_s)
            end = int(end_s)
            if feature in TRANSCRIPT:
                tss_1 = start if strand == "+" else end
                key = (chrom, tss_1, strand, gene_symbol(attrs_of(parts[8] if len(parts) > 8 else "")))
                if key in seen_tss:
                    continue
                seen_tss.add(key)
                bed_start = max(0, tss_1 - 1 - 3000)
                bed_end = min(lengths[chrom], tss_1 - 1 + 3000 + 1)
                if bed_end <= bed_start:
                    continue
                gene = key[3]
                mask.write(f"{chrom}\t{bed_start}\t{bed_end}\n")
                tss.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene}\t{tss_1 - 1}\t{strand}\n")
                n_tx += 1
            elif feature in EXON:
                bed_start = max(0, start - 1 - 200)
                bed_end = min(lengths[chrom], end + 200)
                if bed_end <= bed_start:
                    continue
                mask.write(f"{chrom}\t{bed_start}\t{bed_end}\n")
                n_exon += 1
    merged = dest / "mask.bed"
    tss_bed = dest / "tss.bed"
    subprocess.check_call(
        f"sort -k1,1 -k2,2n {mask_raw} | {BEDTOOLS} merge -i - > {merged}",
        shell=True,
    )
    subprocess.check_call(f"sort -k1,1 -k2,2n {tss_raw} > {tss_bed}", shell=True)
    mask_raw.unlink()
    tss_raw.unlink()
    return n_tx, n_exon


def prepare_ctss(species, dest):
    files = sorted((CTSS_ROOT / species).glob("*.ctss.bed.gz"))
    if not files:
        raise FileNotFoundError(species)
    plain_dir = dest / "ctss"
    plain_dir.mkdir(exist_ok=True)
    names = []
    paths = []
    for src in files:
        name = src.name.replace(".ctss.bed.gz", "")
        plain = plain_dir / f"{name}.bed"
        if not plain.exists() or plain.stat().st_size == 0:
            subprocess.check_call(
                f"zcat {src} | sort -k1,1 -k2,2n > {plain}",
                shell=True,
            )
        names.append(name)
        paths.append(plain)
    listing = dest / "ctss_paths.txt"
    listing.write_text("".join(f"{p}\n" for p in paths))
    (dest / "libraries.txt").write_text("".join(f"{n}\n" for n in names))
    return names, paths, listing


def run_bidir(dest, listing):
    out = dest / "call"
    tmp = dest / "tmp"
    out.mkdir(exist_ok=True)
    tmp.mkdir(exist_ok=True)
    enh = out / "enhancers.bed"
    if enh.exists() and enh.stat().st_size > 0:
        return enh
    env = os.environ.copy()
    env["PATH"] = "/path/to/user-home/miniconda3/envs/cage/bin:" + env.get("PATH", "")
    env["TMPDIR"] = str(tmp)
    log = dest / "bidir.log"
    with log.open("w") as handle:
        proc = subprocess.run(
            [
                "bash", BIDIR,
                "-f", str(listing),
                "-o", str(out),
                "-m", str(dest / "mask.bed"),
                "-t", str(tmp),
            ],
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0 or not enh.exists():
        raise RuntimeError(f"bidir_enhancers failed, see {log}")
    return enh


def read_tpm(path, n_libs):
    """id -> vector. The tool writes id then one column per library."""
    table = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 1 + n_libs:
                continue
            vals = np.array([float(x) if x not in ("NA", "") else np.nan for x in parts[1:1 + n_libs]])
            table[parts[0]] = vals
    return table


def chrom_sizes(species, dest):
    sizes = dest / "chrom.sizes"
    with (GENOME_ROOT / species / "genome.fa.fai").open() as inp, sizes.open("w") as out:
        for line in inp:
            chrom, length = line.split("\t")[:2]
            out.write(f"{chrom}\t{length}\n")
    return sizes


def sort_bed(src, dest, sizes):
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    with dest.open("w") as handle:
        subprocess.check_call([BEDTOOLS, "sort", "-g", str(sizes), "-i", str(src)], stdout=handle)
    return dest


def promoter_table(dest, libraries, ctss_paths, sizes):
    plus = dest / "tmp" / "TCs.plus.bed"
    minus = dest / "tmp" / "TCs.minus.bed"
    if not plus.exists() or not minus.exists():
        # stub is empty, files are TCs.plus.bed
        cands = list((dest / "tmp").glob("*TCs.plus.bed"))
        if not cands:
            raise FileNotFoundError("tag clusters were not kept")
        plus = cands[0]
        minus = Path(str(plus).replace("TCs.plus.bed", "TCs.minus.bed"))
    both = dest / "tcs.bed"
    subprocess.check_call(f"cat {plus} {minus} | sort -k1,1 -k2,2n > {both}", shell=True)
    promoters = dest / "promoters.bed"
    # tss.bed: chrom, winStart, winEnd, gene, tss0, strand
    # keep a TC if it overlaps a same-strand TSS window; gene = nearest TSS
    subprocess.check_call(
        f"{BEDTOOLS} intersect -s -wa -wb -a {both} -b {dest / 'tss.bed'} > {dest / 'tc_tss_hits.txt'}",
        shell=True,
    )
    best = {}
    with (dest / "tc_tss_hits.txt").open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            # TC bed6 + tss bed6
            if len(parts) < 12:
                continue
            tc_id = parts[3]
            tc_start = int(parts[1])
            tc_end = int(parts[2])
            mid = (tc_start + tc_end) / 2
            gene = parts[9]
            tss0 = int(parts[10])
            dist = abs(mid - tss0)
            prev = best.get(tc_id)
            if prev is None or dist < prev[0]:
                best[tc_id] = (dist, parts[:6], gene)
    with promoters.open("w") as handle:
        for tc_id, (_dist, cols, gene) in best.items():
            handle.write("\t".join(cols[:6]) + f"\t{gene}\n")
    sort_bed(promoters, dest / "promoters.genome.bed", sizes)
    promoters = dest / "promoters.genome.bed"
    # quantify
    ids = []
    genes = []
    with promoters.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            ids.append(parts[3])
            genes.append(parts[6] if len(parts) > 6 else "")
    if not ids:
        return [], [], None
    mat = np.zeros((len(ids), len(ctss_paths)))
    index = {tc: i for i, tc in enumerate(ids)}
    for j, ctss in enumerate(ctss_paths):
        total = 0.0
        with ctss.open() as handle:
            for line in handle:
                total += float(line.split("\t")[4])
        mapped = dest / f"pro_map_{j}.txt"
        ctss_sorted = sort_bed(ctss, dest / "ctss_genome" / ctss.name, sizes)
        subprocess.check_call(
            f"{BEDTOOLS} map -g {sizes} -s -c 5 -o sum -null 0 -a {promoters} -b {ctss_sorted} > {mapped}",
            shell=True,
        )
        with mapped.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                tc = parts[3]
                mat[index[tc], j] = float(parts[-1]) / total * 1e6 if total else 0.0
        mapped.unlink()
    return ids, genes, mat


def link(dest, n_libs, pro_ids, pro_genes, pro_mat, oxphos, sizes):
    enh_bed = dest / "call" / "enhancers.bed"
    tpm_path = dest / "call" / "enhancer.expression.tpm.matrix"
    enh_tpm = read_tpm(tpm_path, n_libs)
    # filter as in the R script
    keep = []
    with enh_bed.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            vec = enh_tpm.get(parts[3])
            if vec is None or np.isnan(vec).any():
                continue
            if np.sum(vec > 0) > 1 and np.sum(vec) > 1:
                center = int(parts[6])
                start = max(0, center - 250000)
                end = center + 250000
                keep.append((parts, start, end, vec))
    filtered = dest / "enhancers_1tpm.bed"
    windows = dest / "enh_windows.bed"
    with filtered.open("w") as bed, windows.open("w") as win:
        for parts, start, end, _vec in keep:
            bed.write("\t".join(parts) + "\n")
            win.write(f"{parts[0]}\t{start}\t{end}\t{parts[3]}\n")
    if not keep or pro_mat is None:
        return {"n_1tpm": len(keep), "n_pairs": 0, "n_p05": 0, "n_bh": 0,
                "n_genes": 0, "n_ox_genes": 0, "n_ox_links": 0}
    sort_bed(windows, dest / "enh_windows.genome.bed", sizes)
    hits = dest / "enh_pro_hits.txt"
    subprocess.check_call(
        f"{BEDTOOLS} intersect -g {sizes} -wa -wb -a {dest / 'enh_windows.genome.bed'} -b {dest / 'promoters.genome.bed'} > {hits}",
        shell=True,
    )
    enh_index = {parts[3]: vec for parts, _s, _e, vec in keep}
    pro_index = {tc: i for i, tc in enumerate(pro_ids)}
    rows = []
    pvals = []
    seen = set()
    with hits.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            # window bed4 + promoter bed7
            enh_id = parts[3]
            pro_id = parts[7]
            pair = (enh_id, pro_id)
            if pair in seen:
                continue
            seen.add(pair)
            if pro_id not in pro_index or enh_id not in enh_index:
                continue
            x = enh_index[enh_id]
            y = pro_mat[pro_index[pro_id]]
            if np.std(x) == 0 or np.std(y) == 0:
                continue
            r, p = pearsonr(x, y)
            if not np.isfinite(r) or not np.isfinite(p):
                continue
            gene = pro_genes[pro_index[pro_id]]
            rows.append((enh_id, pro_id, gene, r, p, parts[4], parts[5], parts[6]))
            pvals.append(p)
    n_pairs = len(rows)
    bh = np.ones(n_pairs)
    if n_pairs:
        order = np.argsort(pvals)
        ranked = np.array(pvals)[order]
        adj = np.empty(n_pairs)
        prev = 1.0
        for i in range(n_pairs - 1, -1, -1):
            prev = min(prev, ranked[i] * n_pairs / (i + 1))
            adj[i] = prev
        bh[order] = np.clip(adj, 0, 1)
    out = dest / "links.tsv"
    n_p05 = n_bh = n_ox_links = 0
    genes = set()
    ox_genes = set()
    with out.open("w") as handle:
        handle.write("enh_id\tpro_id\tgene\tr\tp\tbh\tchrom\tpro_start\tpro_end\n")
        for i, (enh_id, pro_id, gene, r, p, chrom, pro_start, pro_end) in enumerate(rows):
            if r > 0 and p < 0.05:
                handle.write(
                    f"{enh_id}\t{pro_id}\t{gene}\t{r:.6g}\t{p:.6g}\t{bh[i]:.6g}\t"
                    f"{chrom}\t{pro_start}\t{pro_end}\n"
                )
                n_p05 += 1
                genes.add(gene)
                if gene in oxphos:
                    n_ox_links += 1
                    ox_genes.add(gene)
            if bh[i] < 0.05 and r > 0:
                n_bh += 1
    return {
        "n_1tpm": len(keep),
        "n_pairs": n_pairs,
        "n_p05": n_p05,
        "n_bh": n_bh,
        "n_genes": len(genes),
        "n_ox_genes": len(ox_genes),
        "n_ox_links": n_ox_links,
    }


def count_lines(path):
    if not Path(path).exists():
        return 0
    n = 0
    with open(path) as handle:
        for _ in handle:
            n += 1
    return n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("species")
    parser.add_argument("--out", default="/path/to/storage/kleiber_birds_gemcat/bidir_enh")
    args = parser.parse_args()
    species = args.species
    dest = Path(args.out) / species
    dest.mkdir(parents=True, exist_ok=True)
    oxphos = {line.strip() for line in OXPHOS.read_text().splitlines() if line.strip()}
    print(f"annotation {species}", flush=True)
    n_tx, n_exon = build_annotation(species, dest)
    print(f"tss_windows {n_tx} exon_or_cds {n_exon}", flush=True)
    libraries, ctss_paths, listing = prepare_ctss(species, dest)
    sizes = chrom_sizes(species, dest)
    (dest / "ctss_genome").mkdir(exist_ok=True)
    print(f"libraries {len(libraries)}", flush=True)
    run_bidir(dest, listing)
    raw = count_lines(dest / "call" / "enhancers.bed")
    bidir_n = count_lines(dest / "call" / "bidir.pairs.filtered.bed")
    print(f"bidir_filtered {bidir_n} enhancers {raw}", flush=True)
    pro_ids, pro_genes, pro_mat = promoter_table(dest, libraries, ctss_paths, sizes)
    print(f"promoter_tcs {len(pro_ids)}", flush=True)
    stats = link(dest, len(libraries), pro_ids, pro_genes, pro_mat, oxphos, sizes)
    summary = dest / "summary.tsv"
    header = [
        "species", "n_libs", "n_tss_windows", "n_exon_mask", "n_bidir_filtered",
        "n_enhancers", "n_promoter_tc", "n_1tpm", "n_pairs", "n_links_p05",
        "n_links_bh05_rpos", "n_genes", "n_oxphos16_genes", "n_oxphos16_links",
    ]
    values = [
        species, len(libraries), n_tx, n_exon, bidir_n, raw, len(pro_ids),
        stats["n_1tpm"], stats["n_pairs"], stats["n_p05"], stats["n_bh"],
        stats["n_genes"], stats["n_ox_genes"], stats["n_ox_links"],
    ]
    summary.write_text("\t".join(header) + "\n" + "\t".join(str(v) for v in values) + "\n")
    print(summary.read_text(), flush=True)


if __name__ == "__main__":
    main()
