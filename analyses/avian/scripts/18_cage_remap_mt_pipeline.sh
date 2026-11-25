#!/usr/bin/env bash
# Remap GSE213253 fib CAGE onto genomes_with_mt (chrM included).
# Resume-safe. Logs to WORK/logs/.
set -euo pipefail

THREADS="${THREADS:-16}"
PARALLEL_JOBS="${PARALLEL_JOBS:-4}"  # concurrent samples (download/map)
WORK=/path/to/storage/kleiber_birds_gemcat/cage_remap_mt
GENOME_ROOT=/path/to/storage/kleiber_birds_gemcat/ref/genomes_with_mt
MANIFEST=/path/to/user-home/chicken/kleiber_birds_gemcat/results/cage_remap_mt/gse213253_cage_sra_manifest.tsv
MT_BED_DIR="$WORK/mt_annotation"
LOG="$WORK/logs"
mkdir -p "$WORK"/{fastq,index,bam,ctss,quant,logs} "$MT_BED_DIR"

ts() { date '+%F %T'; }
log() { echo "[$(ts)] $*" | tee -a "$LOG/pipeline.log"; }

# species -> MT accession used in genome (for gene BED)
declare -A MT_ACC=(
  [Zebra]=NC_007897.1
  [Quail]=NC_003408.1
  [Chicken]=NC_001323.1
  [Duck]=NC_009684.1
  [Peacock]=NC_024533.1
  [Turkey]=NC_010195.2
  [Emu]=NC_002784.1
)

# ---------- MT gene BEDs (chrM) ----------
make_mt_bed() {
  local sp=$1 acc=$2
  local bed="$MT_BED_DIR/${sp}_chrM_genes.bed"
  [[ -s "$bed" ]] && return 0
  local ft="$MT_BED_DIR/${acc}.ft"
  if [[ ! -s "$ft" ]]; then
    curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=ft&retmode=text" -o "$ft" || true
  fi
  python3 - "$ft" "$bed" "$acc" <<'PY'
import re, sys
from pathlib import Path
ft, bed, acc = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
text = ft.read_text() if ft.exists() else ""
# Parse NCBI feature table: gene / CDS with product or gene
rows = []
blocks = re.split(r'\n(?=\S)', text)
cur = None
for b in blocks:
    lines = b.strip().splitlines()
    if not lines:
        continue
    head = lines[0]
    m = re.match(r'^(<?\d+)\t(>?\d+)\t(gene|CDS)$', head)
    if not m:
        # sometimes coords on next style: "     1\t16775\tgene"
        m = re.match(r'^(\d+)\t(\d+)\t(gene|CDS)$', head.strip())
    if not m:
        continue
    a, b2, ftype = m.group(1), m.group(2), m.group(3)
    strand = '+'
    start, end = a, b2
    if start.startswith('<'): start = start[1:]
    if end.startswith('>'): end = end[1:]
    try:
        s, e = int(start), int(end)
    except ValueError:
        continue
    if s > e:
        strand = '-'
        s, e = e, s
    gene = None
    product = None
    for ln in lines[1:]:
        if '\tgene\t' in ln:
            gene = ln.split('\t')[-1].strip()
        if '\tproduct\t' in ln:
            product = ln.split('\t')[-1].strip()
        if '\tgene_synonym\t' in ln and not gene:
            gene = ln.split('\t')[-1].strip()
    # normalize names
    name = (gene or product or '').upper()
    alias = {
        'CYTOCHROME B': 'CYTB', 'CYTOCHROME C OXIDASE SUBUNIT I': 'COX1',
        'CYTOCHROME C OXIDASE SUBUNIT II': 'COX2', 'CYTOCHROME C OXIDASE SUBUNIT III': 'COX3',
        'ATP SYNTHASE F0 SUBUNIT 6': 'ATP6', 'ATP SYNTHASE F0 SUBUNIT 8': 'ATP8',
        'NADH DEHYDROGENASE SUBUNIT 1': 'ND1', 'NADH DEHYDROGENASE SUBUNIT 2': 'ND2',
        'NADH DEHYDROGENASE SUBUNIT 3': 'ND3', 'NADH DEHYDROGENASE SUBUNIT 4': 'ND4',
        'NADH DEHYDROGENASE SUBUNIT 4L': 'ND4L', 'NADH DEHYDROGENASE SUBUNIT 5': 'ND5',
        'NADH DEHYDROGENASE SUBUNIT 6': 'ND6',
    }
    if name in alias:
        name = alias[name]
    for k, v in alias.items():
        if k in name:
            name = v
            break
    # keep standard mito genes + rRNA
    keep = {'ND1','ND2','ND3','ND4','ND4L','ND5','ND6','CYTB','COX1','COX2','COX3','ATP6','ATP8',
            '12S','16S','RNR1','RNR2'}
    if ftype == 'gene' and name:
        # map rRNA
        if '12S' in name or 'RNR1' in name: name = 'RNR1'
        if '16S' in name or 'RNR2' in name: name = 'RNR2'
        if name in keep or name.startswith('ND') or name.startswith('COX') or name.startswith('ATP'):
            rows.append((s-1, e, name, strand, ftype))
# prefer gene features; dedup by name keeping longest
best = {}
for s,e,name,strand,ftype in rows:
    if name not in best or (e-s) > (best[name][1]-best[name][0]):
        best[name] = (s,e,name,strand)
with bed.open('w') as out:
    for s,e,name,strand in sorted(best.values()):
        out.write(f'chrM\t{s}\t{e}\t{name}\t0\t{strand}\n')
print(f'wrote {bed} n={len(best)}')
PY
}

# ---------- indexes ----------
build_index() {
  local sp=$1
  local fa="$GENOME_ROOT/$sp/genome.fa"
  local idx="$WORK/index/$sp"
  mkdir -p "$idx"
  if [[ -f "$idx/genome.1.bt2" || -f "$idx/genome.1.bt2l" ]]; then
    log "INDEX skip $sp"
    return 0
  fi
  [[ -s "$fa" ]] || { log "ERROR no genome $fa"; return 1; }
  log "INDEX build $sp (threads=$THREADS)"
  bowtie2-build --threads "$THREADS" "$fa" "$idx/genome" >"$LOG/index_${sp}.log" 2>&1
  log "INDEX done $sp"
}

# ---------- download ----------
download_fq() {
  local sp=$1 run=$2 url=$3
  local outdir="$WORK/fastq/$sp"
  mkdir -p "$outdir"
  local out="$outdir/${run}.fastq.gz"
  if [[ -s "$out" ]]; then
    # basic size check > 10MB
    local sz; sz=$(stat -c%s "$out")
    if (( sz > 10000000 )); then
      return 0
    fi
  fi
  log "DOWNLOAD $sp $run"
  wget -c -q -O "$out.partial" "$url" && mv "$out.partial" "$out"
  log "DOWNLOAD done $sp $run ($(stat -c%s "$out") bytes)"
}

# ---------- map + CTSS ----------
map_one() {
  local sp=$1 run=$2
  local fq="$WORK/fastq/$sp/${run}.fastq.gz"
  local idx="$WORK/index/$sp/genome"
  local bam="$WORK/bam/$sp/${run}.bam"
  local ctss="$WORK/ctss/$sp/${run}.ctss.bed.gz"
  mkdir -p "$WORK/bam/$sp" "$WORK/ctss/$sp"
  [[ -s "$fq" ]] || { log "ERROR missing fq $fq"; return 1; }
  if [[ -s "$bam" && -s "${bam}.bai" && -s "$ctss" ]]; then
    log "MAP skip $sp $run"
    return 0
  fi
  log "MAP $sp $run"
  bowtie2 -p "$THREADS" -x "$idx" -U "$fq" 2>"$LOG/map_${sp}_${run}.bowtie2.log" \
    | samtools view -bS -F 4 - \
    | samtools sort -@ 4 -o "$bam" -
  samtools index "$bam"
  # CTSS = 5' ends (subprocess so heredoc doesn't steal stdin)
  python3 - "$bam" "$ctss" <<'PY'
import gzip, re, subprocess, sys
bam_path, out_path = sys.argv[1], sys.argv[2]
cig_re = re.compile(r'(\d+)([MIDNSHP=X])')
counts = {}
proc = subprocess.Popen(['samtools', 'view', '-F', '4', bam_path], stdout=subprocess.PIPE, text=True)
assert proc.stdout is not None
for ln in proc.stdout:
    f = ln.split('\t')
    flag = int(f[1]); chrom = f[2]; pos = int(f[3]); cigar = f[5]
    reflen = sum(int(n) for n, op in cig_re.findall(cigar) if op in 'MDN=X')
    if flag & 16:
        p5 = pos + reflen - 1
        strand = '-'
    else:
        p5 = pos
        strand = '+'
    key = (chrom, p5 - 1, strand)
    counts[key] = counts.get(key, 0) + 1
proc.wait()
with gzip.open(out_path, 'wt') as out:
    for (chrom, pos, strand), c in sorted(counts.items()):
        out.write(f'{chrom}\t{pos}\t{pos+1}\t.\t{c}\t{strand}\n')
print(f'CTSS {out_path} sites={len(counts)} tags={sum(counts.values())}')
PY
  log "MAP done $sp $run"
}

quant_mt() {
  local sp=$1 run=$2
  local ctss="$WORK/ctss/$sp/${run}.ctss.bed.gz"
  local bed="$MT_BED_DIR/${sp}_chrM_genes.bed"
  local out="$WORK/quant/${sp}_${run}_mt_genes.tsv"
  [[ -s "$ctss" && -s "$bed" ]] || return 0
  [[ -s "$out" ]] && return 0
  # total tags + chrM tags + per-gene
  python3 - "$ctss" "$bed" "$out" "$sp" "$run" <<'PY'
import gzip, sys
from collections import defaultdict
ctss, bed, out, sp, run = sys.argv[1:6]
genes = []
with open(bed) as f:
    for ln in f:
        c,s,e,name,sc,strand = ln.rstrip().split('\t')[:6]
        genes.append((c,int(s),int(e),name,strand))
total=chrM=0
gene_counts=defaultdict(int)
with gzip.open(ctss,'rt') as f:
    for ln in f:
        chrom,s,e,name,score,strand = ln.rstrip().split('\t')[:6]
        s=int(s); score=int(float(score))
        total += score
        if chrom == 'chrM' or 'mito' in chrom.lower():
            chrM += score
            for gc,gs,ge,gn,gstr in genes:
                if s >= gs and s < ge:
                    gene_counts[gn] += score
with open(out,'w') as o:
    o.write('species\trun\tgene\ttags\ttotal_tags\tchrM_tags\tfrac_of_lib\tfrac_of_chrM\n')
    for gn, tags in sorted(gene_counts.items()):
        o.write(f'{sp}\t{run}\t{gn}\t{tags}\t{total}\t{chrM}\t{tags/max(total,1):.6g}\t{tags/max(chrM,1):.6g}\n')
    o.write(f'{sp}\t{run}\t__chrM_total__\t{chrM}\t{total}\t{chrM}\t{chrM/max(total,1):.6g}\t1\n')
print('quant', out, 'chrM', chrM, 'total', total)
PY
}

MAP_FN=map_one

log "=== START CAGE remap pipeline THREADS=$THREADS PARALLEL_JOBS=$PARALLEL_JOBS ==="

# 1) MT beds + indexes (serial indexes to avoid RAM spike; 7 genomes)
for sp in Zebra Quail Chicken Duck Peacock Turkey Emu; do
  make_mt_bed "$sp" "${MT_ACC[$sp]}"
done

for sp in Zebra Quail Chicken Duck Peacock Turkey Emu; do
  build_index "$sp"
done

# 2) read manifest and process samples with bounded parallelism
mapfile -t LINES < <(tail -n +2 "$MANIFEST")

download_all() {
  local pids=()
  for line in "${LINES[@]}"; do
    IFS=$'\t' read -r sp run title layout spots size ena ena1 sex <<<"$line"
    url="$ena1"
    [[ "$url" == "-" || -z "$url" ]] && url="$ena"
    while (( ${#pids[@]} >= PARALLEL_JOBS )); do
      new=()
      for pid in "${pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then new+=("$pid"); else wait "$pid" || true; fi
      done
      pids=("${new[@]}")
      (( ${#pids[@]} >= PARALLEL_JOBS )) && sleep 2
    done
    download_fq "$sp" "$run" "$url" &
    pids+=($!)
  done
  for pid in "${pids[@]}"; do wait "$pid" || true; done
}

map_all() {
  local pids=()
  # fewer parallel maps — each uses THREADS cores
  local map_parallel=$(( 256 / THREADS )); (( map_parallel < 1 )) && map_parallel=1
  (( map_parallel > 6 )) && map_parallel=6
  log "MAP parallel=$map_parallel"
  for line in "${LINES[@]}"; do
    IFS=$'\t' read -r sp run title layout spots size ena ena1 sex <<<"$line"
    while (( ${#pids[@]} >= map_parallel )); do
      new=()
      for pid in "${pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then new+=("$pid"); else wait "$pid" || true; fi
      done
      pids=("${new[@]}")
      (( ${#pids[@]} >= map_parallel )) && sleep 5
    done
    (
      "$MAP_FN" "$sp" "$run"
      quant_mt "$sp" "$run"
    ) &
    pids+=($!)
  done
  for pid in "${pids[@]}"; do wait "$pid" || true; done
}

log "=== DOWNLOAD phase ==="
download_all
log "=== MAP+CTSS+QUANT phase ==="
map_all

# 3) aggregate
python3 <<'PY'
from pathlib import Path
import pandas as pd
work = Path('/path/to/storage/kleiber_birds_gemcat/cage_remap_mt')
rows = sorted((work / 'quant').glob('*_mt_genes.tsv'))
rows = [p for p in rows if p.name.startswith(('Zebra','Quail','Chicken','Duck','Peacock','Turkey','Emu'))]
if not rows:
    print('No quant tables yet')
else:
    big = pd.concat([pd.read_csv(p, sep='\t') for p in rows], ignore_index=True)
    big.to_csv(work / 'quant' / 'all_runs_mt_genes.tsv', sep='\t', index=False)
    g = big[big.gene != '__chrM_total__']
    sp = g.groupby(['species', 'gene'], as_index=False).agg(
        mean_tags=('tags', 'mean'), mean_frac_lib=('frac_of_lib', 'mean'), n_runs=('run', 'nunique')
    )
    sp.to_csv(work / 'quant' / 'species_mean_mt_genes.tsv', sep='\t', index=False)
    lib = (
        big[big.gene == '__chrM_total__']
        .groupby('species', as_index=False)
        .agg(mean_chrM_frac=('frac_of_lib', 'mean'), mean_chrM_tags=('chrM_tags', 'mean'))
    )
    lib.to_csv(work / 'quant' / 'species_chrM_library_fraction.tsv', sep='\t', index=False)
    print('Wrote aggregates', len(big), 'rows from', len(rows), 'files')
PY

log "=== DONE ==="
log "Outputs: $WORK/quant/  BAMs: $WORK/bam/  CTSS: $WORK/ctss/"
