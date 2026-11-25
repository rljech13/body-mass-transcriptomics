#!/usr/bin/env bash
# Fibroblast CREate for the seven birds.
# Up to MAX_JOBS species at once, PARALLEL threads inside each CREate.
# A species already running (the first chicken call) is left alone.
# Input: bowtie2 5'-end CTSS from the GSE213253 remap (not CREate bam2ctss).
# Call matches the chicken reference mode: window 100, permissive
# (tpm >= 0.05, counts >= 4, length >= 5), PLA/ELA at 2 TPM.
set -euo pipefail

source /path/to/user-home/miniconda3/etc/profile.d/conda.sh
conda activate cage
export PATH="/path/to/storage/chicken/Dosage_compensation_proj/data/tools/bin:${PATH}"

CREATE=/path/to/storage/chicken/Dosage_compensation_proj/data/tools/src/CREate/CREate.sh
CTSS_ROOT=/path/to/storage/kleiber_birds_gemcat/cage_remap_mt/ctss
GENOME_ROOT=/path/to/storage/kleiber_birds_gemcat/ref/genomes_with_mt
OUT=/path/to/user-home/chicken/kleiber_birds_gemcat/results/oxphos_species_thread/create_fib
TMPDIR=/path/to/storage/kleiber_birds_gemcat/tmp_create
export TMPDIR
mkdir -p "$OUT" "$TMPDIR"

PARALLEL=${PARALLEL:-8}
MAX_JOBS=${MAX_JOBS:-3}
SPECIES=(Chicken Quail Duck Zebra Turkey Emu Peacock)

log() { echo "[$(date '+%F %T')] $*" | tee -a "$OUT/run.log"; }

already_running() {
  local sp=$1
  pgrep -f "create_fib/${sp}/${sp}( |$)" >/dev/null 2>&1
}

run_one() {
  local sp=$1
  local dest="$OUT/$sp"
  local prefix="$dest/${sp}"
  mkdir -p "$dest"
  if [[ -s "${prefix}_region_filter_classify.bed.gz" ]]; then
    log "SKIP $sp classified bed exists"
    return 0
  fi
  if already_running "$sp"; then
    log "SKIP $sp already running"
    return 0
  fi
  local files=()
  mapfile -t files < <(find "$CTSS_ROOT/$sp" -name '*.ctss.bed.gz' | sort)
  if [[ ${#files[@]} -eq 0 ]]; then
    log "ERROR $sp no CTSS"
    return 1
  fi
  cut -f1,2 "$GENOME_ROOT/$sp/genome.fa.fai" > "$dest/chrom.sizes"
  local comma
  comma=$(IFS=,; echo "${files[*]}")
  log "START $sp n=${#files[@]} threads=$PARALLEL jobs_cap=$MAX_JOBS"
  if ! "$CREATE" run \
      -i "$comma" \
      -c "$dest/chrom.sizes" \
      -o "$prefix" \
      -w 100 \
      -t 0.05 \
      -n 4 \
      -l 5 \
      -e 2 \
      -p "$PARALLEL" \
      >> "$dest/create.log" 2>&1
  then
    log "ERROR $sp CREate failed"
    return 1
  fi
  log "DONE $sp"
}

for sp in "${SPECIES[@]}"; do
  if [[ -s "$OUT/$sp/${sp}_region_filter_classify.bed.gz" ]] || already_running "$sp"; then
    log "SKIP $sp before launch"
    continue
  fi
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    sleep 20
  done
  run_one "$sp" &
done

wait
log "ALL DONE"
