#!/usr/bin/env bash
# Fetch the exact assemblies the colleagues' CAGE peaks were called on.
#
# Identified 20.09.2026 from the peak tables themselves:
#   Turkey  HG999xxx.1 + CAJPIY01xxxxxx.1  -> MGAL_WU_HG_1.0  GCA_905368555.1
#   Emu     CM0279xx.1 + JABVCD01xxxxxx.1  -> ZJU1.0          GCA_016128335.1
#   Peacock scaffoldN_lenL_covC (local)    -> AIIM_Pcri_1.0   GCA_005519975.1
#           (179332 scaffolds in the assembly vs max scaffold index 179339 in the
#            peak names; ABySS short-read assembly, N50 190 kb — confirmed by
#            matching the lengths encoded in the scaffold names, script 49)
#
# Turkey and Emu are fetched by exact versioned sequence accession from ENA, so
# no assembly-version ambiguity is possible. One accession per request: ENA
# refuses comma-separated batches of whole chromosomes. Resumable — an accession
# whose file already exists is skipped.
set -u

ROOT=/path/to/user-home/chicken/kleiber_birds_gemcat
EXP=$ROOT/data/dropbox_20260918/exp_tab
REF=/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies

fetch_by_accession() {
  species=$1
  dir=$REF/$species
  mkdir -p "$dir/seq"
  if [ -s "$dir/genome.fa" ]; then echo "[$species] genome.fa уже собран, пропускаю"; return 0; fi
  awk -F'\t' 'NR>1{print $2}' "$EXP/${species}_fib_epi_exp.txt" | sort -u > "$dir/accessions.txt"
  n=$(wc -l < "$dir/accessions.txt"); i=0; bad=0
  echo "[$species] последовательностей: $n"
  while read -r acc; do
    i=$((i+1))
    f="$dir/seq/${acc}.fa"
    if [ -s "$f" ] && head -1 "$f" | grep -q '^>'; then continue; fi
    if curl -sf --retry 5 --retry-delay 5 --max-time 900 \
         "https://www.ebi.ac.uk/ena/browser/api/fasta/${acc}?download=false" -o "$f.tmp" \
       && [ -s "$f.tmp" ] && head -1 "$f.tmp" | grep -q '^>'; then
      mv "$f.tmp" "$f"
    else
      rm -f "$f.tmp"; bad=$((bad+1)); echo "[$species] НЕ СКАЧАЛОСЬ: $acc"
    fi
    [ $((i % 10)) -eq 0 ] && echo "[$species] $i/$n"
  done < "$dir/accessions.txt"
  if [ "$bad" -gt 0 ]; then echo "[$species] пропущено $bad последовательностей, genome.fa не собираю"; return 1; fi
  cat "$dir"/seq/*.fa > "$dir/genome.fa"
  rm -f "$dir"/batch.*
  echo "[$species] записей в FASTA: $(grep -c '^>' "$dir/genome.fa")"
}

fetch_by_accession Turkey
fetch_by_accession Emu

if [ ! -s "$REF/Peacock/genome.fa" ]; then
  echo "[Peacock] качаю GCA_005519975.1 (AIIM_Pcri_1.0)"
  mkdir -p "$REF/Peacock"; cd "$REF/Peacock"
  datasets download genome accession GCA_005519975.1 --include genome --filename pavo.zip \
    && unzip -o -q pavo.zip \
    && find ncbi_dataset -name '*.fna' -exec mv {} genome.fa \; \
    && rm -rf pavo.zip ncbi_dataset README.md md5sum.txt \
    && echo "[Peacock] записей в FASTA: $(grep -c '^>' genome.fa)"
fi

echo "=== готово ==="
for s in Turkey Emu Peacock; do
  if [ -s "$REF/$s/genome.fa" ]; then
    printf '%-9s %s  записей %s\n' "$s" "$(du -h "$REF/$s/genome.fa" | cut -f1)" "$(grep -c '^>' "$REF/$s/genome.fa")"
  else
    printf '%-9s НЕ СОБРАН\n' "$s"
  fi
done
