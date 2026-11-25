#!/usr/bin/env bash
# Assemblies and proteomes for the four species that never got a re-annotation.
#
# Script 58 projected a whole proteome only onto Peacock, Turkey and Emu, and gave
# each a different query. The other four kept loci that come from the peak table
# alone — a gene with no peak assigned to it never got a 5' coordinate — so their
# locus vocabulary stops at ~10 000 while the three re-annotated ones reach 16-18k.
# Genes with a locus in all seven number 6567, and that is the hard ceiling on the
# WGCNA panel. Raising it means giving every species the same treatment.
#
# The assemblies below are the ones the CAGE peaks were called on, identified from
# the contig names in the peak tables (NC_006088.5, NC_051772.1, NC_029516.1,
# NC_045000.1). What sits in ref/genomes_with_mt is NEWER for all four and does
# not match. Peacock, Turkey and Emu already have the right genomes under
# ref/source_assemblies and are not touched here.
set -eu

REF=/path/to/storage/kleiber_birds_gemcat/ref
ASM=$REF/source_assemblies
TMP=$REF/_fetch_tmp
mkdir -p "$ASM" "$TMP"

fetch() {
  sp=$1 acc=$2
  out=$ASM/$sp
  mkdir -p "$out"
  if [ -s "$out/genome.fa" ] && [ -s "$out/protein.faa" ] && [ -s "$out/genomic.gff" ]; then
    echo "[$sp] уже на месте"
    return
  fi
  echo "[$sp] $acc: качаю геном, белки и GFF"
  zip=$TMP/$sp.zip
  [ -s "$zip" ] || datasets download genome accession "$acc" \
      --include genome,protein,gff3 --filename "$zip" --no-progressbar
  rm -rf "$TMP/$sp"
  unzip -q -o "$zip" -d "$TMP/$sp"
  d=$TMP/$sp/ncbi_dataset/data/$acc
  cat "$d"/*_genomic.fna > "$out/genome.fa"
  cp "$d/protein.faa" "$out/protein.faa"
  cp "$d/genomic.gff" "$out/genomic.gff"
  echo "$acc" > "$out/accessions.txt"
  printf '[%s] геном %s, белков %s\n' "$sp" \
      "$(grep -c '^>' "$out/genome.fa")" "$(grep -c '^>' "$out/protein.faa")"
}

fetch Chicken GCF_000002315.6   # GRCg6a,  пики на NC_006088.5
fetch Duck    GCF_015476345.1   # ZJU1.0,  пики на NC_051772.1
fetch Quail   GCF_001577835.2   # Cja 2.1, пики на NC_029516.1
fetch Zebra   GCF_008822105.2   # bTaeGut2.pat.W.v2, пики на NC_045000.1

echo "=== проверка: совпадают ли имена контигов с таблицами пиков ==="
for sp in Chicken:NC_006088.5 Duck:NC_051772.1 Quail:NC_029516.1 Zebra:NC_045000.1; do
  s=${sp%%:*}; c=${sp##*:}
  if grep -q "^>$c" "$ASM/$s/genome.fa"; then echo "[$s] $c найден"; else echo "[$s] $c НЕ НАЙДЕН"; fi
done
echo "=== закачка готова ==="
