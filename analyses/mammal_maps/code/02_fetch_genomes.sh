#!/bin/bash
# Download genome FASTA + GFF3 for each species in config/genomes.tsv from NCBI Datasets; a species is DONE when unzip -t passes.
cd "$(dirname "$0")/.."
mkdir -p reference/genomes
tail -n +2 config/genomes.tsv | while IFS=$'\t' read -r species accession note; do
  z="reference/genomes/${species}.zip"
  if [ -f "$z.ok" ]; then echo "SKIP $species"; continue; fi
  for i in 1 2 3 4 5; do
    curl -sS -L -m 7200 -o "$z" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${accession}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF"
    if unzip -tq "$z" >/dev/null 2>&1; then touch "$z.ok"; echo "DONE $(date -u +%FT%TZ) $species $accession $(stat -c %s "$z")"; break; fi
    echo "RETRY $i $species"; sleep 30
  done
done
(cd reference/genomes && sha256sum *.zip > ../../audit/genomes.sha256)
echo GENOMES_DONE
