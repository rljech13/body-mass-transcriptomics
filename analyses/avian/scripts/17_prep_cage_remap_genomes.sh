#!/usr/bin/env bash
# Prepare better genomes (+MT) for CAGE remap of GSE213253 fib libraries.
# Raw CAGE FASTQs are in SRA (library_selection=CAGE), ~26 GB total — not local except chicken.
set -euo pipefail

ROOT=/path/to/storage/kleiber_birds_gemcat/ref/genomes_with_mt
MANIFEST=/path/to/user-home/chicken/kleiber_birds_gemcat/results/cage_remap_mt
mkdir -p "$ROOT" "$MANIFEST/logs"

# Preferred assemblies (chromosome-level / recent) + dedicated MT accessions
# Format: species accession mt_accession
PLAN_FILE="$MANIFEST/genome_accessions.tsv"
cat > "$PLAN_FILE" <<'EOF'
species	assembly	mt
Emu	GCF_036370855.1	NC_002784.1
Peacock	GCA_965225075.1	NC_024533.1
Turkey	GCF_905368555.1	NC_010195.2
Quail	GCF_054131305.1	NC_003408.1
Duck	GCA_059339885.1	NC_009684.1
Zebra	GCA_051427915.1	NC_007897.1
Chicken	LOCAL	NC_001323.1
EOF

tail -n +2 "$PLAN_FILE" | while IFS=$'\t' read -r sp acc mt; do
  [[ -z "${sp:-}" ]] && continue
  echo "=== $sp $acc MT=$mt ==="
  out="$ROOT/$sp"
  mkdir -p "$out"
  if [[ "$acc" == "LOCAL" ]]; then
    ln -sfn /path/to/user-home/chicken/CAGE_atlases/qc_reviewer/bpnet/galGal6.main.fa "$out/genome.fa"
    [[ -f /path/to/user-home/chicken/CAGE_atlases/qc_reviewer/bpnet/galGal6.main.fa.fai ]] && \
      ln -sfn /path/to/user-home/chicken/CAGE_atlases/qc_reviewer/bpnet/galGal6.main.fa.fai "$out/genome.fa.fai"
    [[ -f "$out/genome.fa.fai" ]] || samtools faidx "$out/genome.fa"
    echo "$sp ready (local galGal6 with chrM)" | tee -a "$MANIFEST/logs/genome_prep.log"
    continue
  fi
  if [[ ! -s "$out/genome.fa" ]]; then
    echo "Downloading $acc ..."
    datasets download genome accession "$acc" --include genome --filename "$out/${acc}.zip" 2>&1 | tee -a "$MANIFEST/logs/genome_prep.log"
    unzip -o "$out/${acc}.zip" -d "$out/unzip" >/dev/null
    fa=$(find "$out/unzip" -name '*.fna' -o -name '*.fa' | head -1)
    cp "$fa" "$out/genome.nuclear.fa"
    # append MT if missing
    if ! grep -qE '^>(chrM|MT|mitochondrion)' "$out/genome.nuclear.fa"; then
      curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${mt}&rettype=fasta&retmode=text" \
        | sed "1s/^>.*/>chrM ${mt}/" > "$out/chrM.fa"
      cat "$out/genome.nuclear.fa" "$out/chrM.fa" > "$out/genome.fa"
      echo "Appended $mt as chrM" | tee -a "$MANIFEST/logs/genome_prep.log"
    else
      cp "$out/genome.nuclear.fa" "$out/genome.fa"
    fi
    samtools faidx "$out/genome.fa"
  else
    echo "$sp already has genome.fa"
  fi
  # sanity: chrM present?
  if cut -f1 "$out/genome.fa.fai" | grep -qiE 'chrM|^MT$|mitochondrion'; then
    echo "OK $sp has MT contig" | tee -a "$MANIFEST/logs/genome_prep.log"
  else
    echo "WARN $sp missing MT in fai" | tee -a "$MANIFEST/logs/genome_prep.log"
  fi
done

echo "Done genomes → $ROOT"
echo "CAGE SRA manifest: $MANIFEST/gse213253_cage_sra_manifest.tsv"
