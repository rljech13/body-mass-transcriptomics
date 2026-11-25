#!/bin/bash
# Run 03_promoters_fimo.py for every species once its genome zip is verified; up to 28 in parallel; then step 2b.
cd "$(dirname "$0")/.."
PY=/path/to/user-home/miniconda3/envs/scanpy_env/bin/python
mkdir -p audit/03_species
run_one() {
  s=$1
  while [ ! -f "reference/genomes/$s.zip.ok" ]; do grep -q GENOMES_DONE audit/02_fetch_stdout.log && [ ! -f "reference/genomes/$s.zip.ok" ] && { echo "NO_GENOME $s"; return; }; sleep 20; done
  [ -f "results/motifs/$s.npz" ] && { echo "SKIP $s"; return; }
  [ -f "results/excluded/$s.tsv" ] && { echo "EXCLUDED $s no gene symbols in annotation"; return; }
  $PY code/03_promoters_fimo.py "$s" > "audit/03_species/$s.log" 2>&1 && tail -1 "audit/03_species/$s.log" || echo "FAIL $s"
}
export -f run_one; export PY
tail -n +2 config/genomes.tsv | cut -f1 | xargs -P 28 -I{} bash -c 'run_one {}'
echo "FIMO_ALL_DONE $(date -u +%FT%TZ)"
$PY code/04_enrichment_mass.py > audit/04_stdout.log 2>&1 && echo STEP2_DONE || echo STEP2_FAIL
