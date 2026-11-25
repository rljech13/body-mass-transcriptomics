#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
trap 's=$?; echo "$s" > audit/analysis.exitcode' EXIT
for iteration in $(seq 1 180); do
 if test -f audit/existing_SCENIC.exitcode; then break; fi
 sleep 5
done
test "$(cat audit/existing_SCENIC.exitcode)" = 0
sha256sum code/05_analyze.py code/06_verify.R > audit/analysis_code.before.sha256
cp ../growth_regulators_2026-10-21_v1/inputs/ma2016_full.samples.tsv inputs/species_traits.tsv
/path/to/user-home/miniconda3/envs/scanpy_env/bin/python code/05_analyze.py > logs/analysis.log 2>&1
Rscript code/06_verify.R "$PWD" > logs/independent_R.log 2>&1
sha256sum code/05_analyze.py code/06_verify.R > audit/analysis_code.after.sha256
