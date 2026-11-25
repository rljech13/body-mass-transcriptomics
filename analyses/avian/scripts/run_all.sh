#!/usr/bin/env bash
set -euo pipefail
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libsqlite3.so.0
ROOT=/path/to/user-home/chicken/kleiber_birds_gemcat
cd "$ROOT"
PY=/path/to/user-home/miniconda3/bin/python3

echo "=== A1–A3 ingest + allometry + exploratory ==="
$PY scripts/01_ingest_and_allometry.py

echo "=== A4–A5 ESC/fib + GEMCAT ==="
$PY scripts/02_chicken_esc_fib_gemcat.py

echo "=== A6 CAGE promoter regulation ==="
$PY scripts/03_cage_promoter_regulation.py

echo "=== B1 GWAS CRE inventory ==="
$PY scripts/04_gwas_cre_inventory.py

echo "=== Writing summary note ==="
$PY scripts/05_write_summary.py

echo "ALL DONE"
