#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
trap 's=$?; echo "$s" > audit/existing_SCENIC.exitcode' EXIT
sha256sum code/04_prune_existing.py > audit/existing_code.before.sha256
for seed in 20261023 20261024 20261025; do
 /path/to/user-home/miniconda3/envs/scanpy_env/bin/python code/04_prune_existing.py --seed "$seed" > "logs/existing_SCENIC_${seed}.log" 2>&1
done
sha256sum code/04_prune_existing.py > audit/existing_code.after.sha256
