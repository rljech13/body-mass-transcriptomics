#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
trap 's=$?; echo "$s" > audit/networks.exitcode' EXIT
sha256sum code/02_grn.py code/04_prune.py > audit/network_code.before.sha256
py=/path/to/user-home/miniconda3/envs/scanpy_env/bin/python
for seed in 20261023 20261024 20261025; do
 "$py" code/02_grn.py --seed "$seed" > "logs/GRN_${seed}.log" 2>&1
done
for iteration in $(seq 1 360); do
 if test -f audit/databases.exitcode; then break; fi
 sleep 10
done
test "$(cat audit/databases.exitcode)" = 0
for seed in 20261023 20261024 20261025; do
 "$py" code/04_prune.py --seed "$seed" > "logs/SCENIC_${seed}.log" 2>&1
done
sha256sum code/02_grn.py code/04_prune.py > audit/network_code.after.sha256
