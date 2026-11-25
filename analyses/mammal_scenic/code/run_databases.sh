#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
trap 's=$?; echo "$s" > audit/databases.exitcode' EXIT
sha256sum code/03_build_db.py > audit/build_code.before.sha256
/path/to/user-home/miniconda3/envs/scanpy_env/bin/python code/03_build_db.py
sha256sum code/03_build_db.py > audit/build_code.after.sha256
