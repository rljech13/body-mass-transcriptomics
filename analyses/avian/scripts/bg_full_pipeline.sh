#!/usr/bin/env bash
set -euo pipefail
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libsqlite3.so.0
ROOT=/path/to/user-home/chicken/kleiber_birds_gemcat
PY=/path/to/user-home/miniconda3/bin/python3
LOG=$ROOT/logs/bg_full_pipeline.log
mkdir -p "$ROOT/logs" "$ROOT/data/ncbi" "$ROOT/results"
exec > >(tee -a "$LOG") 2>&1
echo "===== START $(date -Is) ====="
cd "$ROOT"

NCBI=$ROOT/data/ncbi

# Redownload if missing or suspiciously small (<200MB; full gene_info is usually larger)
sz=0
if [ -f "$NCBI/gene_info.gz" ]; then
  sz=$(stat -c%s "$NCBI/gene_info.gz")
fi
if [ ! -f "$NCBI/gene_info.gz" ] || [ "$sz" -lt 200000000 ]; then
  echo "[$(date -Is)] Downloading NCBI gene_info.gz (current size=$sz) ..."
  curl -fL --retry 5 --retry-delay 5 -o "$NCBI/gene_info.gz.partial" https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
  mv "$NCBI/gene_info.gz.partial" "$NCBI/gene_info.gz"
else
  echo "[$(date -Is)] gene_info.gz already present size=$sz"
fi
ls -lh "$NCBI/gene_info.gz"

if [ ! -s "$NCBI/gene_info_chicken.tsv" ]; then
  echo "[$(date -Is)] Filtering chicken (9031) ..."
  zcat "$NCBI/gene_info.gz" | awk -F'\t' '$1==9031 {print $2"\t"$3}' > "$NCBI/gene_info_chicken.tsv"
fi
echo "chicken lines: $(wc -l < "$NCBI/gene_info_chicken.tsv")"

if [ ! -s "$NCBI/gene_info_human_symbol.tsv" ]; then
  echo "[$(date -Is)] Filtering human (9606) symbols ..."
  zcat "$NCBI/gene_info.gz" | awk -F'\t' '$1==9606 {print toupper($3)"\t"$2}' > "$NCBI/gene_info_human_symbol.tsv"
fi
echo "human lines: $(wc -l < "$NCBI/gene_info_human_symbol.tsv")"

echo "[$(date -Is)] Building iES1300 / Recon symbol maps ..."
$PY << 'PY'
import cobra, pandas as pd
from pathlib import Path
DATA = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data")
m = cobra.io.read_sbml_model(str(DATA / "ies1300/S1_File/S1_File in sbml format.xml"))
ents = {g.id for g in m.genes}
info = pd.read_csv(DATA / "ncbi/gene_info_chicken.tsv", sep="\t", header=None, names=["chicken_entrez", "symbol"])
info["chicken_entrez"] = info["chicken_entrez"].astype(str)
info["symbol"] = info["symbol"].astype(str).str.upper()
sub = info[info.chicken_entrez.isin(ents)].drop_duplicates("chicken_entrez")
sub.to_csv(DATA / "ies1300_entrez_to_symbol.tsv", sep="\t", index=False)
sym = sub.drop_duplicates("symbol", keep="first")
sym.to_csv(DATA / "symbol_to_chicken_entrez.tsv", sep="\t", index=False)
hum = pd.read_csv(DATA / "ncbi/gene_info_human_symbol.tsv", sep="\t", header=None, names=["symbol", "human_entrez"])
hum["symbol"] = hum["symbol"].astype(str).str.upper()
hum = hum.drop_duplicates("symbol", keep="first")
merged = sym.merge(hum, on="symbol", how="inner")
merged[["symbol", "human_entrez"]].to_csv(DATA / "symbol_to_human_entrez_ies_overlap.tsv", sep="\t", index=False)
print(f"iES1300 genes={len(ents)} with_symbol={len(sub)} human_overlap={len(merged)}")
PY

if [ ! -s results/chicken_esc_vs_fib_gene_tpm.tsv ]; then
  echo "[$(date -Is)] Regenerating ESC/fib table ..."
  $PY << 'PY'
from pathlib import Path
import numpy as np, pandas as pd
ATLAS = Path("/path/to/user-home/chicken/CAGE_atlases/qc_atlas/results/05_gene_tpm_matrix.tsv")
RES = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/results")
FIB = ["C08_S8_R1_001.tpm","C09_S9_R1_001.tpm","C10_S32_R1_001.tpm","C11_S10_R1_001.tpm","C12_S11_R1_001.tpm"]
ESC = ["E05_S13_R1_001.tpm","E06_S14_R1_001.tpm","E07_S15_R1_001.tpm"]
tpm = pd.read_csv(ATLAS, sep="\t", index_col=0)
tpm.index = tpm.index.astype(str).str.upper()
fib = tpm[FIB].mean(axis=1)
esc = tpm[ESC].mean(axis=1)
logfc = np.log2((esc + 1) / (fib + 1))
pd.DataFrame({"fib_tpm_mean": fib, "esc_tpm_mean": esc, "logFC_ESC_over_fib": logfc}).to_csv(
    RES / "chicken_esc_vs_fib_gene_tpm.tsv", sep="\t"
)
print("wrote esc/fib table")
PY
fi

echo "[$(date -Is)] A5 GEMCAT iES1300 primary + Recon3D baseline ..."
$PY scripts/02b_gemcat_ies_primary_recon_baseline.py

echo "[$(date -Is)] A6 CAGE promoters ..."
$PY scripts/03_cage_promoter_regulation.py

echo "[$(date -Is)] B1 GWAS inventory ..."
$PY scripts/04_gwas_cre_inventory.py

echo "[$(date -Is)] Summary ..."
$PY scripts/05_write_summary.py

echo "===== ALL DONE $(date -Is) ====="
echo DONE > "$ROOT/logs/bg_full_pipeline.done"
