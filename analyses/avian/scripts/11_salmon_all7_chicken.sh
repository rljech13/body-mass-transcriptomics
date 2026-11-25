#!/usr/bin/env bash
# Quantify GSE213253 fibroblast RNA-seq (7 spp × 1F+1M) against chicken cDNA symbols.
# Handles SINGLE and PAIRED libraries (from pysradb). Skips existing quant.sf.
set -euo pipefail
WORKDIR=/path/to/storage/kleiber_birds_gemcat
HOMEPROJ=/path/to/user-home/chicken/kleiber_birds_gemcat
IDX=$WORKDIR/salmon_index/chicken_cdna_symbol
LOG=$WORKDIR/logs/salmon_all7.log
RUNS=$WORKDIR/runs_all7_high_1F1M.tsv
mkdir -p "$WORKDIR"/{fastq,quant,logs}
exec > >(tee -a "$LOG") 2>&1
echo "===== START salmon_all7 $(date -Is) ====="

if [ ! -f "$IDX/info.json" ] && [ ! -f "$IDX/versionInfo.json" ]; then
  echo "ERROR: missing salmon index $IDX"; exit 1
fi

# Refresh run table with correct ENA URLs + layout from pysradb
/path/to/user-home/miniconda3/bin/python3 - << 'PY'
import pandas as pd
from pathlib import Path
W = Path("/path/to/storage/kleiber_birds_gemcat")
org2sp = {
    "Anas platyrhynchos": "Duck",
    "Coturnix japonica": "Quail",
    "Dromaius novaehollandiae": "Emu",
    "Gallus gallus": "Chicken",
    "Meleagris gallopavo": "Turkey",
    "Pavo cristatus": "Peacock",
    "Taeniopygia guttata": "ZebraFinch",
}
base = pd.read_csv(W / "runs_rnaseq_high_1F1M.tsv", sep="\t")
full = pd.read_csv(W / "pysradb_full.tsv", sep="\t", low_memory=False)
cols = ["run_accession", "library_layout", "ena_fastq_http_1", "ena_fastq_http_2", "total_size"]
df = base.merge(full[cols], on="run_accession", how="left", suffixes=("", "_y"))
if "total_size_y" in df.columns:
    df["total_size"] = df["total_size"].fillna(df["total_size_y"])
df["species"] = df["organism_name"].map(org2sp)
df["quant_name"] = df["species"] + "_" + df["sex"].str.lower()
df.to_csv(W / "runs_all7_high_1F1M.tsv", sep="\t", index=False)
print(df[["run_accession", "quant_name", "library_layout", "total_size"]].sort_values("total_size").to_string(index=False))
PY

download_and_quant() {
  local SRR=$1
  local NAME=$2
  local LAYOUT=$3
  local ENA1=$4
  local ENA2=${5:-}
  local QDIR=$WORKDIR/quant/$NAME
  if [ -f "$QDIR/quant.sf" ]; then
    echo "[$(date -Is)] SKIP quant $NAME (exists)"
    return 0
  fi
  echo "[$(date -Is)] === $NAME ($SRR, $LAYOUT) ==="
  cd "$WORKDIR/fastq"
  mkdir -p "$QDIR"

  if [ "$LAYOUT" = "SINGLE" ]; then
    local R="${SRR}.fastq.gz"
    if [ ! -s "$R" ]; then
      echo "[$(date -Is)] curl ENA SE $SRR ..."
      if [ -z "$ENA1" ] || [ "$ENA1" = "-" ] || [ "$ENA1" = "nan" ]; then
        false
      else
        curl -fL --retry 5 --retry-delay 10 -o "$R" "$ENA1"
      fi || {
        echo "[$(date -Is)] ENA SE failed → prefetch"
        rm -f "$R"
        prefetch "$SRR" --max-size 80G
        fasterq-dump -e 8 -O "$WORKDIR/fastq" "$SRR"
        if [ -f "${SRR}.fastq" ]; then pigz -f "${SRR}.fastq"; fi
        # fasterq may emit _1 even for SE in some cases
        if [ ! -s "$R" ] && [ -f "${SRR}_1.fastq" ]; then pigz -f "${SRR}_1.fastq"; mv "${SRR}_1.fastq.gz" "$R"; fi
      }
    fi
    [ -s "$R" ] || { echo "ERROR: no SE FASTQ $SRR"; return 1; }
    echo "[$(date -Is)] salmon quant SE $NAME ..."
    salmon quant -i "$IDX" -l A -r "$R" -p 8 -o "$QDIR" --validateMappings
    rm -f "$R"
  else
    local R1="${SRR}_1.fastq.gz"
    local R2="${SRR}_2.fastq.gz"
    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
      rm -f "$R1" "$R2"
      echo "[$(date -Is)] curl ENA PE $SRR ..."
      if curl -fL --retry 5 --retry-delay 10 -o "$R1" "$ENA1"; then
        curl -fL --retry 5 --retry-delay 10 -o "$R2" "$ENA2"
      else
        echo "[$(date -Is)] ENA PE failed → prefetch"
        rm -f "$R1" "$R2"
        prefetch "$SRR" --max-size 80G
        fasterq-dump --split-files -e 8 -O "$WORKDIR/fastq" "$SRR"
        pigz -f "${SRR}_1.fastq" "${SRR}_2.fastq"
      fi
    fi
    [ -s "$R1" ] && [ -s "$R2" ] || { echo "ERROR: no PE FASTQ $SRR"; return 1; }
    echo "[$(date -Is)] salmon quant PE $NAME ..."
    salmon quant -i "$IDX" -l A -1 "$R1" -2 "$R2" -p 8 -o "$QDIR" --validateMappings
    rm -f "$R1" "$R2"
  fi
  rm -rf "$WORKDIR/fastq/$SRR" 2>/dev/null || true
  echo "[$(date -Is)] done $NAME"
}

# Smaller first
while IFS=$'\t' read -r SRR NAME LAYOUT ENA1 ENA2 SIZE; do
  download_and_quant "$SRR" "$NAME" "$LAYOUT" "$ENA1" "$ENA2"
done < <(/path/to/user-home/miniconda3/bin/python3 - << 'PY'
import pandas as pd
df=pd.read_csv('/path/to/storage/kleiber_birds_gemcat/runs_all7_high_1F1M.tsv', sep='\t')
df=df.sort_values('total_size')
for _,r in df.iterrows():
    e1='' if pd.isna(r.ena_fastq_http_1) else str(r.ena_fastq_http_1)
    e2='' if pd.isna(r.ena_fastq_http_2) else str(r.ena_fastq_http_2)
    print(f"{r.run_accession}\t{r.quant_name}\t{r.library_layout}\t{e1}\t{e2}\t{r.total_size}")
PY
)

echo "[$(date -Is)] Aggregate ..."
/path/to/user-home/miniconda3/bin/python3 "$HOMEPROJ/scripts/11b_aggregate_salmon_all7.py"
echo "===== ALL DONE salmon_all7 $(date -Is) ====="
echo DONE > "$WORKDIR/logs/salmon_all7.done"
