#!/usr/bin/env bash
# Full protein-level re-annotation of the three Augustus-annotated assemblies.
#
# Instead of hunting single genes, project a whole proteome onto the very assembly
# the CAGE peaks were called on. Best query per species:
#   Turkey  -> its own RefSeq proteome GCF_905368555.1 (same assembly lineage)
#   Emu     -> its own RefSeq proteome GCF_036370855.1 (same species, ZJU1.0 vs bDroNov1)
#   Peacock -> quail proteome GCF_001577835.2 (no Pavo RefSeq exists; closest Galliform)
#
# One protein per gene (longest isoform) to keep miniprot at a sane size.
set -eu

REF=/path/to/storage/kleiber_birds_gemcat/ref
P=$REF/proteomes
ASM=$REF/source_assemblies
THREADS=${THREADS:-64}

# protein_id -> gene symbol, from the NCBI GFF; keep the longest protein per gene
build_queries() {
  src=$1 gff=$2 out=$3
  [ -s "$out" ] && { echo "[$(basename "$out")] уже есть"; return; }
  python3 - "$src" "$gff" "$out" <<'PY'
import re, sys
faa, gff, out = sys.argv[1:4]
pid2sym = {}
with open(gff) as fh:
    for line in fh:
        if line[0] == '#':
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = f[8]
        m = re.search(r'protein_id=([^;\n]+)', a)
        g = re.search(r'gene=([^;\n]+)', a)
        if m and g:
            pid2sym.setdefault(m.group(1), g.group(1).upper())
best = {}
pid = None; buf = []
def flush():
    if pid is None or not buf:
        return
    s = pid2sym.get(pid)
    if not s:
        return
    seq = ''.join(buf)
    if s not in best or len(seq) > len(best[s][1]):
        best[s] = (pid, seq)
with open(faa) as fh:
    for line in fh:
        if line[0] == '>':
            flush()
            pid = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
flush()
with open(out, 'w') as fh:
    for sym, (pid, seq) in best.items():
        fh.write(f'>{pid}##{sym}\n')
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + '\n')
print(f'{out}: генов {len(best)} из {len(pid2sym)} белков с символом')
PY
}

build_queries "$P/turkey.faa" "$P/turkey.gff" "$P/turkey_onePerGene.faa"
build_queries "$P/emu.faa"    "$P/emu.gff"    "$P/emu_onePerGene.faa"
build_queries "$P/quail.faa"  "$P/quail.gff"  "$P/quail_onePerGene.faa"

align() {
  sp=$1 q=$2
  genome=$ASM/$sp/genome.fa
  [ "$sp" = "Peacock" ] && genome=$ASM/$sp/genome_localnames.fa
  gff=$ASM/$sp/reanno.gff
  [ -s "$gff" ] && { echo "[$sp] reanno.gff уже есть"; return; }
  echo "[$sp] miniprot: $(grep -c '^>' "$q") белков на $(basename "$genome")"
  /usr/bin/time -f "[$sp] %E прошло, %M КБ памяти" \
    miniprot -t "$THREADS" -Iu --gff "$genome" "$q" > "$gff.part" 2> "$ASM/$sp/reanno.log"
  mv "$gff.part" "$gff"
  echo "[$sp] моделей: $(awk '$3=="mRNA"' "$gff" | wc -l)"
}

align Turkey  "$P/turkey_onePerGene.faa"
align Emu     "$P/emu_onePerGene.faa"
align Peacock "$P/quail_onePerGene.faa"
echo "=== переаннотация готова ==="
