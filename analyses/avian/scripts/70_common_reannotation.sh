#!/usr/bin/env bash
# One proteome projected onto all seven genomes, same query everywhere.
#
# Script 58 gave each of the three Augustus species a DIFFERENT query — turkey got
# its own proteome, emu its own, peacock the quail one — and gave the other four
# nothing at all. That is why locus counts run 10 000 for the natively annotated
# species against 16-18 000 for the re-annotated ones, and why only 6567 genes
# have a locus in all seven. A ceiling set by a species-specific procedure is not
# something the seven-way intersection can be blamed for.
#
# Here the query is a single combined set: one protein per gene from six bird
# proteomes spanning the panel's phylogeny — Galliformes (chicken, quail, turkey),
# Anseriformes (duck), Passeriformes (zebra finch) and Palaeognathae (emu). Every
# genome sees that same set, so each locus has the same chance of being found
# regardless of which species' genome it sits in. Peacock has no RefSeq proteome
# of its own; under this scheme it no longer needs one.
#
# Genomes are the ones the CAGE peaks were called on, NOT what sits in
# ref/genomes_with_mt (newer assemblies for four of the seven).
set -eu

REF=/path/to/storage/kleiber_birds_gemcat/ref
ASM=$REF/source_assemblies
P=$REF/proteomes
THREADS=${THREADS:-96}
QUERY=$P/common_query_named.faa

# protein_id -> gene symbol from the NCBI GFF, longest protein per gene.
#
# LOC-prefixed symbols are dropped. They are per-species placeholders for unnamed
# genes and do not agree between species even for orthologues, so they inflate the
# query from 16 779 real symbols to 38 394 and plant ~21 600 phantom loci in every
# genome. With peaks assigned to the nearest locus those phantoms steal promoters:
# in the first run ATP5MF and UQCRC2 went from detected in all seven species to
# detected in one. Nothing downstream uses unnamed genes anyway.
one_per_gene() {
  faa=$1 gff=$2 out=$3 tagsp=$4
  [ -s "$out" ] && { echo "[$tagsp] onePerGene уже есть: $(grep -c '^>' "$out")"; return; }
  python3 - "$faa" "$gff" "$out" "$tagsp" <<'PY'
import re, sys
faa, gff, out, sp = sys.argv[1:5]
pid2sym = {}
with open(gff) as fh:
    for line in fh:
        if line[0] == '#':
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        m = re.search(r'protein_id=([^;\n]+)', f[8])
        g = re.search(r'gene=([^;\n]+)', f[8])
        if m and g:
            s = g.group(1).upper()
            if re.match(r'^LOC\d', s):
                continue
            pid2sym.setdefault(m.group(1), s)
best = {}
pid, buf = None, []
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
            pid, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
flush()
with open(out, 'w') as fh:
    for sym, (pid, seq) in best.items():
        fh.write(f'>{sp}_{pid}##{sym}\n')
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + '\n')
print(f'[{sp}] генов {len(best)} из {len(pid2sym)} белков с символом')
PY
}

echo "=== 1. один белок на ген из шести протеомов ==="
one_per_gene "$ASM/Chicken/protein.faa" "$ASM/Chicken/genomic.gff" "$P/chicken_named.faa" chicken
one_per_gene "$ASM/Duck/protein.faa"    "$ASM/Duck/genomic.gff"    "$P/duck_named.faa"    duck
one_per_gene "$ASM/Zebra/protein.faa"   "$ASM/Zebra/genomic.gff"   "$P/zebra_named.faa"   zebra
one_per_gene "$ASM/Quail/protein.faa"   "$ASM/Quail/genomic.gff"   "$P/quail_named.faa"   quail
one_per_gene "$P/turkey.faa"            "$P/turkey.gff"            "$P/turkey_named.faa"  turkey
one_per_gene "$P/emu.faa"               "$P/emu.gff"               "$P/emu_named.faa"     emu

echo "=== 2. общий запрос ==="
if [ ! -s "$QUERY" ]; then
  cat "$P/chicken_named.faa" "$P/duck_named.faa" "$P/zebra_named.faa" \
      "$P/quail_named.faa" "$P/turkey_named.faa" "$P/emu_named.faa" > "$QUERY.part"
  mv "$QUERY.part" "$QUERY"
fi
echo "общий запрос: $(grep -c '^>' "$QUERY") белков, символов $(grep '^>' "$QUERY" | sed 's/.*##//' | sort -u | wc -l)"

echo "=== 3. miniprot на семь геномов ==="
align() {
  sp=$1
  genome=$ASM/$sp/genome.fa
  [ "$sp" = "Peacock" ] && genome=$ASM/$sp/genome_localnames.fa
  gff=$ASM/$sp/reanno_named.gff
  [ -s "$gff" ] && { echo "[$sp] reanno_named.gff уже есть: $(awk '$3=="mRNA"' "$gff" | wc -l) моделей"; return; }
  echo "[$sp] старт на $(basename "$genome")"
  /usr/bin/time -f "[$sp] %E прошло, %M КБ" \
    miniprot -t "$THREADS" -Iu --gff "$genome" "$QUERY" \
    > "$gff.part" 2> "$ASM/$sp/reanno_named.log"
  mv "$gff.part" "$gff"
  echo "[$sp] моделей $(awk '$3=="mRNA"' "$gff" | wc -l), символов $(awk '$3=="mRNA"' "$gff" | grep -o 'Target=[^ ]*##[^ ]*' | sed 's/.*##//' | sort -u | wc -l)"
}

for sp in Chicken Duck Quail Zebra Peacock Turkey Emu; do align "$sp"; done
echo "=== переаннотация готова ==="
