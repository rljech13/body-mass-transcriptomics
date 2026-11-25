#!/usr/bin/env python3
"""Step 2a: native promoter windows and FIMO motif presence per species (plan_step2.json).
Usage: 03_promoters_fimo.py prep            -> results/02_panel_backgrounds.json
       03_promoters_fimo.py <Species_name>  -> results/promoters/<species>.{fa,tsv}, results/motifs/<species>.npz"""
import gzip, json, re, shutil, subprocess, sys, zipfile
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN / "code"))
plan = json.loads((RUN / "config/plan_step2.json").read_text())
MEME = RUN.parent / "size_program_regulons_2026-10-21_v1/reference/JASPAR2024_one_per_symbol.meme"
FIMO = "/path/to/user-home/miniconda3/envs/meme/bin/fimo"
UP, DOWN = 1000, 200
OUT = RUN / "results"

if sys.argv[1] == "prep":
    import common_rm
    bg, members = {}, {}
    for panel in ["ma2016", "gse105160", "wagner"]:
        d = common_rm.load(panel)
        Z = common_rm.background(d)
        a = common_rm.adm[(common_rm.adm.panel == panel) & (common_rm.adm.status == "admitted")].set_index("label")
        bg[panel] = sorted(Z.index)
        members[panel] = {lab: a.loc[lab, "scientific_name"] for lab in d["species"]}
    (OUT / "02_panel_backgrounds.json").write_text(json.dumps({"background": bg, "species": members}, indent=1))
    print({p: len(v) for p, v in bg.items()})
    sys.exit()

species = sys.argv[1]
info = json.loads((OUT / "02_panel_backgrounds.json").read_text())
need = sorted({g for p, labs in info["species"].items() if species in labs.values() for g in info["background"][p]})
work = RUN / "tmp" / species
work.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(RUN / "reference/genomes" / f"{species}.zip") as z:
    names = z.namelist()
    fna = [n for n in names if n.endswith(".fna")]
    gff = [n for n in names if n.endswith(".gff")]
    assert len(fna) == 1 and len(gff) == 1, (species, fna, gff)
    for n in fna + gff:
        target = work / Path(n).name
        if not target.exists():
            with z.open(n) as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst, 1 << 24)
fna_p, gff_p = work / Path(fna[0]).name, work / Path(gff[0]).name

attr_re = re.compile(r"([^=;]+)=([^;]*)")
genes, mrna, exlen = {}, {}, defaultdict(int)
with open(gff_p) as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        a = dict(attr_re.findall(c[8]))
        if c[2] == "gene" and a.get("gene_biotype", "protein_coding") == "protein_coding":
            nm = a.get("Name") or a.get("gene")
            if nm:
                genes[a["ID"]] = nm.upper()
        elif c[2] == "mRNA" and "Parent" in a:
            mrna[a["ID"]] = (a["Parent"], c[0], int(c[3]), int(c[4]), c[6])
        elif c[2] == "exon" and "Parent" in a:
            exlen[a["Parent"]] += int(c[4]) - int(c[3]) + 1
count = defaultdict(int)
for nm in genes.values():
    count[nm] += 1
best = {}
for tid, (gid, chrom, s, e, strand) in mrna.items():
    nm = genes.get(gid)
    if nm is None or count[nm] != 1:
        continue
    L = exlen.get(tid, e - s + 1)
    if nm not in best or L > best[nm][0]:
        best[nm] = (L, chrom, s, e, strand)

from Bio import SeqIO
from Bio.Seq import Seq
by_chrom = defaultdict(list)
for nm in need:
    if nm in best:
        by_chrom[best[nm][1]].append(nm)
rows = []
pdir, mdir = OUT / "promoters", OUT / "motifs"
pdir.mkdir(exist_ok=True); mdir.mkdir(exist_ok=True)
fa = pdir / f"{species}.fa"
with open(fa, "w") as out:
  for rec in SeqIO.parse(str(fna_p), "fasta"):
    if rec.id not in by_chrom:
        continue
    chrom_seq = str(rec.seq)
    for nm in by_chrom[rec.id]:
        _, chrom, s, e, strand = best[nm]
        if strand == "+":
            lo, hi = s - 1 - UP, s - 1 + DOWN
        else:
            lo, hi = e - DOWN, e + UP
        if lo < 0 or hi > len(chrom_seq):
            continue
        seq = chrom_seq[lo:hi]
        if strand == "-":
            seq = str(Seq(seq).reverse_complement())
        if len(seq) != UP + DOWN or seq.upper().count("N") > 0.2 * len(seq):
            continue
        up = seq.upper()
        gc = (up.count("G") + up.count("C")) / max(1, len(up) - up.count("N"))
        rows.append(dict(gene=nm, chrom=chrom, strand=strand, window_start0=lo, window_end=hi, gc=gc))
        out.write(f">{nm}\n{seq}\n")
pd.DataFrame(rows).to_csv(pdir / f"{species}.tsv", sep="\t", index=False)

motif_ids = [l.split()[1] + "|" + l.split()[2] for l in MEME.read_text().splitlines() if l.startswith("MOTIF")]
mpos = {m.split("|")[0]: i for i, m in enumerate(motif_ids)}
gpos = {r["gene"]: i for i, r in enumerate(rows)}
P = np.zeros((len(rows), len(motif_ids)), dtype=bool)
lines = fa.read_text().splitlines()
chunks = []
for k in range(4):
    cf = work / f"chunk{k}.fa"
    cf.write_text("\n".join(lines[2 * i + j] for i in range(k, len(lines) // 2, 4) for j in (0, 1)) + "\n")
    chunks.append((cf, open(work / f"chunk{k}.fimo", "w")))
procs = [subprocess.Popen([FIMO, "--text", "--thresh", "1e-4", "--verbosity", "1", str(MEME), str(cf)], stdout=fh) for cf, fh in chunks]
assert all(p.wait() == 0 for p in procs)
for cf, fh in chunks:
    fh.close()
    with open(work / (cf.stem + ".fimo")) as f:
        for line in f:
            if line.startswith("motif_id") or line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) > 3 and c[2] in gpos:
                P[gpos[c[2]], mpos[c[0]]] = True
np.savez_compressed(mdir / f"{species}.npz", presence=P, genes=np.array([r["gene"] for r in rows]),
                    motifs=np.array(motif_ids), gc=np.array([r["gc"] for r in rows]))
shutil.rmtree(work)
print(species, "needed", len(need), "unique_pc_genes", len(best), "windows", len(rows), "positive_pairs", int(P.sum()))
