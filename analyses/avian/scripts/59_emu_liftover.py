#!/usr/bin/env python3
"""Lift the emu CAGE peaks from ZJU1.0 onto bDroNov1 and use the real annotation.

Why. Emu lagged behind the other species after re-annotation (83 of 136 OXPHOS
names against 89-97 in the natively annotated birds). Two reasons, and neither is
the query proteome, which is already the species' own newest RefSeq:

  1. the peaks sit on ZJU1.0, which has no annotation of its own, so gene
     positions had to be projected there from protein alignments;
  2. a protein alignment gives the start of the CODING region, not the
     transcription start. Where the first exon is non-coding, the true TSS can be
     kilobases upstream, and a peak sitting at the real TSS then falls outside the
     window.

bDroNov1 (GCF_036370855.1) is the same species and carries a full RefSeq
annotation with 23 985 named genes and proper transcript starts. So: take each
peak's sequence from ZJU1.0, map it to bDroNov1 with minimap2, and score the
distance to the nearest annotated transcript start there.

99 JABVCD contigs are not served individually by ENA, so 586 of the 14 361 peaks
cannot be lifted; that loss is reported, not hidden.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES, load_geneset, normalize_symbol  # noqa: E402

EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
ZJU = Path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies/Emu/genome.fa")
BD = Path("/path/to/storage/kleiber_birds_gemcat/ref/emu_bdronov")
OUT = RES / "emu_liftover"
WINDOWS = (1000, 2000)


def sh(cmd: list[str], **kw):
    print("  " + " ".join(map(str, cmd))[:160])
    return subprocess.run(cmd, check=True, **kw)


def step_peaks_fasta() -> pd.DataFrame:
    d = pd.read_csv(EXP / "Emu_fib_epi_exp.txt", sep="\t", low_memory=False)
    fib = [c for c in d.columns if c.startswith("norm.C")]
    d["fibv"] = (
        d[fib].apply(pd.to_numeric, errors="coerce").fillna(0.0).where(lambda x: x >= 1.0, 0.0).mean(axis=1)
    )
    bed = OUT / "peaks_zju.bed"
    fa = OUT / "peaks_zju.fa"
    if not fa.exists():
        # ENA headers are 'ENA|CM027945|CM027945.1 ...'; make a plain-name index copy
        plain = OUT / "zju_plain.fa"
        if not plain.exists():
            with ZJU.open() as fi, plain.open("w") as fo:
                for line in fi:
                    fo.write(">" + line[1:].split()[0].split("|")[-1] + "\n" if line[0] == ">" else line)
        sh(["samtools", "faidx", str(plain)])
        have = {l.split("\t")[0] for l in (OUT / "zju_plain.fa.fai").read_text().splitlines()}
        sub = d[d["seqnames"].isin(have)]
        with bed.open("w") as fh:
            for _, r in sub.iterrows():
                fh.write(f"{r['seqnames']}\t{int(r['start'])-1}\t{int(r['end'])}\t{r['id']}\t0\t{r['strand']}\n")
        sh(["bedtools", "getfasta", "-fi", str(plain), "-bed", str(bed), "-nameOnly", "-fo", str(fa)])
        print(f"  пиков в FASTA: {sum(1 for l in fa.open() if l[0] == '>')} из {len(d)}")
    return d


def step_map(threads: int) -> pd.DataFrame:
    paf = OUT / "peaks_to_bdronov.paf"
    if not paf.exists():
        with paf.open("w") as fh:
            sh(["minimap2", "-t", str(threads), "-x", "asm5", "-k", "15", "-w", "5",
                "--secondary=no", "-N", "1", str(BD / "genome.fa"), str(OUT / "peaks_zju.fa")], stdout=fh)
    rows = []
    for line in paf.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        rows.append(dict(peak=f[0], qlen=int(f[1]), strand=f[4], tseq=f[5],
                         tstart=int(f[7]), tend=int(f[8]), matches=int(f[9]), mapq=int(f[11])))
    p = pd.DataFrame(rows)
    if p.empty:
        return p
    p["cov"] = p["matches"] / p["qlen"]
    p = p.sort_values(["peak", "matches"], ascending=[True, False]).drop_duplicates("peak")
    p = p[(p["mapq"] >= 20) & (p["cov"] >= 0.5)]
    p.to_csv(OUT / "lifted_peaks.tsv", sep="\t", index=False)
    return p


def step_tss() -> pd.DataFrame:
    """gene symbol -> transcript starts from the bDroNov1 RefSeq GFF"""
    f = OUT / "bdronov_tss.tsv"
    if f.exists():
        return pd.read_csv(f, sep="\t")
    rows = []
    with (BD / "genomic.gff").open() as fh:
        for line in fh:
            if line[0] == "#":
                continue
            c = line.split("\t")
            if len(c) < 9 or c[2] not in ("mRNA", "transcript"):
                continue
            g = re.search(r"gene=([^;\n]+)", c[8])
            if not g:
                continue
            st = c[6]
            rows.append(dict(gene=normalize_symbol(g.group(1)), tseq=c[0], strand=st,
                             tss=int(c[3]) if st == "+" else int(c[4])))
    t = pd.DataFrame(rows).drop_duplicates()
    t = t[t["gene"] != ""]
    t.to_csv(f, sep="\t", index=False)
    return t


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=64)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    peaks = step_peaks_fasta()
    lifted = step_map(a.threads)
    print(f"перенесено пиков: {len(lifted)} из {len(peaks)} ({len(lifted)/len(peaks):.1%})")
    tss = step_tss()
    print(f"транскриптных стартов в bDroNov1: {len(tss)}, генов {tss['gene'].nunique()}")

    m = lifted.merge(peaks[["id", "fibv"]], left_on="peak", right_on="id", how="left")
    m["mid"] = (m["tstart"] + m["tend"]) / 2.0
    m = m[m["fibv"] > 0]
    by = {k: v for k, v in m.groupby("tseq", observed=True)}

    sets = {"hsa00190": "kegg_00190", "hsa04714": "kegg_04714"}
    res = []
    for nm, key in sets.items():
        want = {normalize_symbol(g) for g in load_geneset(key)} - {""}
        sub = tss[tss["gene"].isin(want)]
        found = {w: 0 for w in WINDOWS}
        for g, grp in sub.groupby("gene", observed=True):
            best = None
            for _, r in grp.iterrows():
                cand = by.get(r["tseq"])
                if cand is None:
                    continue
                dmin = (cand["mid"] - r["tss"]).abs().min()
                best = dmin if best is None else min(best, dmin)
            if best is None:
                continue
            for w in WINDOWS:
                if best <= w:
                    found[w] += 1
        res.append(dict(набор=nm, n_карта=len(want), генов_в_bDroNov1=sub["gene"].nunique(),
                        **{f"с_пиком_{w}": found[w] for w in WINDOWS}))
    t = pd.DataFrame(res)
    t.to_csv(OUT / "emu_after_liftover.tsv", sep="\t", index=False)
    print("\n" + t.to_string(index=False))
    print("\nдля сравнения: у эму по проекции протеома было 83 (hsa00190) и 123 (hsa04714);"
          "\nнативная полоса 89–97 и 149–158")
    print(f"\nwrote {OUT / 'emu_after_liftover.tsv'}")


if __name__ == "__main__":
    main()
