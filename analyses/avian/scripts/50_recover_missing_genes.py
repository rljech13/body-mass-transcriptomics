#!/usr/bin/env python3
"""Find the unannotated copies of the missing OXPHOS / thermogenesis genes.

Scope set by the user: only hsa00190 and hsa04714, no whole-genome re-annotation.

Idea. A gene is "missing" in Peacock / Turkey / Emu because the Augustus gene set
delivered with the CAGE tables has no model for it, so no peak could be labelled
with its symbol. The CAGE peaks themselves are there — those three species have
as many peaks as the others (14361 in Emu against 13503 in Chicken), they are
just classified 'not promoter' or tied to an unnamed model. So we place the gene
ourselves, in the coordinate system the peaks already live in, and ask whether a
peak sits at its 5' end.

Procedure per species:
  1. chicken protein (longest isoform) of each target gene
  2. miniprot onto the very assembly the peaks were called on -> gene models
  3. 5' end of the best model = start if '+' else end
  4. is there a peak on the same sequence and strand within W bp of that 5' end
  5. windows 500 / 1000 / 2000 bp reported separately

Negative control. 15 genes that are absent from the CHICKEN peak table too, not
because of annotation but because they are not transcribed in fibroblasts or
epiblast (UCP1, PPARGC1A, ADRB3, FGF21, PRDM16, NPPA, NPPB, CPT1B, COX4I2,
COX6A2, COX7A1, ATP4A, ATP6V1B1, ADCY1, PLIN1). miniprot must find their loci,
and there must be NO peak there. If controls come back "recovered", the window
is too permissive and the whole recovery is noise.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES  # noqa: E402

EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
REF = Path("/path/to/storage/kleiber_birds_gemcat/ref")
ASM = REF / "source_assemblies"
PROT = REF / "chicken_prot" / "query_proteins_all.faa"
OUT = RES / "gene_recovery"
WINDOWS = (500, 1000, 2000)
CONTROLS = {
    "UCP1", "PPARGC1A", "ADRB3", "FGF21", "PRDM16", "NPPA", "NPPB", "CPT1B",
    "COX4I2", "COX6A2", "COX7A1", "ATP4A", "ATP6V1B1", "ADCY1", "PLIN1",
}
HDR = re.compile(r"^>(\S+)\s+(\S+)\s+\[organism=")


def longest_per_symbol(faa: Path) -> dict[str, tuple[str, str]]:
    """symbol -> (protein_id, sequence), longest isoform."""
    best: dict[str, tuple[str, str]] = {}
    pid = sym = None
    buf: list[str] = []

    def flush() -> None:
        if not sym or not buf:
            return
        seq = "".join(buf)
        if sym not in best or len(seq) > len(best[sym][1]):
            best[sym] = (pid, seq)

    for line in faa.read_text().splitlines():
        if line.startswith(">"):
            flush()
            m = HDR.match(line)
            pid, sym, buf = (m.group(1), m.group(2), []) if m else (None, None, [])
        elif sym:
            buf.append(line.strip())
    flush()
    return best


def write_queries(symbols: list[str], prots: dict[str, tuple[str, str]], path: Path) -> dict[str, str]:
    pid2sym = {}
    with path.open("w") as fh:
        for s in symbols:
            if s not in prots:
                continue
            pid, seq = prots[s]
            pid2sym[pid] = s
            fh.write(f">{pid}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    return pid2sym


def run_miniprot(genome: Path, queries: Path, gff: Path, threads: int) -> None:
    if gff.exists() and gff.stat().st_size > 0:
        print(f"  miniprot: {gff.name} уже есть, пропускаю")
        return
    cmd = ["miniprot", "-t", str(threads), "-Iu", "--gff", str(genome), str(queries)]
    print("  " + " ".join(cmd))
    with gff.open("w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


def clean_seqname(s: str) -> str:
    """ENA serves FASTA as 'ENA|HG999680|HG999680.1 ...'; peak tables use 'HG999680.1'."""
    s = s.split()[0]
    if s.startswith("ENA|"):
        s = s.split("|")[-1]
    return s


def parse_gff(gff: Path, pid2sym: dict[str, str]) -> pd.DataFrame:
    rows = []
    for line in gff.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
        tgt = attrs.get("Target", "").split()[0]
        sym = pid2sym.get(tgt)
        if not sym:
            continue
        rows.append(
            dict(
                gene=sym,
                seqname=clean_seqname(f[0]),
                start=int(f[3]),
                end=int(f[4]),
                strand=f[6],
                score=float(f[5]) if f[5] not in (".", "") else 0.0,
                identity=float(attrs.get("Identity", 0) or 0),
            )
        )
    return pd.DataFrame(rows)


def peaks_of(species: str) -> pd.DataFrame:
    d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    fib = [c for c in d.columns if c.startswith("norm.C")]
    epi = [c for c in d.columns if c.startswith("norm.E")]
    d["fib_mean"] = d[fib].mean(axis=1) if fib else 0.0
    d["epi_mean"] = d[epi].mean(axis=1) if epi else 0.0
    keep = ["id", "seqnames", "start", "end", "strand", "tss", "class", "fib_mean", "epi_mean"]
    return d[[c for c in keep if c in d.columns]].rename(columns={"seqnames": "seqname"})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="+", default=["Turkey", "Emu", "Peacock"])
    ap.add_argument("--threads", type=int, default=32)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    prots = longest_per_symbol(PROT)
    print(f"белков курицы с символом: {len(prots)}")

    tg = pd.read_csv(RES / "missing_gene_targets" / "recovery_targets.tsv", sep="\t")
    targets = {r["species"]: r["genes"].split(";") for _, r in tg.iterrows()}

    all_rows = []
    for sp in args.species:
        genome = ASM / sp / ("genome_localnames.fa" if sp == "Peacock" else "genome.fa")
        if not genome.exists():
            print(f"\n[{sp}] нет {genome}, пропускаю")
            continue
        print(f"\n[{sp}] цели {len(targets[sp])} + контроли {len(CONTROLS)}")
        syms = sorted(set(targets[sp]) | CONTROLS)
        q = OUT / f"{sp}_queries.faa"
        pid2sym = write_queries(syms, prots, q)
        missing_prot = [s for s in syms if s not in prots]
        if missing_prot:
            print(f"  нет белка курицы для: {', '.join(missing_prot)}")
        gff = OUT / f"{sp}_miniprot.gff"
        run_miniprot(genome, q, gff, args.threads)
        hits = parse_gff(gff, pid2sym)
        if hits.empty:
            print("  miniprot не дал ни одного попадания")
            continue
        best = hits.sort_values("score", ascending=False).drop_duplicates("gene")
        print(f"  локус найден у {len(best)} из {len(syms)} генов")

        pk = peaks_of(sp)
        for _, h in best.iterrows():
            five = h["start"] if h["strand"] == "+" else h["end"]
            same = pk[(pk["seqname"] == h["seqname"]) & (pk["strand"] == h["strand"])].copy()
            same["dist"] = (same["start"] + same["end"]) / 2 - five
            same["absdist"] = same["dist"].abs()
            row = dict(
                species=sp,
                gene=h["gene"],
                is_control=h["gene"] in CONTROLS,
                seqname=h["seqname"],
                five_prime=int(five),
                strand=h["strand"],
                miniprot_identity=round(h["identity"], 3),
                miniprot_score=h["score"],
            )
            for w in WINDOWS:
                near = same[same["absdist"] <= w]
                row[f"n_peaks_{w}"] = len(near)
                if w == max(WINDOWS) and len(near):
                    b = near.sort_values("absdist").iloc[0]
                    row.update(
                        best_peak=b["id"],
                        best_dist=int(b["dist"]),
                        best_class=b["class"],
                        best_fib=round(float(b["fib_mean"]), 2),
                        best_epi=round(float(b["epi_mean"]), 2),
                    )
            all_rows.append(row)

    if not all_rows:
        print("\nнечего сводить")
        return
    tab = pd.DataFrame(all_rows)
    tab.to_csv(OUT / "recovery_hits.tsv", sep="\t", index=False)

    print("\n" + "=" * 96)
    print("ИТОГ: сколько генов получают пик на своём 5'-конце")
    print("=" * 96)
    for sp, sub in tab.groupby("species"):
        t = sub[~sub["is_control"]]
        c = sub[sub["is_control"]]
        print(f"\n{sp}: локусов {len(sub)} (целей {len(t)}, контролей {len(c)})")
        for w in WINDOWS:
            col = f"n_peaks_{w}"
            print(
                f"  ±{w:>5} bp: целей с пиком {int((t[col] > 0).sum()):>3}/{len(t):<3} "
                f"({(t[col] > 0).mean():.0%})   КОНТРОЛЕЙ с пиком {int((c[col] > 0).sum()):>2}/{len(c)} "
                f"({(c[col] > 0).mean():.0%})"
            )
    print(f"\nwrote {OUT / 'recovery_hits.tsv'}")


if __name__ == "__main__":
    main()
