#!/usr/bin/env python3
"""Second recovery round, driven by thermogenesis.

Round one queried only the longest CHICKEN protein per symbol and accepted a peak
within 1000 bp. It left, per Augustus species, 53 / 33 / 50 thermogenesis genes
and 29 / 17 / 23 OXPHOS genes that at least three native birds do detect — so
those are annotation gaps, not biology.

Three changes here:
  1. target = detected in >= 3 of the four natively annotated birds, missing in
     this species. Chicken alone is no longer the reference, so genes chicken
     happens to miss are still hunted.
  2. queries = proteins of the same symbol from four birds (chicken, quail,
     zebra finch, duck), not chicken only. Symbols that NCBI does not serve under
     one name are usually served under another species' naming.
  3. window sensitivity 1000 / 2000 / 5000 bp, reported separately instead of one
     fixed choice.

Stage is selected with --stage: targets | proteins | align | match.
"""
from __future__ import annotations

import argparse
import importlib
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
s53 = importlib.import_module("53_three_checks")
from oxphos_thread_lib import RES, SPECIES_ORDER, load_geneset, normalize_symbol  # noqa: E402

OUT = RES / "recovery_round2"
ASM = Path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies")
PROTDIR = Path("/path/to/storage/kleiber_birds_gemcat/ref/bird_prot")
NATIVE = ["Zebra", "Quail", "Chicken", "Duck"]
AUG = ["Peacock", "Turkey", "Emu"]
WINDOWS = (1000, 2000, 5000)
TAXA = {"9031": "chicken", "93934": "quail", "59729": "zebrafinch", "8839": "duck"}
HDR = re.compile(r"^>(\S+)\s+(\S+)\s+\[organism=")


def detection_matrix() -> pd.DataFrame:
    f = OUT / "detection_fib.tsv"
    if f.exists():
        return pd.read_csv(f, sep="\t", index_col=0)
    mats = {sp: s53.measure(sp, promoter_only=False, use_miniprot=True) for sp in SPECIES_ORDER}
    w = pd.DataFrame({sp: mats[sp]["fib"] for sp in SPECIES_ORDER})
    OUT.mkdir(parents=True, exist_ok=True)
    w.to_csv(f, sep="\t")
    return w


def stage_targets() -> dict[str, list[str]]:
    w = detection_matrix()
    det = w.notna() & (w > 0)
    genes = sorted(
        {normalize_symbol(g) for g in load_geneset("kegg_04714")}
        | {normalize_symbol(g) for g in load_geneset("kegg_00190")}
        - {""}
    )
    genes = [g for g in genes if g]
    d = det.reindex(genes).fillna(False)
    n_nat = d[NATIVE].sum(axis=1)
    out = {}
    for sp in AUG:
        t = sorted(d.index[(~d[sp]) & (n_nat >= 3)])
        out[sp] = t
        print(f"[{sp}] цели второго круга: {len(t)}")
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"species": k, "n": len(v), "genes": ";".join(v)} for k, v in out.items()]).to_csv(
        OUT / "targets_round2.tsv", sep="\t", index=False
    )
    allsym = sorted(set().union(*out.values()))
    (OUT / "symbols_round2.txt").write_text("\n".join(allsym) + "\n")
    print(f"объединение символов: {len(allsym)} -> {OUT / 'symbols_round2.txt'}")
    return out


def stage_proteins() -> None:
    PROTDIR.mkdir(parents=True, exist_ok=True)
    syms = OUT / "symbols_round2.txt"
    for taxid, label in TAXA.items():
        faa = PROTDIR / f"{label}.faa"
        if faa.exists() and faa.stat().st_size > 0:
            print(f"[{label}] уже есть, {sum(1 for l in faa.open() if l.startswith('>'))} белков")
            continue
        z = PROTDIR / f"{label}.zip"
        cmd = ["datasets", "download", "gene", "symbol", "--inputfile", str(syms),
               "--taxon", taxid, "--include", "protein", "--filename", str(z)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"[{label}] не скачалось: {r.stderr.strip()[:200]}")
            continue
        subprocess.run(["unzip", "-o", "-q", str(z), "-d", str(PROTDIR / label)], check=True)
        src = next((PROTDIR / label).rglob("protein.faa"), None)
        if src is None:
            print(f"[{label}] нет protein.faa")
            continue
        # tag the species into the header so best-hit picking can prefer chicken
        with faa.open("w") as fh:
            for line in src.read_text().splitlines():
                if line.startswith(">"):
                    p = line.split()
                    fh.write(f">{label}|{p[0][1:]} {p[1] if len(p) > 1 else 'NA'} [organism={label}]\n")
                else:
                    fh.write(line + "\n")
        subprocess.run(["rm", "-rf", str(z), str(PROTDIR / label)], check=False)
        print(f"[{label}] {sum(1 for l in faa.open() if l.startswith('>'))} белков")


def longest_per_symbol_species(faa: Path) -> dict[tuple[str, str], str]:
    best: dict[tuple[str, str], str] = {}
    pid = sym = None
    buf: list[str] = []

    def flush() -> None:
        if not sym or not buf:
            return
        seq = "".join(buf)
        k = (sym, pid)
        if k not in best or len(seq) > len(best[k]):
            best[k] = seq

    for line in faa.read_text().splitlines():
        if line.startswith(">"):
            flush()
            m = HDR.match(line)
            pid, sym, buf = (m.group(1), m.group(2), []) if m else (None, None, [])
        elif sym:
            buf.append(line.strip())
    flush()
    return best


def stage_align(threads: int) -> None:
    tg = pd.read_csv(OUT / "targets_round2.tsv", sep="\t")
    targets = {r["species"]: str(r["genes"]).split(";") for _, r in tg.iterrows()}
    prots: dict[str, list[tuple[str, str]]] = {}
    for label in TAXA.values():
        f = PROTDIR / f"{label}.faa"
        if not f.exists():
            continue
        for (sym, pid), seq in longest_per_symbol_species(f).items():
            prots.setdefault(sym, []).append((pid, seq))
    print(f"символов с хотя бы одним белком: {len(prots)}")
    for sp in AUG:
        want = [g for g in targets[sp] if g in prots]
        q = OUT / f"{sp}_q2.faa"
        with q.open("w") as fh:
            for g in want:
                for pid, seq in prots[g]:
                    fh.write(f">{pid}##{g}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i : i + 60] + "\n")
        miss = [g for g in targets[sp] if g not in prots]
        print(f"[{sp}] запросов на {len(want)} генов; без белка остались {len(miss)}: {', '.join(miss[:20])}")
        gff = OUT / f"{sp}_q2.gff"
        if gff.exists() and gff.stat().st_size:
            print(f"[{sp}] {gff.name} уже есть")
            continue
        genome = ASM / sp / ("genome_localnames.fa" if sp == "Peacock" else "genome.fa")
        with gff.open("w") as fh:
            subprocess.run(["miniprot", "-t", str(threads), "-Iu", "--gff", str(genome), str(q)],
                           stdout=fh, check=True)


def stage_match() -> None:
    EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
    rows = []
    for sp in AUG:
        gff = OUT / f"{sp}_q2.gff"
        if not gff.exists():
            continue
        best: dict[str, tuple[float, dict]] = {}
        for line in gff.read_text().splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            at = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
            tgt = at.get("Target", "").split()[0]
            if "##" not in tgt:
                continue
            sym = tgt.split("##")[1]
            seq = f[0].split()[0]
            if seq.startswith("ENA|"):
                seq = seq.split("|")[-1]
            sc = float(f[5]) if f[5] != "." else 0.0
            st = f[6]
            cand = dict(gene=sym, seqname=seq, strand=st,
                        five=float(f[3] if st == "+" else f[4]),
                        identity=float(at.get("Identity", 0) or 0))
            if sym not in best or sc > best[sym][0]:
                best[sym] = (sc, cand)
        d = pd.read_csv(EXP / f"{sp}_fib_epi_exp.txt", sep="\t", low_memory=False)
        fib = [c for c in d.columns if c.startswith("norm.C")]
        d["mid"] = (d["start"] + d["end"]) / 2.0
        d["fibv"] = d[fib].apply(pd.to_numeric, errors="coerce").fillna(0.0).where(
            lambda x: x >= 1.0, 0.0).sum(axis=1) / max(len(fib), 1)
        by = {k: v for k, v in d.groupby("seqnames", observed=True)}
        for sym, (_, c) in best.items():
            sub = by.get(c["seqname"])
            row = dict(species=sp, gene=sym, seqname=c["seqname"], identity=round(c["identity"], 3))
            for w in WINDOWS:
                hit = 0
                if sub is not None:
                    near = sub[(sub["strand"] == c["strand"]) & ((sub["mid"] - c["five"]).abs() <= w)]
                    hit = int((near["fibv"] > 0).sum())
                row[f"peaks_{w}"] = hit
            rows.append(row)
        print(f"[{sp}] локусов найдено {len(best)}")
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "round2_hits.tsv", sep="\t", index=False)
    print("\nвосстановлено дополнительно, по окну:")
    for sp, s in t.groupby("species"):
        line = "  ".join(f"±{w}: {int((s[f'peaks_{w}'] > 0).sum())}/{len(s)}" for w in WINDOWS)
        print(f"  {sp:8s} {line}")
    print(f"\nwrote {OUT / 'round2_hits.tsv'}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["targets", "proteins", "align", "match"])
    ap.add_argument("--threads", type=int, default=64)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {"targets": stage_targets, "proteins": stage_proteins,
     "align": lambda: stage_align(a.threads), "match": stage_match}[a.stage]()


if __name__ == "__main__":
    main()
