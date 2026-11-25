#!/usr/bin/env python3
"""The 19 genes that carried the old ±1000 decline, measured via CAGE.

The old window called them detected in all seven because a peak sat near a
miniprot 5'. For Augustus species that 5' was often a neighbour (RGS22, ADCK2,
GATC, CABP7). This script does not change the set: it asks, for each of those
19, where the gene actually is and what CAGE belongs to it.

Two independent placements, then CAGE at that place:

  1. CDS from the CHICKEN protein only (not the highest-scoring turkey/emu
     isoform — those can extend 5' onto a neighbour, which is how emu COX6C
     ate RGS22). Promoter = CAGE on the transcribed strand in the 5' catchment
     up to the next gene, window [-5000, +300] bp of that CDS 5'. A peak already
     named as a different gene is kept only if it overlaps the lifted chicken
     promoter (the name was wrong, the place is right).

  2. The chicken CAGE peak that already carries this gene's name, sequence
     lifted with minimap2 onto each genome. CAGE overlapping the lift (±200 bp)
     is the homologous promoter.

A species counts as having the gene if the chicken-protein CDS is there AND
there is CAGE from (1) or (2). Native-annotation species almost always have
both; Augustus are the test.
"""
from __future__ import annotations

import importlib
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
s62 = importlib.import_module("62_build_library_matrices")
from oxphos_thread_lib import (  # noqa: E402
    RES,
    RUSLAN,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    normalize_symbol,
    peak_to_symbol_de,
)

ASM = Path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies")
OUT = RES / "nineteen_via_cage"
NINETEEN = [
    "ATP5F1E", "ATP6V0E1", "ATP6V0E2", "COX6A1", "COX6C", "COX7A2",
    "NDUFA5", "NDUFA7", "NDUFAB1", "NDUFB1", "NDUFB2", "NDUFB5", "NDUFB6",
    "NDUFB8", "NDUFC2", "NDUFS1", "SDHD", "UQCR10", "UQCRC1",
]
WANT = set(NINETEEN)
LOGM = np.log10([SPECIES[s]["mass_g"] for s in SPECIES_ORDER])
THR = 1.0
UP, DOWN = 5000, 300
LIFT_PAD = 200
GENOME = {
    "Chicken": ASM / "Chicken" / "genome.fa",
    "Duck": ASM / "Duck" / "genome.fa",
    "Quail": ASM / "Quail" / "genome.fa",
    "Zebra": ASM / "Zebra" / "genome.fa",
    "Peacock": ASM / "Peacock" / "genome_localnames.fa",
    "Turkey": ASM / "Turkey" / "genome.fa",
    "Emu": ASM / "Emu" / "genome.fa",
}
PRODUCT = {
    "ATP5F1E": "ATP synthase F1 ε, ~50 aa",
    "ATP6V0E1": "V-ATPase c″ subunit",
    "ATP6V0E2": "V-ATPase c″ isoform 2",
    "COX6A1": "Complex IV subunit 6A liver isoform",
    "COX6C": "Complex IV subunit 6C",
    "COX7A2": "Complex IV subunit 7A2",
    "NDUFA5": "Complex I supernumerary A5",
    "NDUFA7": "Complex I supernumerary A7",
    "NDUFAB1": "Complex I ACP / SDAP",
    "NDUFB1": "Complex I supernumerary B1 (MNLL)",
    "NDUFB2": "Complex I supernumerary B2 (AGGG)",
    "NDUFB5": "Complex I supernumerary B5",
    "NDUFB6": "Complex I supernumerary B6",
    "NDUFB8": "Complex I supernumerary B8",
    "NDUFC2": "Complex I supernumerary C2 (B14.5b)",
    "NDUFS1": "Complex I core Fe-S 75 kDa",
    "SDHD": "Complex II cytochrome b small subunit",
    "UQCR10": "Complex III subunit 10 (QCR9)",
    "UQCRC1": "Complex III core subunit 1",
}


def contig(s: str) -> str:
    s = str(s).split()[0]
    return s.split("|")[-1]


def complex_of(g: str) -> str:
    if g.startswith("NDUF"):
        return "CI"
    if g.startswith("SDH"):
        return "CII"
    if g.startswith("UQCR"):
        return "CIII"
    if g.startswith("COX"):
        return "CIV"
    if g.startswith("ATP5"):
        return "CV"
    if g.startswith("ATP6V"):
        return "V-ATPase"
    return "other"


def sh(cmd: list[str], **kw):
    return subprocess.run(cmd, check=True, **kw)


def chicken_genes(gff: Path) -> pd.DataFrame:
    rows = []
    with gff.open() as fh:
        for line in fh:
            if line[0] == "#" or "\tgene\t" not in line:
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            if "gene_biotype=protein_coding" not in f[8] and "biotype=protein_coding" not in f[8]:
                # NCBI uses gene_biotype=
                if "gene_biotype=" in f[8] and "protein_coding" not in f[8]:
                    continue
            m = re.search(r"(?:gene=|Name=)([^;]+)", f[8])
            if not m:
                continue
            sym = normalize_symbol(m.group(1))
            if not sym:
                continue
            st = f[6]
            five = int(f[3] if st == "+" else f[4])
            rows.append(dict(symbol=sym, seqname=contig(f[0]), strand=st,
                             start=int(f[3]), end=int(f[4]), five=five,
                             product=re.search(r"description=([^;]+)", f[8])))
    g = pd.DataFrame(rows).drop_duplicates("symbol", keep="first")
    g["product"] = g["product"].map(lambda m: m.group(1) if m else "")
    return g


def flanks_of(genes: pd.DataFrame, sym: str) -> tuple[str, str]:
    row = genes[genes["symbol"] == sym]
    if row.empty:
        return "", ""
    r = row.iloc[0]
    same = genes[genes["seqname"] == r["seqname"]].sort_values("start").reset_index(drop=True)
    hits = same.index[same["symbol"] == sym].tolist()
    if not hits:
        return "", ""
    i = int(hits[0])
    left = same.iloc[i - 1]["symbol"] if i > 0 else ""
    right = same.iloc[i + 1]["symbol"] if i < len(same) - 1 else ""
    return str(left), str(right)


def chicken_protein_loci(gff: Path, symbols: set[str]) -> dict[str, dict]:
    """One locus per symbol: the chicken_* query with the highest identity."""
    best: dict[str, tuple[float, dict]] = {}
    if not gff.exists():
        return {}
    with gff.open() as fh:
        for line in fh:
            if line[0] == "#" or "\tmRNA\t" not in line:
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            m = re.search(r"Target=(\S+)", f[8])
            if not m or "##" not in m.group(1):
                continue
            src, raw = m.group(1).split("##", 1)
            if not src.lower().startswith("chicken"):
                continue
            sym = normalize_symbol(raw)
            if sym not in symbols:
                continue
            ident = re.search(r"Identity=([0-9.]+)", f[8])
            idv = float(ident.group(1)) if ident else 0.0
            sc = float(f[5]) if f[5] != "." else 0.0
            st = f[6]
            seq = contig(f[0])
            rec = dict(seqname=seq, strand=st,
                       start=int(f[3]), end=int(f[4]),
                       five=float(f[3] if st == "+" else f[4]),
                       ident=idv, score=sc, query=src)
            key = (idv, sc)
            prev = best.get(sym)
            if prev is None or key > (prev[0], prev[1]["score"]):
                best[sym] = (idv, rec)
    return {k: v[1] for k, v in best.items()}


def all_fives(gff: Path) -> dict[str, list[tuple[str, float, str]]]:
    """seqname -> list of (symbol, five, strand) for catchment barriers."""
    best: dict[str, tuple[float, tuple]] = {}
    with gff.open() as fh:
        for line in fh:
            if line[0] == "#" or "\tmRNA\t" not in line:
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            m = re.search(r"Target=(\S+)", f[8])
            if not m or "##" not in m.group(1):
                continue
            src = m.group(1)
            if not src.lower().startswith("chicken"):
                continue
            sym = normalize_symbol(src.split("##", 1)[1])
            if not sym:
                continue
            sc = float(f[5]) if f[5] != "." else 0.0
            ident = re.search(r"Identity=([0-9.]+)", f[8])
            idv = float(ident.group(1)) if ident else 0.0
            st = f[6]
            seq = contig(f[0])
            five = float(f[3] if st == "+" else f[4])
            rec = (idv, sc, (sym, five, st, seq))
            if sym not in best or (idv, sc) > (best[sym][0], best[sym][1]):
                best[sym] = rec
    by: dict[str, list[tuple[str, float, str]]] = {}
    for _, _, (sym, five, st, seq) in best.values():
        by.setdefault(seq, []).append((sym, five, st))
    for seq in by:
        by[seq].sort(key=lambda t: t[1])
    return by


def oriented(mid: float, five: float, strand: str) -> float:
    d = mid - five
    return -d if strand == "-" else d


def catchment_lohi(five: float, strand: str, barriers: list[tuple[str, float, str]],
                   self_sym: str) -> tuple[float, float]:
    """Window [-UP, +DOWN] of CDS 5', clipped by the next gene 5' upstream/down."""
    if strand == "+":
        up_b = [b[1] for b in barriers if b[0] != self_sym and b[1] < five]
        down_b = [b[1] for b in barriers if b[0] != self_sym and b[1] > five]
        lo = max(five - UP, max(up_b) if up_b else five - UP)
        hi = min(five + DOWN, min(down_b) if down_b else five + DOWN)
    else:
        # minus: upstream = higher coords
        up_b = [b[1] for b in barriers if b[0] != self_sym and b[1] > five]
        down_b = [b[1] for b in barriers if b[0] != self_sym and b[1] < five]
        lo = max(five - DOWN, max(down_b) if down_b else five - DOWN)
        hi = min(five + UP, min(up_b) if up_b else five + UP)
    return lo, hi


def load_peaks(species: str) -> pd.DataFrame:
    d = pd.read_csv(s62.EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
    cols = lib_cols(d, "fib")
    num = d[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    num = num.where(num >= THR, 0.0)
    d = d.assign(mid=(d["start"] + d["end"]) / 2.0,
                 fib=num.mean(axis=1),
                 seq=d["seqnames"].map(contig))
    return d[d["fib"] > 0].copy()


def ensure_faidx(fa: Path) -> None:
    if not (Path(str(fa) + ".fai")).exists():
        sh(["samtools", "faidx", str(fa)])


def chicken_promoter_fasta(peaks: pd.DataFrame, fa: Path, out_fa: Path) -> pd.DataFrame:
    """One strongest named CAGE peak per gene on chicken."""
    ensure_faidx(fa)
    rows = []
    bed = OUT / "chicken_promoters.bed"
    with bed.open("w") as fh:
        for g in NINETEEN:
            hit = peaks[peaks["symbol"] == g].sort_values("fib", ascending=False)
            if hit.empty:
                continue
            r = hit.iloc[0]
            start = max(0, int(r["start"]) - 50)
            end = int(r["end"]) + 50
            name = f"{g}|{r['id']}"
            fh.write(f"{r['seqnames']}\t{start}\t{end}\t{name}\t0\t{r['strand']}\n")
            rows.append(dict(gene=g, peak_id=r["id"], seqname=r["seqnames"],
                             start=int(r["start"]), end=int(r["end"]),
                             strand=r["strand"], fib=float(r["fib"])))
    sh(["bedtools", "getfasta", "-fi", str(fa), "-bed", str(bed),
        "-s", "-nameOnly", "-fo", str(out_fa)])
    return pd.DataFrame(rows)


def lift_one(genome: Path, query: Path, paf: Path) -> pd.DataFrame:
    ensure_faidx(genome)
    with paf.open("w") as fh:
        sh(["minimap2", "-t", "8", "-x", "asm20", "-k", "11", "-w", "5",
            "--secondary=no", "-N", "3", str(genome), str(query)],
           stdout=fh, stderr=subprocess.DEVNULL)
    cols = ["q", "qlen", "qs", "qe", "strand", "t", "tlen", "ts", "te",
            "nmatch", "alen", "mapq"]
    if paf.stat().st_size == 0:
        return pd.DataFrame(columns=cols + ["gene"])
    raw = pd.read_csv(paf, sep="\t", header=None, usecols=range(12), names=cols)
    raw["gene"] = raw["q"].str.split("|").str[0]
    raw["t"] = raw["t"].map(contig)
    raw["ident"] = raw["nmatch"] / raw["alen"].clip(lower=1)
    raw = raw.sort_values(["gene", "mapq", "ident"], ascending=[True, False, False])
    return raw.drop_duplicates("gene", keep="first")


def score_peaks(near: pd.DataFrame) -> float:
    if near is None or near.empty:
        return float("nan")
    v = float(near["fib"].sum())
    return v if v > 0 else float("nan")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("курица: гены NCBI + CAGE промоторы", flush=True)
    ch_genes = chicken_genes(ASM / "Chicken" / "genomic.gff")
    peaks = {sp: load_peaks(sp) for sp in SPECIES_ORDER}
    ch_prom = chicken_promoter_fasta(peaks["Chicken"], GENOME["Chicken"],
                                     OUT / "chicken_promoters.fa")
    ch_prom.to_csv(OUT / "chicken_promoters.tsv", sep="\t", index=False)

    cards = []
    for g in NINETEEN:
        left, right = flanks_of(ch_genes, g)
        prow = ch_prom[ch_prom["gene"] == g]
        grow = ch_genes[ch_genes["symbol"] == g]
        cards.append(dict(
            gene=g, complex=complex_of(g), product=PRODUCT.get(g, ""),
            chicken_chr="" if grow.empty else grow.iloc[0]["seqname"],
            chicken_tss="" if grow.empty else int(grow.iloc[0]["five"]),
            chicken_strand="" if grow.empty else grow.iloc[0]["strand"],
            flank_left=left, flank_right=right,
            cage_tpm=None if prow.empty else round(float(prow.iloc[0]["fib"]), 1),
            cage_peak="" if prow.empty else prow.iloc[0]["peak_id"],
        ))
    pd.DataFrame(cards).to_csv(OUT / "identity_cards.tsv", sep="\t", index=False)
    print("карточки:", OUT / "identity_cards.tsv", flush=True)

    lifts = {}
    for sp in SPECIES_ORDER:
        if sp == "Chicken":
            continue
        print(f"lift → {sp}", flush=True)
        lifts[sp] = lift_one(GENOME[sp], OUT / "chicken_promoters.fa",
                             OUT / f"lift_{sp}.paf")
        print(f"  mapped {lifts[sp]['gene'].nunique() if len(lifts[sp]) else 0} / 19",
              flush=True)

    loci = {sp: chicken_protein_loci(ASM / sp / "reanno_named.gff", WANT)
            for sp in SPECIES_ORDER}
    barriers = {sp: all_fives(ASM / sp / "reanno_named.gff") for sp in SPECIES_ORDER}

    # chicken promoters as trivial lifts
    lifts["Chicken"] = pd.DataFrame([
        dict(gene=r.gene, t=r.seqname, ts=r.start - 50, te=r.end + 50,
             strand=r.strand, mapq=60, ident=1.0)
        for r in ch_prom.itertuples()
    ])

    rows = []
    for g in NINETEEN:
        for sp in SPECIES_ORDER:
            rec = dict(gene=g, species=sp, complex=complex_of(g))
            loc = loci[sp].get(g)
            rec["cds"] = False
            if loc:
                rec.update(cds=True, cds_ident=round(loc["ident"], 3),
                           cds_seq=loc["seqname"], cds_strand=loc["strand"],
                           cds_five=int(loc["five"]))
            pk = peaks[sp]
            catch_val = np.nan
            catch_n = 0
            catch_names = ""
            if loc:
                barr = barriers[sp].get(loc["seqname"], [])
                lo, hi = catchment_lohi(loc["five"], loc["strand"], barr, g)
                rec["catch_lo"] = int(lo)
                rec["catch_hi"] = int(hi)
                sub = pk[(pk["seq"] == loc["seqname"]) & (pk["strand"] == loc["strand"])
                         & (pk["mid"] >= lo) & (pk["mid"] <= hi)].copy()
                lf = lifts.get(sp, pd.DataFrame())
                lift_hit = lf[lf["gene"] == g] if len(lf) else lf
                lift_iv = None
                if lift_hit is not None and len(lift_hit):
                    lr = lift_hit.iloc[0]
                    lift_iv = (contig(lr["t"]), int(lr["ts"]) - LIFT_PAD, int(lr["te"]) + LIFT_PAD)
                keep_idx = []
                for i, r in sub.iterrows():
                    nom = r["symbol"] if isinstance(r["symbol"], str) else ""
                    if nom == g or nom == "":
                        keep_idx.append(i)
                        continue
                    if lift_iv and r["seq"] == lift_iv[0] and r["end"] >= lift_iv[1] and r["start"] <= lift_iv[2]:
                        keep_idx.append(i)
                        continue
                    # foreign name, not on the lifted promoter: refuse
                kept = sub.loc[keep_idx] if keep_idx else sub.iloc[0:0]
                catch_val = score_peaks(kept)
                catch_n = int(len(kept))
                catch_names = ",".join(sorted({str(x) for x in kept["symbol"] if isinstance(x, str) and x})[:5])
            rec["cage_catchment"] = catch_val
            rec["catch_n"] = catch_n
            rec["catch_names"] = catch_names

            lf = lifts.get(sp, pd.DataFrame())
            rec["lifted"] = False
            lift_val = np.nan
            if lf is not None and len(lf) and g in set(lf["gene"]):
                lr = lf[lf["gene"] == g].iloc[0]
                rec.update(lifted=True, lift_seq=contig(lr["t"]),
                           lift_start=int(lr["ts"]), lift_end=int(lr["te"]),
                           lift_strand=lr["strand"], lift_mapq=int(lr["mapq"]),
                           lift_ident=round(float(lr["ident"]), 3))
                sub = pk[(pk["seq"] == contig(lr["t"])) &
                         (pk["end"] >= int(lr["ts"]) - LIFT_PAD) &
                         (pk["start"] <= int(lr["te"]) + LIFT_PAD)]
                # same strand as lift if we can read it as +/- 
                st = str(lr["strand"])
                if st in "+-":
                    sub = sub[sub["strand"] == st]
                lift_val = score_peaks(sub)
                rec["lift_peak_names"] = ",".join(sorted({str(x) for x in sub["symbol"] if isinstance(x, str) and x})[:5])
            rec["cage_lift"] = lift_val

            # chosen measurement: lift if it landed next to this gene's CDS, else catchment
            chosen = np.nan
            how = "нет"
            if loc and rec.get("lifted") and rec.get("lift_seq") == loc["seqname"]:
                dist = abs(((int(rec["lift_start"]) + int(rec["lift_end"])) / 2) - loc["five"])
                rec["lift_vs_cds_bp"] = int(dist)
                if dist <= 20000 and np.isfinite(lift_val):
                    chosen, how = lift_val, "lift+CDS"
            if not np.isfinite(chosen) and np.isfinite(catch_val):
                chosen, how = catch_val, "catchment"
            if not np.isfinite(chosen) and loc and rec.get("lifted") and np.isfinite(lift_val):
                chosen, how = lift_val, "lift_only"
            rec["tpm"] = chosen
            rec["how"] = how if np.isfinite(chosen) else ("CDS_без_CAGE" if loc else "нет_локуса")
            rows.append(rec)

    t = pd.DataFrame(rows)
    t.to_csv(OUT / "per_gene_species.tsv", sep="\t", index=False)

    M = t.pivot(index="gene", columns="species", values="tpm").reindex(index=NINETEEN, columns=SPECIES_ORDER)
    M.to_csv(OUT / "M_nineteen_tpm.tsv", sep="\t")
    det = M.notna() & (M > 0)
    print("\nобнаружено (CAGE при курином CDS или lift):")
    print(" ", " ".join(f"{sp} {int(det[sp].sum())}" for sp in SPECIES_ORDER))
    n7 = [g for g in NINETEEN if det.loc[g].all()]
    print(f" во всех семи: {len(n7)} — {', '.join(n7) if n7 else '—'}")

    def fit_of(genes: list[str], label: str) -> None:
        if not genes:
            print(f"\n{label}: пусто")
            return
        y = np.log10(M.reindex(genes)[SPECIES_ORDER].mean(axis=0).to_numpy(float))
        r = stats.linregress(LOGM, y)
        print(f"\n{label}: n={len(genes)} r={r.rvalue:+.3f} p={r.pvalue:.4f}")
        print(" ", " ".join(f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y)))

    fit_of(n7, "19, только гены с CAGE во всех семи")
    # holes as missing skipped in mean — forbidden, but report n per species mean of found
    y_found = []
    ns = []
    for sp in SPECIES_ORDER:
        v = M[sp].to_numpy(float)
        ok = np.isfinite(v) & (v > 0)
        ns.append(int(ok.sum()))
        y_found.append(np.log10(v[ok].mean()) if ok.any() else np.nan)
    r = stats.linregress(LOGM, np.array(y_found, float))
    print(f"\nсреднее только найденных (дыры завышают): n={ns} r={r.rvalue:+.3f} p={r.pvalue:.4f}")
    print(" ", " ".join(f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y_found)))

    # compare to old window library means for the same 19
    old = pd.read_csv(RES / "wgcna_input" / "M_fib_w1000.tsv", sep="\t", index_col=0)
    spcol = pd.Series([c.split("|")[0] for c in old.columns], index=old.columns)
    per = old.T.groupby(spcol.to_numpy()).mean().T.reindex(columns=SPECIES_ORDER)
    y = np.log10(per.reindex(NINETEEN).mean(axis=0).to_numpy(float))
    r = stats.linregress(LOGM, y)
    print(f"\nте же 19 на СТАРОМ окне ±1000: r={r.rvalue:+.3f} p={r.pvalue:.4f}")
    print(" ", " ".join(f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y)))

    print("\nкак измерен каждый ген × вид (how):")
    print(t.groupby(["species", "how"]).size().unstack(fill_value=0).reindex(SPECIES_ORDER).to_string())
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
