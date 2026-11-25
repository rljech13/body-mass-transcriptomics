#!/usr/bin/env python3
"""Gene x library matrices rebuilt on a symmetric footing.

Two things are changed against script 62, and both were forced by measurement.

Window. Script 62 took peaks within 1000 bp of a gene's 5' end in either
direction. Among genes detected in exactly four of the seven species, the three
missing are Peacock, Turkey and Emu in 1036, 1016 and 1030 cases out of 1073 —
the loss runs along the axis under test. The cause is that their annotations mark
where translation starts while CAGE marks where transcription starts: for genes
whose nearest peak falls outside the old window the median offset is -1766 bp in
Peacock, -2066 in Turkey and -2126 in Emu, against +1399 in Chicken with only 41%
upstream. So the window is opened upstream — for all seven species alike, which
is what check B of script 53 requires. Chicken has 627 genes outside the old
window and gains little; Peacock has 2805 and gains most.

A window reaching 5 kb upstream can cover a neighbouring gene's promoter, so a
peak is no longer counted for every locus that covers it. Two rules, in order:

  1. If the peak table already names a gene that we have a locus for, the peak
     stays with that gene, on any contig. Distance to the miniprot CDS start is
     not consulted: RGS22's promoter sits 75 kb from its own CDS, and a 20 kb
     cap would refuse the name and also refuse to give the peak to COX6C, which
     is how emu COX6C went to zero. A peak named as a different gene is never
     given to a neighbour (COX6C/RGS22, SDHD/C24H11ORF57, NDUFB1/TRIP11).
  2. Unnamed peaks go to the nearest locus inside the window, provided no other
     gene's 5' lies between them.

Loci. Script 58 projected a proteome onto three genomes out of seven and gave
each a different query, which is why locus counts ran 10 953 in Chicken against
18 618 in Emu. Script 70 replaces that with one combined query — one protein per
gene from chicken, duck, quail, zebra, turkey and emu, 79 570 proteins over
16 779 named symbols — projected onto all seven genomes. Symbols matter: the
first pass kept LOC-prefixed placeholders, which do not agree between species and
planted some 21 600 phantom loci per genome; with peaks going to the nearest
locus those phantoms stole real promoters, and ATP5MF and UQCRC2 fell from
detected in all seven species to detected in one. With LOC dropped the layer
comes out even: 16 054 to 16 523 symbols per species, under 3% spread.

Two ways of using it are offered and the choice is made on evidence, not taste:

  annot_first  the species' own annotation first, reanno_common filling the rest.
               More accurate per species (real TSS where the annotation is good),
               but the accuracy itself is unequal across species.
  common_only  reanno_common for everyone. Less accurate — every locus is a
               coding start — but identical in kind for all seven.

Universe is unchanged: a gene enters if non-zero in at least one library of that
tissue in EVERY species. That requirement is not an artefact; it is the only
threshold at which detection rate does not track body mass (r = 0.000 at k = 7
against -0.670 at k >= 6 and -0.790 at k >= 4 under the old window).
"""
from __future__ import annotations

import argparse
import importlib
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
s57 = importlib.import_module("57_working_sets_round2")
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
OUT = RES / "wgcna_input_v2"
THR = 1.0


def loci_from_reanno(gff: Path) -> dict[str, dict]:
    """Best-scoring model per NORMALISED symbol.

    Normalising is not cosmetic here. The six source annotations disagree on the
    ATP synthase nomenclature, so the query carries ATP5C1 next to ATP5F1C,
    ATP5H next to ATP5PD and so on. Left as they are, one gene gets two loci at
    the same place which then compete for the same peak, and the gene drops out
    of most species: ATP5F1E and ATP5PF survived in one species each before this
    was applied. PREFERRED_ALIAS is the table the rest of the project already
    uses (oxphos_thread_lib.normalize_symbol).
    """
    best: dict[str, tuple[float, dict]] = {}
    if not gff.exists():
        return {}
    with gff.open() as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            m = re.search(r"Target=(\S+)", f[8])
            if not m or "##" not in m.group(1):
                continue
            sym = normalize_symbol(m.group(1).split("##")[1])
            if not sym:
                continue
            seq = f[0].split()[0]
            if seq.startswith("ENA|"):
                seq = seq.split("|")[-1]
            sc = float(f[5]) if f[5] != "." else 0.0
            st = f[6]
            if sym not in best or sc > best[sym][0]:
                best[sym] = (sc, dict(seqname=seq, strand=st,
                                      five=float(f[3] if st == "+" else f[4])))
    return {k: v[1] for k, v in best.items()}


def loci_of(species: str, d: pd.DataFrame, mode: str) -> tuple[dict[str, dict], dict[str, int]]:
    """symbol -> 5' locus, plus how many each layer contributed."""
    fp: dict[str, dict] = {}
    src: dict[str, int] = {}

    if mode == "annot_first":
        ann = d.dropna(subset=["geneStart", "geneEnd"])
        ann = ann[ann["symbol"].ne("")].copy()
        ann["sg"] = ann["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(ann["strand"])
        for sym, g in ann.groupby("symbol", observed=True):
            r = g.iloc[0]
            st = r["sg"]
            fp[sym] = dict(seqname=r["seqnames"], strand=st,
                           five=float(r["geneStart"] if st == "+" else r["geneEnd"]))
        src["аннотация"] = len(fp)

    common = loci_from_reanno(ASM / species / "reanno_named.gff")
    added = 0
    for sym, c in common.items():
        if sym not in fp:
            fp[sym] = c
            added += 1
    src["reanno_named"] = added

    if mode == "common_then_annot":
        ann = d.dropna(subset=["geneStart", "geneEnd"])
        ann = ann[ann["symbol"].ne("")].copy()
        ann["sg"] = ann["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"}).fillna(ann["strand"])
        back = 0
        for sym, g in ann.groupby("symbol", observed=True):
            if sym in fp:
                continue
            r = g.iloc[0]
            st = r["sg"]
            fp[sym] = dict(seqname=r["seqnames"], strand=st,
                           five=float(r["geneStart"] if st == "+" else r["geneEnd"]))
            back += 1
        src["аннотация добором"] = back
    return fp, src


def _oriented_offset(mid: float, five: float, strand: str) -> float:
    d = mid - five
    return -d if strand == "-" else d


def _table_fives(peaks: pd.DataFrame) -> dict[tuple[str, str], np.ndarray]:
    """Every named gene's 5' from the peak table, used only as a jump barrier."""
    if peaks.empty or "geneStart" not in peaks.columns:
        return {}
    sg = peaks["geneStrand"].map({1: "+", 2: "-", "+": "+", "-": "-"})
    if "strand" in peaks.columns:
        sg = sg.fillna(peaks["strand"])
    rows: dict[tuple[str, str], list[float]] = {}
    gs = pd.to_numeric(peaks["geneStart"], errors="coerce").to_numpy(float)
    ge = (pd.to_numeric(peaks["geneEnd"], errors="coerce").to_numpy(float)
          if "geneEnd" in peaks.columns else np.full(len(peaks), np.nan))
    syms = peaks["symbol"].to_numpy(object) if "symbol" in peaks.columns else np.array([""] * len(peaks), dtype=object)
    for seq, st, gsv, gev, sym in zip(
        peaks["seqnames"].to_numpy(object), sg.to_numpy(object), gs, ge, syms,
    ):
        if not isinstance(sym, str) or not sym or st not in ("+", "-"):
            continue
        if not np.isfinite(gsv):
            continue
        five = float(gsv if st == "+" else (gev if np.isfinite(gev) else gsv))
        rows.setdefault((str(seq), str(st)), []).append(five)
    return {k: np.unique(np.asarray(v, float)) for k, v in rows.items()}


def _has_between(mid: float, five: float, others: np.ndarray) -> bool:
    if others.size == 0:
        return False
    lo, hi = (mid, five) if mid <= five else (five, mid)
    i0 = int(np.searchsorted(others, lo, side="right"))
    i1 = int(np.searchsorted(others, hi, side="left"))
    return i1 > i0


def assign_peaks(fp: dict[str, dict], peaks: pd.DataFrame, up: int, down: int
                 ) -> tuple[pd.Series, dict[str, int]]:
    """Named peak stays with its gene, regardless of the miniprot CDS start.

    A peak named as some other gene is never given to a neighbour. Remaining
    unnamed peaks go to the nearest locus in [-up, +down] that does not skip
    over another gene's 5'.
    """
    owner = pd.Series(index=peaks.index, dtype=object)
    locked = np.zeros(len(peaks), dtype=bool)
    n_named = n_locked = 0
    named_i: list = []
    named_s: list[str] = []
    if "symbol" in peaks.columns and not peaks.empty:
        for k, (i, seq, strand, mid, raw) in enumerate(zip(
            peaks.index.to_numpy(),
            peaks["seqnames"].to_numpy(object),
            peaks["strand"].to_numpy(object),
            peaks["mid"].to_numpy(float),
            peaks["symbol"].to_numpy(object),
        )):
            if not isinstance(raw, str) or not raw:
                continue
            sym = normalize_symbol(raw)
            if not sym:
                continue
            c = fp.get(sym)
            if c is not None:
                named_i.append(i)
                named_s.append(sym)
                locked[k] = True
                n_named += 1
            else:
                locked[k] = True
                n_locked += 1
    if named_i:
        owner.loc[named_i] = named_s

    free = peaks.loc[~locked]
    n_near = 0
    if not free.empty:
        loci = pd.DataFrame(
            [dict(symbol=s, seqname=c["seqname"], strand=c["strand"], five=c["five"])
             for s, c in fp.items()]
        )
        extra = _table_fives(peaks)
        by_peak = {k: v for k, v in free.groupby(["seqnames", "strand"], observed=True)}
        lab: list = []
        won: list[str] = []
        for (seq, strand), lg in loci.groupby(["seqname", "strand"], observed=True):
            pg = by_peak.get((seq, strand))
            if pg is None or pg.empty:
                continue
            lg = lg.sort_values("five")
            five = lg["five"].to_numpy(float)
            syms = lg["symbol"].to_numpy(object)
            mid = pg["mid"].to_numpy(float)
            idx = pg.index.to_numpy()
            block = extra.get((str(seq), str(strand)), np.array([], float))
            if five.size:
                block = np.unique(np.concatenate([block, five])) if block.size else five
            lo = mid - (down if strand == "+" else up)
            hi = mid + (up if strand == "+" else down)
            i0 = np.searchsorted(five, lo, side="left")
            i1 = np.searchsorted(five, hi, side="right")
            for k in np.flatnonzero(i1 > i0):
                a, b = int(i0[k]), int(i1[k])
                dd = mid[k] - five[a:b]
                if strand == "-":
                    dd = -dd
                order = np.argsort(np.abs(dd))
                chosen = None
                for j in order:
                    cand_five = five[a:b][j]
                    others = block[block != cand_five]
                    if _has_between(float(mid[k]), float(cand_five), others):
                        continue
                    chosen = syms[a + int(j)]
                    break
                if chosen is None:
                    continue
                lab.append(idx[k])
                won.append(chosen)
        if lab:
            owner.loc[lab] = won
            n_near = len(lab)
    stats_n = dict(якорь=n_named, чужое_имя=n_locked, nearest=n_near)
    return owner, stats_n


def fib_species_means(M: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    sp = pd.Series([c.split("|")[0] for c in M.columns], index=M.columns)
    per = M.T.groupby(sp.to_numpy()).mean().T
    hit = [g for g in genes if g in per.index]
    y = np.log10(per.reindex(hit)[SPECIES_ORDER].mean(axis=0).to_numpy(float))
    return y, hit


def report_oxphos(M: pd.DataFrame, core: list[str], uf: list[str],
                  logm: np.ndarray, tag: str) -> None:
    """Compare the working set on this matrix vs the previous nearest-only run."""
    watch = ["COX6C", "SDHD", "NDUFB1"]
    old_tag = tag.replace("_anchored", "")
    old_path = OUT / f"M_fib_{old_tag}.tsv"
    w1000 = RES / "wgcna_input" / "M_fib_w1000.tsv"

    def one(label: str, mat: pd.DataFrame, genes: list[str]) -> None:
        y, hit = fib_species_means(mat, genes)
        r = stats.linregress(logm, y)
        means = " ".join(f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y))
        print(f"  {label}: n={len(hit)} r={r.rvalue:+.3f} p={r.pvalue:.4f}")
        print(f"    {means}")

    print("\nOXPHOS, log10(среднее по генам), фибробласты:")
    one(f"якорь {tag}, 55∩вселенная", M, [g for g in core if g in uf])
    one(f"якорь {tag}, все 55 (дыры как есть)", M, core)
    if old_path.exists():
        old = pd.read_csv(old_path, sep="\t", index_col=0)
        u_old = set(old.index)
        g47 = [g for g in core if g in u_old]
        print(f"\n  те же {len(g47)} генов, что были на nearest-only:")
        one("nearest-only", old, g47)
        one("якорь", M, g47)
        print("\n  три гена, которые украли соседний пик:")
        for g in watch:
            print(f"    {g}")
            for lab, mat in (("nearest", old), ("якорь", M)):
                if g not in mat.index:
                    print(f"      {lab}: нет в матрице")
                    continue
                sp = pd.Series([c.split("|")[0] for c in mat.columns])
                y = mat.loc[g].groupby(sp.to_numpy()).mean().reindex(SPECIES_ORDER)
                print("      " + lab + ": " + " ".join(
                    f"{s}={v:.3f}" for s, v in zip(SPECIES_ORDER, y.to_numpy(float))))
    if w1000.exists():
        w = pd.read_csv(w1000, sep="\t", index_col=0)
        print("\n  старое окно ±1000 (скрипт 62):")
        one("55 генов", w, core)
        g47 = [g for g in core if g in set(pd.read_csv(old_path, sep="\t", index_col=0).index)] if old_path.exists() else [g for g in core if g in uf]
        one("те же 47", w, g47)


def per_library(species: str, up: int, down: int, mode: str) -> pd.DataFrame:
    d = pd.read_csv(s62.EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    d = annotate_table(d, peak_to_symbol_de(RUSLAN / SPECIES[species]["de"]))
    fp, src = loci_of(species, d, mode)

    cols = lib_cols(d, "fib") + lib_cols(d, "epi")
    num = d[cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    num = num.where(num >= THR, 0.0)
    pool = d.assign(mid=(d["start"] + d["end"]) / 2.0)
    pool = pool.assign(**{c: num[c] for c in cols})
    pool = pool[num.to_numpy().sum(axis=1) > 0].reset_index(drop=True)

    owner, how = assign_peaks(fp, pool, up, down)
    pool = pool.assign(owner=owner).dropna(subset=["owner"])
    m = pool.groupby("owner", observed=True)[cols].sum()
    m = m[m.to_numpy().sum(axis=1) > 0]
    m.columns = [f"{species}|{'fib' if c.startswith('norm.C') else 'epi'}|{c.split('.')[1]}"
                 for c in cols]
    m.index.name = None
    layers = ", ".join(f"{k} {v}" for k, v in src.items())
    print(f"[{species:8s}] локусов {len(fp):>6} ({layers}); "
          f"якорь {how['якорь']} / чужое имя {how['чужое_имя']} / nearest {how['nearest']}; "
          f"генов с сигналом {len(m):>6}", flush=True)
    return np.log10(m + 1.0)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--up", type=int, default=5000)
    ap.add_argument("--down", type=int, default=1000)
    ap.add_argument("--loci", choices=["annot_first", "common_only", "common_then_annot"],
                    default="common_then_annot")
    ap.add_argument("--tag-suffix", default="anchored",
                    help="appended to the output tag so the previous nearest-only "
                         "matrices are kept")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"{a.loci}_u{a.up}d{a.down}"
    if a.tag_suffix:
        tag = f"{tag}_{a.tag_suffix}"

    per = {sp: per_library(sp, a.up, a.down, a.loci) for sp in SPECIES_ORDER}
    for sp, m in per.items():
        m.to_csv(OUT / f"M_{sp}_{tag}.tsv", sep="\t")

    def universe(tissue: str) -> list[str]:
        sets = []
        for sp in SPECIES_ORDER:
            c = [x for x in per[sp].columns if f"|{tissue}|" in x]
            m = per[sp][c]
            sets.append(set(m.index[(m > 0).any(axis=1)]))
        return sorted(set.intersection(*sets))

    uf, ue = universe("fib"), universe("epi")
    ua = sorted(set(uf) | set(ue))
    print(f"\nвселенная: фибробласты {len(uf)}, эпибласт {len(ue)}, объединение {len(ua)}")
    for nm, u in (("fib", uf), ("epi", ue), ("all", ua)):
        (OUT / f"universe_{nm}_{tag}.txt").write_text("\n".join(u) + "\n")

    logm = np.log10([SPECIES[s]["mass_g"] for s in SPECIES_ORDER])
    nz = [int((per[sp][[x for x in per[sp].columns if "|fib|" in x]] > 0).any(axis=1).sum())
          for sp in SPECIES_ORDER]
    r = stats.linregress(logm, np.array(nz, float))
    print("обнаружено в фибробластах: " + ", ".join(f"{s} {n}" for s, n in zip(SPECIES_ORDER, nz)))
    print(f"перекос обнаружения по массе: r={r.rvalue:+.3f} p={r.pvalue:.3f}")

    def combine(tissue: str | None, genes: list[str]) -> pd.DataFrame:
        parts = []
        for sp in SPECIES_ORDER:
            c = [x for x in per[sp].columns if tissue is None or f"|{tissue}|" in x]
            parts.append(per[sp].reindex(genes)[c])
        return pd.concat(parts, axis=1)

    mats = {"fib": combine("fib", uf), "epi": combine("epi", ue), "all": combine(None, ua)}
    for nm, m in mats.items():
        m = m.dropna(how="all")
        m.to_csv(OUT / f"M_{nm}_{tag}.tsv", sep="\t")
        print(f"M_{nm}: {m.shape[0]} генов x {m.shape[1]} библиотек, "
              f"пропусков {int(m.isna().sum().sum())}")

    meta = [dict(sample=c, species=c.split("|")[0], tissue=c.split("|")[1],
                 library=c.split("|")[2], mass_g=SPECIES[c.split("|")[0]]["mass_g"],
                 log10_mass=np.log10(SPECIES[c.split("|")[0]]["mass_g"]))
            for c in mats["all"].columns]
    pd.DataFrame(meta).to_csv(OUT / f"samples_{tag}.tsv", sep="\t", index=False)

    ox = RES / "working_sets_round2" / "working_kegg_00190_w1000.txt"
    if ox.exists():
        core = [l.strip() for l in ox.read_text().splitlines() if l.strip()]
        n_in = sum(1 for g in core if g in uf)
        print(f"\nрабочий набор OXPHOS: {n_in} из {len(core)} "
              f"в новой фибробластной вселенной")
        report_oxphos(mats["fib"], core, uf, logm, tag)
    print(f"\nwrote {OUT}/*_{tag}.tsv")


if __name__ == "__main__":
    main()
