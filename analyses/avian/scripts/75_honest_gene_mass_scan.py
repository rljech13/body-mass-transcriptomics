#!/usr/bin/env python3
"""Per-gene correlation with adult mass under the homologous-CAGE rule.

Same assignment as the 19-gene rescue (nearest CDS ≤ 2 kb, chicken protein
only, peak goes to the closer 5' if two genes compete). A gene enters the
scan only if that rule gives a fibroblast CAGE number in all seven species.
Empty cells are not filled with zero and are not averaged over.

Null: exact permutation of the seven mass labels (7! = 5040). BH over genes.
"""
from __future__ import annotations

import importlib
import itertools
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
s74 = importlib.import_module("74_nineteen_via_cage")
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    RES,
    SPECIES,
    SPECIES_ORDER,
    load_geneset,
    normalize_symbol,
)

ASM = Path("/path/to/storage/kleiber_birds_gemcat/ref/source_assemblies")
OUT = RES / "honest_gene_scan"
MAX_BP = 2000
LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in SPECIES_ORDER], float))
NINETEEN = list(s74.NINETEEN)
SIXTEEN = [
    "ATP5F1E", "ATP6V0E1", "ATP6V0E2", "COX6A1", "COX6C", "NDUFA5",
    "NDUFA7", "NDUFAB1", "NDUFB2", "NDUFB5", "NDUFB6", "NDUFB8",
    "NDUFC2", "SDHD", "UQCR10", "UQCRC1",
]


def contig(s: str) -> str:
    return s74.contig(s)


def all_chicken_loci(gff: Path) -> pd.DataFrame:
    """One locus per symbol: chicken_* query with highest Identity, then score."""
    best: dict[str, tuple[float, float, dict]] = {}
    with gff.open() as fh:
        for line in fh:
            if "\tmRNA\t" not in line:
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
            if not sym:
                continue
            ident = re.search(r"Identity=([0-9.]+)", f[8])
            idv = float(ident.group(1)) if ident else 0.0
            sc = float(f[5]) if f[5] != "." else 0.0
            st = f[6]
            rec = dict(
                symbol=sym,
                seqname=contig(f[0]),
                strand=st,
                five=float(f[3] if st == "+" else f[4]),
                ident=idv,
                score=sc,
            )
            prev = best.get(sym)
            if prev is None or (idv, sc) > (prev[0], prev[1]):
                best[sym] = (idv, sc, rec)
    if not best:
        return pd.DataFrame(columns=["symbol", "seqname", "strand", "five", "ident", "score"])
    return pd.DataFrame([v[2] for v in best.values()])


def assign_species(peaks: pd.DataFrame, loci: pd.DataFrame) -> pd.DataFrame:
    """Sum fibroblast TPM of peaks whose midpoint is closer to this CDS 5'
    than to any other, and within MAX_BP."""
    if loci.empty or peaks.empty:
        return pd.DataFrame(columns=["gene", "tpm", "n_peaks", "dist", "peak_names", "frac_foreign"])
    index: dict[tuple[str, str], tuple[np.ndarray, list[str]]] = {}
    for (seq, st), g in loci.groupby(["seqname", "strand"], sort=False):
        g = g.sort_values("five")
        index[(seq, st)] = (g["five"].to_numpy(float), g["symbol"].tolist())

    acc: dict[str, list] = {}
    for r in peaks.itertuples():
        key = (r.seq, r.strand)
        packed = index.get(key)
        if packed is None:
            continue
        fives, syms = packed
        mid = float(r.mid)
        i = int(np.searchsorted(fives, mid))
        cands: list[int] = []
        if i > 0:
            cands.append(i - 1)
        if i < len(fives):
            cands.append(i)
        if not cands:
            continue
        ds = [abs(mid - fives[j]) for j in cands]
        jbest = cands[int(np.argmin(ds))]
        d = abs(mid - fives[jbest])
        if d > MAX_BP:
            continue
        if len(cands) == 2 and abs(ds[0] - ds[1]) < 1e-6:
            continue
        gene = syms[jbest]
        nom = r.symbol if isinstance(r.symbol, str) else ""
        acc.setdefault(gene, []).append((float(r.fib), d, nom))

    rows = []
    for gene, hits in acc.items():
        tpm = float(sum(h[0] for h in hits))
        if tpm <= 0:
            continue
        names = sorted({h[2] for h in hits if h[2]})
        foreign = sum(h[0] for h in hits if h[2] and h[2] != gene)
        rows.append(dict(
            gene=gene,
            tpm=tpm,
            n_peaks=len(hits),
            dist=min(h[1] for h in hits),
            peak_names=",".join(names[:6]),
            frac_foreign=foreign / tpm,
        ))
    return pd.DataFrame(rows)


def r_matrix(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    yc = y - y.mean()
    Xc = X - X.mean(axis=1, keepdims=True)
    num = Xc @ yc
    den = np.sqrt((Xc ** 2).sum(axis=1) * (yc ** 2).sum())
    return num / np.where(den == 0, np.nan, den)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ox = {normalize_symbol(g) for g in load_geneset("kegg_00190")} - {""}
    print("локусы куриного белка + CAGE ≤ 2 т.п.н.", flush=True)
    wide = {}
    meta_rows = []
    ident_wide = {}
    for sp in SPECIES_ORDER:
        print(f"  {sp}", flush=True)
        loci = all_chicken_loci(ASM / sp / "reanno_named.gff")
        peaks = s74.load_peaks(sp)
        assigned = assign_species(peaks, loci)
        wide[sp] = assigned.set_index("gene")["tpm"] if len(assigned) else pd.Series(dtype=float)
        ident_wide[sp] = loci.set_index("symbol")["ident"] if len(loci) else pd.Series(dtype=float)
        if len(assigned):
            a = assigned.copy()
            a["species"] = sp
            meta_rows.append(a)
        print(f"    локусов {len(loci)}  генов с CAGE {len(assigned)}", flush=True)

    M = pd.DataFrame(wide).reindex(columns=SPECIES_ORDER)
    I = pd.DataFrame(ident_wide).reindex(index=M.index, columns=SPECIES_ORDER)
    M.to_csv(OUT / "M_fib_nearest2kb.tsv", sep="\t")
    pd.concat(meta_rows, ignore_index=True).to_csv(OUT / "peak_assignment.tsv", sep="\t", index=False)

    det = M.notna() & (M > 0)
    n7 = det.all(axis=1)
    genes = M.index[n7].tolist()
    print(f"\nво всех семи с CAGE: {len(genes)}  (локусов хотя бы у одного: {int(M.notna().any(axis=1).sum())})")

    # check the 19 against the previous file
    prev = RES / "nineteen_via_cage" / "M_nineteen_nearest_cds.tsv"
    if prev.exists():
        P = pd.read_csv(prev, sep="\t", index_col=0)
        print("сверка с 19 nearest_cds:")
        for g in NINETEEN:
            if g not in M.index:
                print(f"  {g}: нет в скане")
                continue
            a = M.loc[g].to_numpy(float)
            b = P.loc[g].reindex(SPECIES_ORDER).to_numpy(float) if g in P.index else np.full(7, np.nan)
            ok = np.isfinite(a) & np.isfinite(b) & (b > 0)
            if ok.any():
                rel = np.nanmax(np.abs(a[ok] - b[ok]) / np.clip(b[ok], 1e-6, None))
            else:
                rel = np.nan
            print(f"  {g}: max отн. разница {rel:.4f}  скан={np.array2string(np.nan_to_num(a), precision=0)}")

    X = np.log10(M.loc[genes, SPECIES_ORDER].to_numpy(float))
    r_obs = r_matrix(X, LOGM)
    p_obs = 2 * stats.t.sf(np.abs(r_obs) * np.sqrt(5 / np.clip(1 - r_obs ** 2, 1e-12, None)), df=5)
    bh = stats.false_discovery_control(np.nan_to_num(p_obs, nan=1.0), method="bh")

    orders = list(itertools.permutations(range(7)))
    n_ge = np.zeros(len(genes), dtype=int)
    thr = (0.70, 0.75, 0.80, 0.85, 0.90)
    null_counts = {t: [] for t in thr}
    for perm in orders:
        rp = r_matrix(X, LOGM[list(perm)])
        n_ge += (np.abs(rp) >= np.abs(r_obs)).astype(int)
        for t in thr:
            null_counts[t].append(int((np.abs(rp) >= t).sum()))
    p_perm = n_ge / len(orders)

    min_ident = I.loc[genes].min(axis=1)
    foreign = (
        pd.concat(meta_rows, ignore_index=True)
        .groupby("gene")["frac_foreign"].mean()
        .reindex(genes)
    )
    fit = pd.DataFrame({
        "gene": genes,
        "r": r_obs,
        "p": p_obs,
        "bh": bh,
        "p_perm": p_perm,
        "min_ident": min_ident.to_numpy(float),
        "frac_foreign": foreign.fillna(0).to_numpy(float),
        "in_hsa00190": [g in ox for g in genes],
        "in_16": [g in SIXTEEN for g in genes],
        "in_19": [g in NINETEEN for g in genes],
        **{sp: M.loc[genes, sp].to_numpy(float) for sp in SPECIES_ORDER},
    }).sort_values("r")
    fit["rank"] = np.arange(1, len(fit) + 1)
    fit.to_csv(OUT / "per_gene_mass_fit.tsv", sep="\t", index=False)

    rows = []
    for t in thr:
        obs = int((np.abs(r_obs) >= t).sum())
        nc = np.array(null_counts[t], float)
        # two-sided count of permutations with at least as many |r|>=t
        pn = float((nc >= obs).mean())
        rows.append(dict(
            порог=t, наблюдено=obs,
            медиана_нуля=float(np.median(nc)),
            перц95_нуля=float(np.percentile(nc, 95)),
            p_перестановки=pn,
        ))
    pd.DataFrame(rows).to_csv(OUT / "threshold_vs_null.tsv", sep="\t", index=False)

    print(f"\nгенов {len(genes)}  p<0.05 {(p_obs < 0.05).sum()}  BH<0.05 {(bh < 0.05).sum()}")
    print("перестановка |r|:")
    for rec in rows:
        print(f"  |r|≥{rec['порог']}  наблюдано {rec['наблюдено']}  "
              f"нуль {rec['медиана_нуля']:.0f} (95% {rec['перц95_нуля']:.0f})  p={rec['p_перестановки']:.3f}")

    print("\n20 самых отрицательных r:")
    show = fit.head(20)
    for r in show.itertuples():
        flag = []
        if r.in_16:
            flag.append("16")
        elif r.in_hsa00190:
            flag.append("00190")
        if r.frac_foreign >= 0.5:
            flag.append("чужое_имя")
        print(f"  {int(r.rank):4d} {r.gene:12s} r={r.r:+.3f} p={r.p:.4f} bh={r.bh:.3f}  "
              f"{','.join(flag) if flag else ''}")

    print("\nгде 16:")
    sub = fit[fit.in_16]
    for r in sub.itertuples():
        print(f"  rank {int(r.rank):4d} {r.gene:10s} r={r.r:+.3f} p={r.p:.4f} bh={r.bh:.3f}  "
              f"чужое_имя {r.frac_foreign:.2f}")

    y16 = np.log10(M.reindex(SIXTEEN)[SPECIES_ORDER].mean(axis=0).to_numpy(float))
    rr = stats.linregress(LOGM, y16)
    print(f"\n16 как набор: r={rr.rvalue:+.3f} p={rr.pvalue:.4f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
