#!/usr/bin/env python3
"""Epiblast promoters with the same 2 kb rule, without the fibroblast filter.

Script 80 kept a peak only when fibroblast CAGE was above the threshold and
then read epiblast on that peak. Here a peak enters if its own epiblast mean
is above the threshold. The catchment is script 75: chicken-protein CDS, same
strand, closer 5' end, within 2 kb. Missing gene is dropped. y = log10(sum).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

s74 = importlib.import_module("74_nineteen_via_cage")
s75 = importlib.import_module("75_honest_gene_mass_scan")
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    RES,
    SPECIES,
    SPECIES_ORDER,
    annotate_table,
    lib_cols,
    load_geneset,
    normalize_symbol,
    peak_to_symbol_de,
)

OUT = RES / "epi_own_promoters"
ORDER = SPECIES_ORDER
LOGM = np.log10(np.array([SPECIES[s]["mass_g"] for s in ORDER], float))
SIXTEEN = s75.SIXTEEN


def load_peaks_epi(species: str) -> pd.DataFrame:
    d = pd.read_csv(s74.s62.EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    d = annotate_table(d, peak_to_symbol_de(
        Path("/path/to/storage/deviatiiarov/chicken/bird_sizes") / SPECIES[species]["de"]))
    ecols = lib_cols(d, "epi")
    fcols = lib_cols(d, "fib")
    enum = d[ecols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    enum = enum.where(enum >= s74.THR, 0.0)
    fnum = d[fcols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    fnum = fnum.where(fnum >= s74.THR, 0.0)
    d = d.assign(
        mid=(d["start"] + d["end"]) / 2.0,
        fib=enum.mean(axis=1),
        had_fib=fnum.mean(axis=1) > 0,
        seq=d["seqnames"].map(s74.contig),
    )
    return d[d["fib"] > 0].copy()


def fit_sum(wide: pd.DataFrame, genes: list[str]) -> tuple[int, float, float, list[str]]:
    hit = []
    for g in genes:
        if g not in wide.index:
            continue
        y = wide.loc[g, ORDER].to_numpy(float)
        if np.all(np.isfinite(y)) and np.all(y > 0):
            hit.append(g)
    if len(hit) < 3:
        return len(hit), float("nan"), float("nan"), hit
    total = np.log10(wide.loc[hit, ORDER].sum(axis=0).to_numpy(float))
    r, p = stats.pearsonr(LOGM, total)
    return len(hit), float(r), float(p), hit


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ox = [normalize_symbol(g) for g in load_geneset("kegg_00190")]
    ox = [g for g in dict.fromkeys(ox) if g]
    wide = {}
    only_epi_peaks = {}
    for sp in ORDER:
        print(f"{sp}", flush=True)
        loci = s75.all_chicken_loci(s74.ASM / sp / "reanno_named.gff")
        peaks = load_peaks_epi(sp)
        only_epi_peaks[sp] = int((~peaks["had_fib"]).sum())
        print(f"  epi peaks {len(peaks)}  of which no fibroblast {only_epi_peaks[sp]}  loci {len(loci)}", flush=True)
        assigned = s75.assign_species(peaks, loci)
        wide[sp] = assigned.set_index("gene")["tpm"] if len(assigned) else pd.Series(dtype=float)
        print(f"  genes {len(assigned)}", flush=True)
    M = pd.DataFrame(wide).reindex(columns=ORDER)
    M.to_csv(OUT / "M_epi_own_nearest2kb.tsv", sep="\t")

    old = pd.read_csv(RES / "story_epi_traits" / "M_epi_nearest2kb.tsv", sep="\t", index_col=0)
    old.index = old.index.astype(str).str.upper()
    M.index = M.index.astype(str).str.upper()

    both = [g for g in ox if g in M.index and g in old.index]
    same = 0
    changed = []
    for g in both:
        a = M.loc[g, ORDER].to_numpy(float)
        b = old.loc[g, ORDER].to_numpy(float)
        if np.allclose(np.nan_to_num(a), np.nan_to_num(b), rtol=1e-6, atol=1e-4):
            same += 1
        else:
            changed.append(g)
    new_genes = sorted(set(M.index) - set(old.index))
    new_ox = [g for g in ox if g in M.index and g not in old.index]
    print(f"\nOXPHOS in both matrices {len(both)} identical {same} changed {len(changed)}")
    print("changed", " ".join(changed[:40]))
    print(f"genes only in the epiblast call {len(new_genes)}  of them OXPHOS {len(new_ox)} {' '.join(new_ox)}")
    print("epi-only peaks", only_epi_peaks)

    for label, genes in (
        ("16", SIXTEEN),
        ("00190", ox),
    ):
        n, r, p, hit = fit_sum(M, genes)
        n0, r0, p0, hit0 = fit_sum(old, genes)
        print(f"SUM {label:8s} own-epi n={n:3d} r={r:+.3f} p={p:.3g}   old-epi n={n0:3d} r={r0:+.3f} p={p0:.3g}")
        gained = sorted(set(hit) - set(hit0))
        lost = sorted(set(hit0) - set(hit))
        if gained or lost:
            print(f"  membership gained {gained} lost {lost}")

    # nominal per-gene slopes on the new all-7 OXPHOS list
    n, r, p, hit = fit_sum(M, ox)
    rows = []
    for g in hit:
        y = np.log10(M.loc[g, ORDER].to_numpy(float))
        rr, pp = stats.pearsonr(LOGM, y)
        rows.append((g, rr, pp, g in SIXTEEN))
    tab = pd.DataFrame(rows, columns=["gene", "r", "p", "in16"]).sort_values("p")
    ps = tab["p"].to_numpy(float)
    order = np.argsort(ps)
    q = np.empty(len(ps))
    prev = 1.0
    m = len(ps)
    for rank, i in enumerate(order[::-1]):
        k = m - rank
        prev = min(prev, ps[i] * m / k)
        q[i] = prev
    tab["bh"] = q
    tab.to_csv(OUT / "oxphos_epi_own_per_gene.tsv", sep="\t", index=False)
    nom = tab[tab.p < 0.05]
    print(f"nominal p<0.05 among {len(tab)} all-7 genes: {len(nom)}  BH<0.05 {(tab.bh<0.05).sum()}")
    print(nom.to_string(index=False))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
