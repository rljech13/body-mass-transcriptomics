#!/usr/bin/env python3
"""Is the RNA-seq layer an independent confirmation of the CAGE result?

The colleagues' metric — sum the module's expression over the genes the species
has, average over that species' libraries — reproduces the decline on CAGE and on
RNA-seq, and two platforms agreeing is normally strong evidence. This script asks
whether the two layers are actually independent where it matters.

They are not necessarily: for Peacock, Turkey and Emu both layers resolve Augustus
model identifiers to gene symbols through the same DE spreadsheet, CAGE through the
peak id column and RNA through transcriptId. If the resulting symbol vocabularies
coincide, the count term of the summed metric is shared between the platforms and
the agreement is guaranteed rather than earned.

Checks:
1. Jaccard overlap of the per-species symbol universe between CAGE and RNA.
2. The colleagues' metric on RNA, decomposed into count and level.
3. The same metric on RNA restricted to genes quantified in all seven species.
4. A direct simulation: degrade the chicken gene universe to the size the Augustus
   species have and see whether the decline appears out of annotation alone.
"""
from __future__ import annotations

import json
import pickle
import re
import sys
from itertools import permutations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    RES, SPECIES, SPECIES_ORDER, normalize_symbol, nuclear_oxphos)

RMD = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/sizes.Rmd")
CACHE = RES / "reverify_confounds" / "gene_matrices.pkl"
OUT = RES / "rna_layer_independence"
FIG = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/figures")
RNG = np.random.default_rng(20260824)
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))
PERMS = np.array(list(permutations(range(7))))
AUG = ["Peacock", "Turkey", "Emu"]


def fit(y):
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan
    lr = stats.linregress(LOGM[ok], y[ok])
    xc, yc = LOGM - LOGM.mean(), y - y.mean()
    den = np.sqrt((xc ** 2).sum() * (yc ** 2).sum())
    rn = (yc[PERMS] * xc).sum(axis=1) / den
    p_exact = float((np.abs(rn) >= abs(lr.rvalue) - 1e-12).mean())
    return float(lr.rvalue), float(lr.pvalue), p_exact


def rmd_list(var: str) -> list[str]:
    m = re.search(rf"^{var}\s*<-\s*c\((.*?)\)\s*$", RMD.read_text(errors="ignore"), re.S | re.M)
    if not m:
        return []
    return sorted({normalize_symbol(t) for t in re.findall(r"[A-Za-z][A-Za-z0-9\-]*", m.group(1))} - {""})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with CACHE.open("rb") as fh:
        layers = pickle.load(fh)
    cage = pd.DataFrame({sp: layers["CAGE_fib"][sp][0] for sp in SPECIES_ORDER})
    rna = pd.DataFrame({sp: layers["GEO_RNA"][sp][0] for sp in SPECIES_ORDER})
    report: dict = {}

    # ---------- 1. do the two platforms share a gene vocabulary ----------
    ov = []
    for sp in SPECIES_ORDER:
        a = set(cage[sp].dropna().index)
        b = set(rna[sp].dropna().index)
        ov.append(dict(species=sp, cage_symbols=len(a), rna_symbols=len(b),
                       shared=len(a & b), jaccard=len(a & b) / len(a | b),
                       rna_is_subset_of_cage=bool(b <= a), cage_is_subset_of_rna=bool(a <= b)))
        print(f"1. {sp:8s} CAGE {len(a):6d}  RNA {len(b):6d}  shared {len(a & b):6d}  "
              f"jaccard {ov[-1]['jaccard']:.3f}  RNA subset of CAGE: {ov[-1]['rna_is_subset_of_cage']}")
    report["vocabulary_overlap"] = ov

    r_c, p_c, _ = fit(np.array([o["cage_symbols"] for o in ov], float))
    r_r, p_r, _ = fit(np.array([o["rna_symbols"] for o in ov], float))
    report["universe_vs_mass"] = dict(cage_r=r_c, cage_p=p_c, rna_r=r_r, rna_p=p_r)
    print(f"\n   symbol universe vs log10 mass: CAGE r={r_c:+.3f} p={p_c:.4f}, "
          f"RNA r={r_r:+.3f} p={p_r:.4f}")

    # ---------- 2-3. the colleagues' metric on each layer, variable and frozen ----------
    mito = rmd_list("mt_genes")
    canon = sorted({normalize_symbol(g) for g in nuclear_oxphos()} - {""})
    cc = rmd_list("cell_cycle_genes")
    rows = []
    for lay, mat in (("CAGE fib", cage), ("RNA", rna)):
        for name, gl in (("mitochondrial (deck)", mito), ("OXPHOS canon", canon),
                         ("cell cycle (deck)", cc)):
            hit = mat.reindex([g for g in gl if g in mat.index])
            if len(hit) < 8:
                continue
            n_det = hit.notna().sum(axis=0).values.astype(float)
            r_s, p_s, pe_s = fit(np.log10(hit.sum(axis=0).values))
            r_n, p_n, _ = fit(np.log10(n_det))
            r_m, p_m, _ = fit(np.log10(hit.mean(axis=0).values))
            core = hit.dropna()
            r_f, p_f, pe_f = fit(np.log10(core.sum(axis=0).values)) if len(core) >= 8 else (np.nan,) * 3
            rows.append(dict(layer=lay, gene_set=name, n_list=len(hit),
                             r_sum=r_s, p_sum=p_s, p_exact_sum=pe_s,
                             r_count=r_n, p_count=p_n, r_mean=r_m, p_mean=p_m,
                             n_core=len(core), r_frozen=r_f, p_frozen=p_f, p_exact_frozen=pe_f,
                             n_detected=n_det.astype(int).tolist()))
            print(f"\n2. {lay:9s} {name:22s} n={len(hit):4d}")
            print(f"   their metric  r={r_s:+.3f} p={p_s:.4f}  = count r={r_n:+.3f} p={p_n:.4f}"
                  f"  +  level r={r_m:+.3f} p={p_m:.4f}")
            print(f"   frozen core   n={len(core):4d}  r={r_f:+.3f} p={p_f:.4f} (exact {pe_f:.4f})")
            print(f"   genes found:  {n_det.astype(int).tolist()}")
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "their_metric_by_layer.tsv", sep="\t", index=False)
    report["their_metric"] = rows

    # ---------- 4. can annotation alone manufacture the decline? ----------
    # Take chicken only. Give every species the chicken numbers, then delete from each
    # species exactly the genes that species actually fails to detect. Body mass now has
    # no effect on expression by construction: every gene carries its chicken value.
    sim_rows = []
    for lay, mat in (("CAGE fib", cage), ("RNA", rna)):
        chick = mat["Chicken"]
        for name, gl in (("mitochondrial (deck)", mito), ("OXPHOS canon", canon)):
            present = [g for g in gl if g in mat.index]
            vals = []
            for sp in SPECIES_ORDER:
                seen = mat.loc[present, sp].notna()
                vals.append(float(chick.reindex(pd.Index(present)[seen.values]).dropna().sum()))
            r, p, pe = fit(np.log10(np.array(vals)))
            sim_rows.append(dict(layer=lay, gene_set=name, r=r, p=p, p_exact=pe,
                                 values=np.round(np.log10(vals), 4).tolist()))
            print(f"\n4. {lay:9s} {name:22s} chicken values, real detection pattern: "
                  f"r={r:+.3f} p={p:.4f} (exact {pe:.4f})")
    pd.DataFrame(sim_rows).to_csv(OUT / "annotation_only_simulation.tsv", sep="\t", index=False)
    report["annotation_only_simulation"] = sim_rows

    # ---------- figure ----------
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

    ax = axes[0]
    x = np.arange(7)
    w = 0.38
    ax.bar(x - w / 2, [o["cage_symbols"] for o in ov], w, color="#1d3557", label="CAGE")
    ax.bar(x + w / 2, [o["rna_symbols"] for o in ov], w, color="#e07a5f", label="RNA-seq")
    ax.set_xticks(x)
    ax.set_xticklabels(SPECIES_ORDER, rotation=45, ha="right", fontsize=8)
    ax.set(ylabel="gene symbols in the layer",
           title="1. the two platforms share one vocabulary\nin exactly the three heaviest species")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    sub = tab[tab.gene_set != "cell cycle (deck)"]
    lab = [f"{r.layer}\n{r.gene_set.split(' (')[0]}" for _, r in sub.iterrows()]
    xx = np.arange(len(sub))
    ax.bar(xx - w / 2, sub.r_sum, w, color="#c1121f", label="their metric, variable list")
    ax.bar(xx + w / 2, sub.r_frozen, w, color="#457b9d", label="same metric, frozen list")
    ax.axhline(0, color="#333", lw=1)
    ax.axhline(-0.754, color="#adb5bd", ls="--", lw=1)
    ax.annotate("p=0.05 threshold at n=7", (len(sub) - 0.5, -0.754), fontsize=7,
                ha="right", va="bottom", color="#6c757d")
    ax.set_xticks(xx)
    ax.set_xticklabels(lab, fontsize=7)
    ax.set(ylim=(-1.05, 1.05), ylabel="r vs log10 body mass",
           title="2. on our RNA assembly their metric\ndoes not reproduce the decline at all")
    ax.legend(fontsize=7, frameon=False)

    ax = axes[2]
    s = pd.DataFrame(sim_rows)
    xx = np.arange(len(s))
    cols = ["#c1121f" if v < -0.754 else "#adb5bd" for v in s.r]
    ax.bar(xx, s.r, color=cols)
    for i, (rv, pv) in enumerate(zip(s.r, s.p)):
        ax.annotate(f"p={pv:.3f}", (i, rv), fontsize=7, ha="center", va="top")
    ax.axhline(-0.754, color="#adb5bd", ls="--", lw=1)
    ax.axhline(0, color="#333", lw=1)
    ax.set_xticks(xx)
    ax.set_xticklabels([f"{r.layer}\n{r.gene_set.split(' (')[0]}" for _, r in s.iterrows()], fontsize=7)
    ax.set(ylim=(-1.05, 0.3), ylabel="r vs log10 body mass",
           title="4. every species given the CHICKEN values,\nonly the detection pattern kept")

    fig.suptitle("Do CAGE and RNA-seq confirm each other, or share one annotation?", fontsize=12)
    for p in (FIG / "fig_rna_layer_independence.png", OUT / "fig_rna_layer_independence.png"):
        fig.savefig(p, dpi=160)
    plt.close()

    (OUT / "summary.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
