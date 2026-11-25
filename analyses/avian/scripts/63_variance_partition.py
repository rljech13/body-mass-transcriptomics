#!/usr/bin/env python3
"""How much of the variance is species, how much tissue, how much left over.

This is the gate for reading every network built later. If species identity eats
most of the variance in the 65-library matrix, then a network over all libraries
mostly recovers species clusters, and only the per-tissue networks (42 and 23
libraries, no tissue contrast inside them) and the consensus run are safe to read
as co-expression.

No VST: the input is log10(sum of peak TPM + 1) as produced by script 62, because
species are being compared and a within-dataset variance stabilisation would be
fitted on the very differences under test.

Decomposition per gene, type-II sums of squares from a two-factor model with
interaction, species x tissue, using plain least squares on indicator matrices —
no statsmodels dependency, which is broken in scanpy_env.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES  # noqa: E402

IN = RES / "wgcna_input"
OUT = RES / "variance_partition"


def design(labels: list[str]) -> np.ndarray:
    lv = sorted(set(labels))
    return np.stack([[1.0 if l == v else 0.0 for l in labels] for v in lv], axis=1)


def ss_explained(Y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """residual sum of squares after projecting Y (genes x samples) on X (samples x k)"""
    Q, _ = np.linalg.qr(X)
    fitted = Y @ Q @ Q.T
    return ((Y - fitted) ** 2).sum(axis=1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=1000)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    M = pd.read_csv(IN / f"M_all_w{a.window}.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(IN / f"samples_w{a.window}.tsv", sep="\t").set_index("sample").loc[M.columns]
    Y = M.to_numpy(float)
    Y = Y - Y.mean(axis=1, keepdims=True)
    n = Y.shape[1]

    sp = design(meta["species"].tolist())
    tis = design(meta["tissue"].tolist())
    inter = design([f"{s}:{t}" for s, t in zip(meta["species"], meta["tissue"])])
    ones = np.ones((n, 1))

    tot = (Y**2).sum(axis=1)
    rss_full = ss_explained(Y, inter)                      # species x tissue cells
    rss_sp = ss_explained(Y, np.hstack([ones, sp]))
    rss_tis = ss_explained(Y, np.hstack([ones, tis]))
    rss_both = ss_explained(Y, np.hstack([ones, sp, tis]))

    frac = pd.DataFrame({
        "gene": M.index,
        "вид": (rss_tis - rss_both) / tot,                 # added by species over tissue
        "ткань": (rss_sp - rss_both) / tot,                # added by tissue over species
        "общая_часть": (tot - rss_both - (rss_tis - rss_both) - (rss_sp - rss_both)) / tot,
        "взаимодействие": (rss_both - rss_full) / tot,
        "остаток": rss_full / tot,
    }).set_index("gene")
    frac.to_csv(OUT / f"per_gene_w{a.window}.tsv", sep="\t")

    print(f"матрица: {M.shape[0]} генов x {n} библиотек, окно {a.window}\n")
    print("доля дисперсии, медиана по генам (и межквартильный размах):")
    for c in ["вид", "ткань", "общая_часть", "взаимодействие", "остаток"]:
        q1, med, q3 = np.percentile(frac[c], [25, 50, 75])
        print(f"  {c:16s} {med:6.3f}   [{q1:.3f}; {q3:.3f}]")

    ox = [l.strip() for l in (RES / "working_sets_round2" / "working_kegg_00190_w1000.txt").read_text().splitlines() if l.strip()]
    inside = frac.reindex([g for g in ox if g in frac.index])
    print(f"\nтот же расклад только на 55 генах OXPHOS (найдено {len(inside)}):")
    for c in ["вид", "ткань", "общая_часть", "взаимодействие", "остаток"]:
        print(f"  {c:16s} {inside[c].median():6.3f}")

    print("\nто же по тканям отдельно (внутри ткани остаётся только вид и остаток):")
    for tname, mfile in (("фибробласты", f"M_fib_w{a.window}.tsv"), ("эпибласт", f"M_epi_w{a.window}.tsv")):
        Mt = pd.read_csv(IN / mfile, sep="\t", index_col=0)
        mt = pd.read_csv(IN / f"samples_w{a.window}.tsv", sep="\t").set_index("sample").loc[Mt.columns]
        Yt = Mt.to_numpy(float)
        Yt = Yt - Yt.mean(axis=1, keepdims=True)
        tt = (Yt**2).sum(axis=1)
        r = ss_explained(Yt, design(mt["species"].tolist()))
        f = pd.Series((tt - r) / tt, index=Mt.index)
        print(f"  {tname:12s} {Mt.shape[1]:>2} библиотек, доля дисперсии на вид: медиана {f.median():.3f} "
              f"[{np.percentile(f,25):.3f}; {np.percentile(f,75):.3f}]")

    print(f"\nwrote {OUT}/per_gene_w{a.window}.tsv")


if __name__ == "__main__":
    main()
