#!/usr/bin/env python3
"""Find a size-linked gene set on honest promoters, confirmed across tissues.

Genes are ranked by correlation of log10(promoter TPM) with log10(mass) in
one tissue. The candidate set is the k strongest genes in one direction.
The sum of that set is then scored in the other tissue: y = log10(sum TPM)
versus log10(mass). The null repeats the whole selection under every
permutation of the seven masses, so the confirmation tissue is not used to
choose the genes and the p-value includes the search.

Primary test, fixed before looking: k = 16, decline, select on fibroblast,
confirm on epiblast. The other k and the rise are the same procedure,
reported beside it.
"""

from __future__ import annotations

import math
from itertools import permutations
from pathlib import Path

import numpy as np
from scipy import stats

BASE = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
ORDER = ["Zebra", "Quail", "Chicken", "Duck", "Peacock", "Turkey", "Emu"]
MASS = np.array([11.93, 96.6, 912.5, 1121.0, 4093.75, 5811.0, 36200.0])
LOGM = np.log10(MASS)
FIB = BASE / "results/oxphos_species_thread/honest_gene_scan/M_fib_nearest2kb.tsv"
EPI = BASE / "results/oxphos_species_thread/story_epi_traits/M_epi_nearest2kb.tsv"
OUT = BASE / "results/oxphos_species_thread/story_epi_traits/cross_tissue_size_set.tsv"
KS = (16, 50)


def load(path: Path) -> dict[str, np.ndarray]:
    lines = path.read_text().splitlines()
    hdr = lines[0].split("\t")
    cols = [hdr.index(s) for s in ORDER]
    out: dict[str, np.ndarray] = {}
    for line in lines[1:]:
        parts = line.split("\t")
        vals = []
        ok = True
        for i in cols:
            if i >= len(parts) or parts[i] == "":
                ok = False
                break
            v = float(parts[i])
            if not v > 0:
                ok = False
                break
            vals.append(v)
        if ok:
            out.setdefault(parts[0].upper(), np.array(vals, float))
    return out


def gene_r(raw: np.ndarray, logm: np.ndarray) -> np.ndarray:
    x = np.log10(raw)
    xc = x - x.mean(axis=1, keepdims=True)
    c = logm - logm.mean()
    return (xc @ c) / (
        np.sqrt((xc ** 2).sum(axis=1)) * math.sqrt(float(c @ c))
    )


def sum_r(raw: np.ndarray, idx: np.ndarray, logm: np.ndarray) -> float:
    y = np.log10(raw[idx].sum(axis=0))
    return float(stats.pearsonr(logm, y).statistic)


def main() -> None:
    fib = load(FIB)
    epi = load(EPI)
    genes = sorted(set(fib) & set(epi))
    raw_f = np.array([fib[g] for g in genes])
    raw_e = np.array([epi[g] for g in genes])
    r_f = gene_r(raw_f, LOGM)
    r_e = gene_r(raw_e, LOGM)
    rho = float(stats.spearmanr(r_f, r_e).statistic)
    print(f"genes in both tissues, TPM>0 in all 7: {len(genes)}")
    print(f"Spearman of per-gene r, fib vs epi: {rho:+.3f}")

    perms = np.array(list(permutations(range(7))))
    # null sum-r for each (direction, k, select_on_fib)
    # direction -1 = decline (smallest r), +1 = rise
    records = []
    for select_fib in (True, False):
        src = raw_f if select_fib else raw_e
        dst = raw_e if select_fib else raw_f
        src_name = "fib" if select_fib else "epi"
        dst_name = "epi" if select_fib else "fib"
        src_r = r_f if select_fib else r_e
        for sign, tail in ((-1, "decline"), (1, "rise")):
            for k in KS:
                if sign < 0:
                    obs_idx = np.argpartition(src_r, k - 1)[:k]
                else:
                    obs_idx = np.argpartition(-src_r, k - 1)[:k]
                obs = sum_r(dst, obs_idx, LOGM)
                obs_src = sum_r(src, obs_idx, LOGM)
                null = np.empty(len(perms))
                for i, perm in enumerate(perms):
                    shuffled = LOGM[perm]
                    rr = gene_r(src, shuffled)
                    if sign < 0:
                        idx = np.argpartition(rr, k - 1)[:k]
                    else:
                        idx = np.argpartition(-rr, k - 1)[:k]
                    null[i] = sum_r(dst, idx, shuffled)
                if sign < 0:
                    p = float((null <= obs).mean())
                else:
                    p = float((null >= obs).mean())
                names = [genes[i] for i in obs_idx]
                # order names by the selection-tissue r
                names.sort(key=lambda g: src_r[genes.index(g)] * sign, reverse=True)
                records.append(
                    (src_name, dst_name, tail, k, obs_src, obs, p, names)
                )
                print(
                    f"select {src_name:3} {tail:7} k={k:3}  "
                    f"sum r_{src_name}={obs_src:+.3f}  "
                    f"sum r_{dst_name}={obs:+.3f}  null p={p:.4f}"
                )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as handle:
        handle.write(
            "select\tconfirm\ttail\tk\tr_select\tr_confirm\tp_null\tgenes\n"
        )
        for src, dst, tail, k, rs, rc, p, names in records:
            handle.write(
                f"{src}\t{dst}\t{tail}\t{k}\t{rs:.6f}\t{rc:.6f}\t{p:.6f}\t"
                + ",".join(names)
                + "\n"
            )
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
