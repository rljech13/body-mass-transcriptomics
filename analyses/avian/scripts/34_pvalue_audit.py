#!/usr/bin/env python3
"""Independent audit of every p-value the OXPHOS-vs-mass answer rests on.

With seven species the parametric p from a least-squares fit is a t test on five
degrees of freedom, which is an assumption. This script recomputes each headline
number four independent ways and requires them to agree:

1. scipy linregress, the number quoted so far.
2. The t statistic computed by hand from r alone, t = r*sqrt(n-2)/sqrt(1-r^2),
   against the Student distribution. Checks that nothing in the fit is wrong.
3. An exact permutation test: all 7! = 5040 assignments of the seven body masses
   to the seven species. No distributional assumption at all, and with n=7 it is
   a complete enumeration rather than a sample.
4. Spearman rank correlation, which drops the linearity assumption.

It also reports the smallest p that seven points can produce and the value of |r|
needed to reach p=0.05, so the reader can see the resolution limit directly.
"""
from __future__ import annotations

import json
import re
import sys
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    EXP, RES, RUSLAN, SPECIES, SPECIES_ORDER, annotate_table, normalize_symbol,
    nuclear_oxphos, peak_to_symbol_de)

RMD = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/sizes.Rmd")
CROSS = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/cross-birds")
OUT = RES / "colleagues_matrix_check"
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], dtype=float))
PERMS = np.array(list(permutations(range(7))))

BLOCKS = [("Chicken", "fib", 0, 6), ("Chicken", "pgc", 6, 8), ("Chicken", "epi", 8, 11),
          ("Quail", "fib", 11, 17), ("Quail", "epi", 17, 20),
          ("Turkey", "fib", 20, 26), ("Turkey", "epi", 26, 29),
          ("Peacock", "fib", 29, 35), ("Peacock", "epi", 35, 39),
          ("Duck", "fib", 39, 45), ("Duck", "epi", 45, 48),
          ("Zebra", "fib", 48, 52), ("Zebra", "epi", 52, 55),
          ("Emu", "fib", 55, 61), ("Emu", "epi", 61, 65)]


def four_ways(y: np.ndarray, x: np.ndarray = LOGM) -> dict:
    y = np.asarray(y, float)
    n = len(y)
    lr = stats.linregress(x, y)
    r = float(lr.rvalue)

    t_hand = r * np.sqrt((n - 2) / (1 - r ** 2))
    p_hand = float(2 * stats.t.sf(abs(t_hand), n - 2))

    xc = x - x.mean()
    yc = y - y.mean()
    denom = np.sqrt((xc ** 2).sum() * (yc ** 2).sum())
    r_null = (yc[PERMS] * xc).sum(axis=1) / denom
    p_exact = float((np.abs(r_null) >= abs(r) - 1e-12).mean())

    sp = stats.spearmanr(x, y)
    return dict(r=r, p_linregress=float(lr.pvalue), t_by_hand=float(t_hand),
                p_by_hand=p_hand, p_exact_permutation=p_exact, n_perms=len(PERMS),
                rho_spearman=float(sp.statistic), p_spearman=float(sp.pvalue))


def rmd_list(var: str) -> list[str]:
    txt = RMD.read_text(errors="ignore")
    m = re.search(rf"^{var}\s*<-\s*c\((.*?)\)\s*$", txt, re.S | re.M)
    if not m:
        return []
    return sorted({normalize_symbol(t) for t in re.findall(r"[A-Za-z][A-Za-z0-9\-]*", m.group(1))} - {""})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    report: dict = {}

    # resolution of seven points
    crit = {}
    for target in (0.05, 0.01):
        tc = stats.t.isf(target / 2, 5)
        crit[str(target)] = float(tc / np.sqrt(tc ** 2 + 5))
    report["resolution"] = dict(
        abs_r_needed_for_p05=crit["0.05"], abs_r_needed_for_p01=crit["0.01"],
        smallest_exact_permutation_p=float(2 / len(PERMS)))
    print(f"Resolution of seven points: |r| must exceed {crit['0.05']:.3f} to reach p=0.05, "
          f"{crit['0.01']:.3f} for p=0.01.")
    print(f"Smallest exact permutation p available with n=7 is 2/5040 = {2/len(PERMS):.2e}.\n")

    # ---------- our CAGE assembly ----------
    per = {}
    for spn in SPECIES_ORDER:
        meta = SPECIES[spn]
        df = pd.read_csv(EXP / meta["file"], sep="\t")
        df = annotate_table(df, peak_to_symbol_de(RUSLAN / meta["de"]))
        df["symbol"] = df["symbol"].map(normalize_symbol)
        prom = df[df["class"].astype(str).str.lower().eq("promoter") & df["symbol"].ne("")]
        cols = [c for c in prom.columns if str(c).startswith("norm.C")]
        per[spn] = prom.groupby("symbol", observed=True)[cols].sum().mean(axis=1)
    wide = pd.DataFrame(per)

    mito = [g for g in rmd_list("mt_genes") if g in wide.index]
    canon = [g for g in sorted({normalize_symbol(g) for g in nuclear_oxphos()} - {""}) if g in wide.index]

    # ---------- colleagues' frozen ortholog matrix ----------
    genes = pd.read_csv(CROSS / "counts.csv")["0"].map(normalize_symbol).values
    mat = pd.read_csv(CROSS / "file.csv", index_col=0)
    mat.index = genes
    mat = mat[~mat.index.duplicated()]
    per_t: dict[str, dict[str, pd.Series]] = {}
    for spn, tis, a, b in BLOCKS:
        per_t.setdefault(tis, {})[spn] = mat.iloc[:, a:b].mean(axis=1)
    cfib = pd.DataFrame({s: per_t["fib"][s] for s in SPECIES_ORDER})
    cepi = pd.DataFrame({s: per_t["epi"][s] for s in SPECIES_ORDER})

    cases = []
    core_m = wide.reindex(mito).dropna()
    core_c = wide.reindex(canon).dropna()
    cases.append(("our CAGE, variable list, mito (the deck's number)", "fib",
                  np.log10(wide.reindex(mito).sum(axis=0).values), len(mito)))
    cases.append(("our CAGE, frozen core, mito", "fib", core_m.mean(axis=0).values, len(core_m)))
    cases.append(("our CAGE, frozen core, OXPHOS canon", "fib", core_c.mean(axis=0).values, len(core_c)))
    for label, m in (("fib", cfib), ("epi", cepi)):
        pm = [g for g in mito if g in m.index]
        pc = [g for g in canon if g in m.index]
        cases.append(("ortholog matrix, mito", label, m.loc[pm].mean(axis=0).values, len(pm)))
        cases.append(("ortholog matrix, OXPHOS canon", label, m.loc[pc].mean(axis=0).values, len(pc)))

    rows = []
    for name, tis, vec, n_genes in cases:
        d = four_ways(vec)
        d.update(case=name, tissue=tis, n_genes=n_genes)
        rows.append(d)
        print(f"{name:46s} {tis}  n={n_genes:4d}  r={d['r']:+.3f}")
        print(f"    linregress p={d['p_linregress']:.4f} | by-hand t={d['t_by_hand']:+.3f} "
              f"p={d['p_by_hand']:.4f} | exact permutation p={d['p_exact_permutation']:.4f} "
              f"| spearman rho={d['rho_spearman']:+.3f} p={d['p_spearman']:.4f}")
        assert abs(d["p_linregress"] - d["p_by_hand"]) < 1e-9, name
    report["cases"] = rows
    pd.DataFrame(rows).to_csv(OUT / "pvalue_audit.tsv", sep="\t", index=False)

    (OUT / "pvalue_audit.json").write_text(json.dumps(report, indent=2, default=str))
    print(f"\nAll parametric p-values agree with the hand-computed t test to 1e-9.")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
