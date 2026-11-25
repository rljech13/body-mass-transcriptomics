#!/usr/bin/env python3
"""Which genes track body mass, without asking a pathway map first.

Two questions.

A. Unbiased scan. Over the genes detected in all seven birds, correlate the
   fibroblast CAGE level with log10 adult body mass. With n = 7 and thousands of
   genes, |r| > 0.75 is p < 0.05 by itself, so the honest null is an EXACT
   permutation of the mass labels: all 7! = 5040 orderings, of which 2520 are
   distinct up to sign. For every threshold we compare the observed number of
   genes against that null, and we also do BH over genes.

B. Avian thermogenesis. Birds have no functional UCP1 and no brown fat, so the
   KEGG hsa04714 map is mammal-centric by construction. The avian mechanism runs
   through SERCA and sarcolipin, plus the uncoupling-protein paralogues birds do
   have and the thyroid axis. That set is tested explicitly with the same
   matched null used for the KEGG sets.
"""
from __future__ import annotations

import importlib
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
s53 = importlib.import_module("53_three_checks")
s57 = importlib.import_module("57_working_sets_round2")
from oxphos_thread_lib import RES, SPECIES, SPECIES_ORDER, load_geneset, normalize_symbol  # noqa: E402

OUT = RES / "unbiased_scan"
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], float))

AVIAN_THERMO = {
    # SERCA / sarcolipin futile calcium cycling — the avian non-shivering route
    "ATP2A1", "ATP2A2", "ATP2A3", "SLN", "PLN", "RYR1", "RYR2", "RYR3",
    "CASQ1", "CASQ2", "MYOZ1", "TRDN",
    # uncoupling-protein paralogues birds do have (no functional UCP1)
    "UCP2", "UCP3", "SLC25A27", "SLC25A14",
    # adenine nucleotide translocase, the other proposed leak route
    "SLC25A4", "SLC25A5", "SLC25A6",
    # thyroid axis, the main systemic driver of avian heat production
    "DIO1", "DIO2", "DIO3", "THRA", "THRB", "TRH", "TSHB",
    # avian-relevant lipid handling in muscle
    "PPARA", "PPARD", "CPT1A", "CPT2", "ACADM", "ACADVL",
}


def matrices():
    f = OUT / "cage_fib.tsv"
    g = OUT / "cage_epi.tsv"
    if f.exists() and g.exists():
        return pd.read_csv(f, sep="\t", index_col=0), pd.read_csv(g, sep="\t", index_col=0)
    mats = {sp: s57.measure(sp, 1000) for sp in SPECIES_ORDER}
    fib = pd.DataFrame({sp: mats[sp]["fib"] for sp in SPECIES_ORDER})
    epi = pd.DataFrame({sp: mats[sp]["epi"] for sp in SPECIES_ORDER})
    OUT.mkdir(parents=True, exist_ok=True)
    fib.to_csv(f, sep="\t")
    epi.to_csv(g, sep="\t")
    return fib, epi


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fib, epi = matrices()
    idx = sorted(set(fib.index) | set(epi.index))
    fib = fib.reindex(idx)
    epi = epi.reindex(idx)
    det = (fib.notna() & (fib > 0)).all(axis=1) & (epi.notna() & (epi > 0)).all(axis=1)
    genes = np.asarray(idx)[det.to_numpy()]
    X = np.log10(fib.loc[genes].to_numpy(float))
    print(f"вселенная: {len(genes)} генов, детектированных у всех семи птиц в обеих тканях")

    def rvec(order: np.ndarray) -> np.ndarray:
        y = LOGM[order]
        yc = y - y.mean()
        Xc = X - X.mean(axis=1, keepdims=True)
        num = Xc @ yc
        den = np.sqrt((Xc**2).sum(axis=1) * (yc**2).sum())
        return num / np.where(den == 0, np.nan, den)

    r_obs = rvec(np.arange(7))
    p_obs = 2 * stats.t.sf(np.abs(r_obs) * np.sqrt(5 / np.clip(1 - r_obs**2, 1e-12, None)), df=5)
    bh = stats.false_discovery_control(np.nan_to_num(p_obs, nan=1.0), method="bh")

    perms = list(itertools.permutations(range(7)))
    thr = [0.80, 0.85, 0.90, 0.95]
    obs_cnt = {t: int((np.abs(r_obs) >= t).sum()) for t in thr}
    null = {t: [] for t in thr}
    for pm in perms:
        rv = rvec(np.array(pm))
        for t in thr:
            null[t].append(int((np.abs(rv) >= t).sum()))
    print("\nA. Беспристрастный поиск: сколько генов превышают порог |r|, против точной перестановки массы")
    print(f"{'|r| >=':>8} {'наблюдено':>10} {'медиана нуля':>13} {'95-й перц. нуля':>16} {'доля перестановок >= набл.':>28}")
    rows = []
    for t in thr:
        nl = np.array(null[t])
        pe = float((nl >= obs_cnt[t]).mean())
        print(f"{t:>8.2f} {obs_cnt[t]:>10d} {np.median(nl):>13.0f} {np.percentile(nl, 95):>16.0f} {pe:>28.4f}")
        rows.append(dict(порог=t, наблюдено=obs_cnt[t], медиана_нуля=float(np.median(nl)),
                         перц95_нуля=float(np.percentile(nl, 95)), p_перестановки=pe))
    pd.DataFrame(rows).to_csv(OUT / "threshold_vs_null.tsv", sep="\t", index=False)
    print(f"\nпроходят BH 0.05 по генам: {int((bh < 0.05).sum())}")
    print(f"минимальная p по генам: {np.nanmin(p_obs):.2e}")

    tab = pd.DataFrame({"gene": genes, "r": r_obs.round(3), "p": p_obs, "bh": bh.round(3)})
    tab = tab.sort_values("r")
    tab.to_csv(OUT / "per_gene_mass_fit.tsv", sep="\t", index=False)
    print("\n25 самых отрицательных генов:")
    print(", ".join(f"{g}({r:+.2f})" for g, r in zip(tab["gene"].head(25), tab["r"].head(25))))
    print("\n25 самых положительных:")
    t2 = tab.sort_values("r", ascending=False)
    print(", ".join(f"{g}({r:+.2f})" for g, r in zip(t2["gene"].head(25), t2["r"].head(25))))

    ox = {normalize_symbol(g) for g in load_geneset("kegg_00190")} - {""}
    inuni = set(genes)
    print("\nB. Наборы на той же вселенной, средний уровень против массы")
    for nm, s in (("OXPHOS hsa00190", ox),
                  ("птичий термогенез (SERCA/UCP-паралоги/тиреоид)", AVIAN_THERMO)):
        present = sorted(s & inuni)
        if len(present) < 5:
            print(f"  {nm}: в вселенной только {len(present)} генов — не считаю")
            continue
        for layer, w in (("фибробласты", fib), ("эпибласт", epi)):
            r, p = s53.fit(np.log10(w.reindex(present).mean(axis=0).to_numpy(float)))
            pn, nm2 = s53.matched_null_p(w, present, r)
            print(f"  {nm:48s} {layer:12s} n={len(present):>3} r={r:+.3f} p={p:.3f} против нуля "
                  f"{'—' if np.isnan(pn) else round(pn, 4)}")
        print(f"     состав: {', '.join(present)}")
    print(f"\nwrote {OUT / 'per_gene_mass_fit.tsv'}")


if __name__ == "__main__":
    main()
