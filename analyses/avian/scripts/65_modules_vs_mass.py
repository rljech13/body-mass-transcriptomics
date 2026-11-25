#!/usr/bin/env python3
"""Modules against body mass, and where the 55 OXPHOS genes landed.

Module-trait correlation is the standard WGCNA one: the module eigengene against
the trait over samples, Pearson, BH across modules. Two versions of the trait:

  по библиотекам  eigengene vs log10 mass over all libraries (42 in fib, 23 in
                  epi). This is what the WGCNA tutorial does, but the libraries
                  of one species share a mass value, so the effective number of
                  independent points is 7, not 42 — the p is optimistic and is
                  reported as such.
  по видам        eigengene averaged per species first, then correlated over the
                  seven species. Comparable with every other number in this
                  project.

Then the actual question: which modules contain the 55 genes of the working
OXPHOS set, and what else is in those modules. For every such module the other
members are listed with their own per-gene correlation with mass, and the module
is tested as a set against breadth-and-abundance matched nulls.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
s53 = importlib.import_module("53_three_checks")
from oxphos_thread_lib import RES, SPECIES_ORDER  # noqa: E402

IN = RES / "wgcna_input"
WG = RES / "wgcna_out"
OUT = RES / "modules_vs_mass"
W = 1000
OXFILE = RES / "working_sets_round2" / "working_kegg_00190_w1000.txt"


def bh(p: np.ndarray) -> np.ndarray:
    return stats.false_discovery_control(np.nan_to_num(p, nan=1.0), method="bh")


def load(net: str):
    gm = pd.read_csv(WG / f"gene_modules_{net}.tsv", sep="\t")
    me = pd.read_csv(WG / f"module_eigengenes_{net}.tsv", sep="\t").set_index("sample")
    return gm, me


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(IN / f"samples_w{W}.tsv", sep="\t").set_index("sample")
    ox = [l.strip() for l in OXFILE.read_text().splitlines() if l.strip()]
    print(f"рабочий набор OXPHOS: {len(ox)} генов\n")

    summary = []
    for net in ("fib", "epi", "all_tc"):
        gm, me = load(net)
        s = samples.loc[me.index]
        mass_lib = s["log10_mass"].to_numpy(float)
        rows = []
        for mod in me.columns:
            v = me[mod].to_numpy(float)
            r1, p1 = stats.pearsonr(mass_lib, v)
            per = pd.DataFrame({"sp": s["species"].to_numpy(), "v": v}).groupby("sp")["v"].mean()
            per = per.reindex(SPECIES_ORDER)
            r2, p2 = s53.fit(per.to_numpy(float))
            name = mod.replace("ME", "")
            n_in = int((gm["module"] == name).sum())
            n_ox = int(((gm["module"] == name) & gm["gene"].isin(ox)).sum())
            rows.append(dict(сеть=net, модуль=name, n_генов=n_in, n_OXPHOS=n_ox,
                             r_по_библиотекам=round(r1, 3), p_по_библиотекам=p1,
                             r_по_видам=round(r2, 3), p_по_видам=p2))
        t = pd.DataFrame(rows)
        t["bh_по_библиотекам"] = bh(t["p_по_библиотекам"].to_numpy()).round(4)
        t["bh_по_видам"] = bh(t["p_по_видам"].to_numpy()).round(4)
        t["p_по_библиотекам"] = t["p_по_библиотекам"].round(5)
        t["p_по_видам"] = t["p_по_видам"].round(4)
        t = t.sort_values("r_по_видам")
        t.to_csv(OUT / f"module_trait_{net}.tsv", sep="\t", index=False)
        summary.append(t)
        print("=" * 104)
        print(f"сеть {net}: модулей {len(t)}, из них с генами OXPHOS {int((t['n_OXPHOS'] > 0).sum())}")
        print("=" * 104)
        print(t[t["n_OXPHOS"] > 0].to_string(index=False))
        sig = t[(t["bh_по_видам"] < 0.05) & (t["модуль"] != "grey")]
        print(f"\nпроходят BH 0.05 по видам: {len(sig)}"
              + ("" if sig.empty else "\n" + sig[["модуль", "n_генов", "n_OXPHOS", "r_по_видам", "bh_по_видам"]].to_string(index=False)))

    # где лежат 55 генов в фибробластной сети и что рядом
    gm, _ = load("fib")
    fib = pd.read_csv(IN / f"M_fib_w{W}.tsv", sep="\t", index_col=0)
    sp_of = pd.read_csv(IN / f"samples_w{W}.tsv", sep="\t").set_index("sample").loc[fib.columns, "species"]
    permean = fib.T.groupby(sp_of.to_numpy()).mean().T.reindex(columns=SPECIES_ORDER)
    logm = np.log10([__import__("oxphos_thread_lib").SPECIES[s]["mass_g"] for s in SPECIES_ORDER])

    def gene_r(g: str) -> float:
        if g not in permean.index:
            return np.nan
        y = permean.loc[g].to_numpy(float)
        return float(stats.pearsonr(logm, y)[0]) if np.isfinite(y).all() else np.nan

    inmod = gm[gm["gene"].isin(ox)]
    dist = inmod["module"].value_counts()
    print("\n" + "=" * 104)
    print("КУДА ПОПАЛИ 55 ГЕНОВ OXPHOS (сеть фибробластов)")
    print("=" * 104)
    print(dist.to_string())

    det = []
    for mod, k in dist.items():
        members = gm[gm["module"] == mod]["gene"].tolist()
        others = [g for g in members if g not in ox]
        rr = {g: gene_r(g) for g in others}
        oth = pd.Series(rr).dropna().sort_values()
        r_set, p_set = s53.fit(np.log10(permean.reindex(members).mean(axis=0).to_numpy(float)))
        pn, nm = s53.matched_null_p(permean, [g for g in members if g in permean.index], r_set)
        det.append(dict(модуль=mod, n_всего=len(members), n_OXPHOS=k, n_прочих=len(others),
                        r_модуля_по_видам=round(r_set, 3), p=round(p_set, 4),
                        p_против_нуля=None if np.isnan(pn) else round(pn, 4),
                        медиана_r_прочих=round(float(oth.median()), 3) if len(oth) else None))
        print(f"\n--- модуль {mod}: {len(members)} генов, из них OXPHOS {k}, прочих {len(others)} ---")
        print(f"    модуль как набор: r={r_set:+.3f} p={p_set:.4f} против нуля "
              f"{'—' if np.isnan(pn) else round(pn, 4)}")
        if len(oth):
            print(f"    прочие члены, 20 самых отрицательных по массе: "
                  + ", ".join(f"{g}({v:+.2f})" for g, v in oth.head(20).items()))
            print(f"    медиана r прочих {oth.median():+.3f}; отрицательных {int((oth < 0).sum())} из {len(oth)}")
            pd.DataFrame({"gene": oth.index, "r_mass": oth.values}).to_csv(
                OUT / f"module_{mod}_other_members.tsv", sep="\t", index=False)
    pd.DataFrame(det).to_csv(OUT / "oxphos_modules_detail.tsv", sep="\t", index=False)
    print(f"\nwrote {OUT}/module_trait_*.tsv, oxphos_modules_detail.tsv")


if __name__ == "__main__":
    main()
