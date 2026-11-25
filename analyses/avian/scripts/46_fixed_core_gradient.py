#!/usr/bin/env python3
"""Fixed-core test on the colleagues' 133-gene OXPHOS list.

Decisions taken BEFORE the run (20.09.2026), so that no threshold is chosen
after seeing the answer:

  primary list  = his oxphos_genes, 133 names
  "present in all" = a gene counts for a species only if it is detected in that
                     species in ALL THREE layers at once (CAGE fib, CAGE epi,
                     GEO RNA). Strict core = genes that clear this in 7/7 species.
  reported      = the whole breadth ladder k = 7,6,5,4,3, not the single k=7
                  point, because restricting to the common core is itself a
                  biased choice: it discards exactly the small, highly expressed
                  subunits that Augustus annotations lose.

With the gene list frozen the count term is a constant, so log10 SUM and
log10 MEAN differ by a constant and give an identical correlation. One number
per layer per k, plus an empirical p against nulls matched on breadth and on
chicken abundance decile.

Sensitivity lists are computed alongside but are not the headline.
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    RES,
    SPECIES,
    SPECIES_ORDER,
    load_geneset,
    normalize_symbol,
    nuclear_oxphos,
)

CACHE = RES / "reverify_confounds" / "gene_matrices.pkl"
OUT = RES / "fixed_core_gradient"
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], float))
LAYERS = {"CAGE_fib": "CAGE фибробласты", "CAGE_epi": "CAGE эпибласт", "GEO_RNA": "RNA-seq фибробласты"}
PRIMARY = "его OXPHOS 133"
N_NULL = 2000
RNG = np.random.default_rng(20260920)


def fit(y: np.ndarray) -> tuple[float, float]:
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan
    res = stats.linregress(LOGM[ok], y[ok])
    return float(res.rvalue), float(res.pvalue)


def norm(it) -> list[str]:
    return sorted({normalize_symbol(g) for g in it} - {""})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with CACHE.open("rb") as fh:
        layers = pickle.load(fh)
    wide = {lay: pd.DataFrame({sp: layers[lay][sp][0] for sp in SPECIES_ORDER}) for lay in LAYERS}

    sets = {
        PRIMARY: norm(load_geneset("ruslan_oxphos_133")),
        "KEGG hsa00190": norm(load_geneset("kegg_00190")),
        "наш ядерный OXPHOS 77": norm(nuclear_oxphos()),
        "KEGG hsa04714": norm(load_geneset("kegg_04714")),
        "контроль: гликолиз": norm(load_geneset("GLYCOLYSIS")),
    }

    # breadth under the chosen rule: detected in a species in all three layers
    ok_all = None
    for lay in LAYERS:
        m = wide[lay]
        good = m.notna() & (m > 0)
        good = good.reindex(index=sorted(set().union(*[wide[l].index for l in LAYERS])), fill_value=False)
        ok_all = good if ok_all is None else (ok_all & good)
    breadth = ok_all.sum(axis=1)

    # abundance reference for null matching: chicken CAGE fibroblasts
    chick = wide["CAGE_fib"]["Chicken"].reindex(breadth.index)
    pool = breadth.index[breadth > 0]
    decile = pd.qcut(chick.reindex(pool).rank(method="first"), 10, labels=False, duplicates="drop")
    key = pd.DataFrame({"breadth": breadth.reindex(pool), "decile": decile}, index=pool)

    def null_p(lay: str, genes: list[str], obs_r: float) -> tuple[float, float]:
        hit = [g for g in genes if g in key.index]
        if len(hit) < 5 or not np.isfinite(obs_r):
            return np.nan, np.nan
        want = key.loc[hit].groupby(["breadth", "decile"]).size()
        buckets = {
            k: key.index[(key["breadth"] == k[0]) & (key["decile"] == k[1])].to_numpy() for k in want.index
        }
        m = wide[lay]
        null = np.empty(N_NULL)
        for i in range(N_NULL):
            pick: list[str] = []
            for k, need in want.items():
                av = buckets[k]
                if len(av):
                    pick.extend(RNG.choice(av, size=min(need, len(av)), replace=False))
            null[i] = fit(np.log10(m.reindex(pick).mean(axis=0).to_numpy(float)))[0]
        null = null[np.isfinite(null)]
        return float((np.abs(null) >= abs(obs_r)).mean()), float(np.median(null))

    rows = []
    for name, genes in sets.items():
        present = [g for g in genes if g in breadth.index]
        for k in (7, 6, 5, 4, 3):
            core = [g for g in present if breadth.get(g, 0) >= k]
            for lay, lab in LAYERS.items():
                m = wide[lay]
                vals = np.log10(m.reindex(core).mean(axis=0).to_numpy(float))
                r, p = fit(vals)
                pn, nm = (np.nan, np.nan)
                if name == PRIMARY:
                    pn, nm = null_p(lay, core, r)
                rows.append(
                    dict(
                        gene_set=name,
                        n_list=len(genes),
                        k=k,
                        n_core=len(core),
                        layer=lab,
                        r=None if np.isnan(r) else round(r, 3),
                        p=None if np.isnan(p) else round(p, 4),
                        p_matched_null=None if np.isnan(pn) else round(pn, 4),
                        null_median_r=None if np.isnan(nm) else round(nm, 3),
                    )
                )
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "fixed_core_gradient.tsv", sep="\t", index=False)

    strict = [g for g in sets[PRIMARY] if breadth.get(g, 0) >= 7]
    lost = [g for g in sets[PRIMARY] if breadth.get(g, 0) < 7]
    pd.DataFrame(
        {"gene": sets[PRIMARY], "breadth_all_layers": [int(breadth.get(g, 0)) for g in sets[PRIMARY]]}
    ).sort_values("breadth_all_layers", ascending=False).to_csv(
        OUT / "his133_breadth.tsv", sep="\t", index=False
    )

    pd.set_option("display.width", 220)
    print(f"Его список: {len(sets[PRIMARY])} имён. Строгое ядро 7/7 во всех трёх слоях: {len(strict)}")
    print(f"В ядро вошли: {', '.join(strict)}")
    print(f"\nВыпало {len(lost)} имён; первые 40: {', '.join(lost[:40])}")

    print("\n" + "=" * 110)
    print(f"ГЛАВНЫЙ ТЕСТ — {PRIMARY}, лестница по широте (k = у скольких видов ген есть во всех трёх слоях)")
    print("=" * 110)
    prim = tab[tab.gene_set == PRIMARY].drop(columns=["gene_set", "n_list"])
    for lab in LAYERS.values():
        print(f"\n{lab}")
        print(prim[prim.layer == lab].drop(columns=["layer"]).to_string(index=False))

    print("\n" + "=" * 110)
    print("ЧУВСТВИТЕЛЬНОСТЬ — другие списки, только CAGE фибробласты")
    print("=" * 110)
    sens = tab[(tab.layer == "CAGE фибробласты") & (tab.gene_set != PRIMARY)]
    print(sens.drop(columns=["layer", "p_matched_null", "null_median_r"]).to_string(index=False))

    print(f"\nwrote {OUT / 'fixed_core_gradient.tsv'}")
    print(f"wrote {OUT / 'his133_breadth.tsv'}")


if __name__ == "__main__":
    main()
