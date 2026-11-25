#!/usr/bin/env python3
"""One metric, every gene set, all three bird layers.

Built for the MBSJ story, so that section 3 and 4 quote numbers produced by a
single recipe instead of by five different scripts.

Metric, exactly as in the colleagues' figures ('Ruslan metric'):
    per species, per gene of the list -> value in that layer
    SUM   = log10( sum over the genes the species detects )
    MEAN  = log10( mean over the genes the species detects )
    COUNT = how many genes of the list the species detects
Fit is OLS of each quantity on log10 adult body mass, n = 7.

SUM = COUNT + MEAN in log space, so the three rows of each block are a complete
decomposition, not three alternative analyses.

MEAN additionally gets an empirical p against nulls matched on detection breadth
(in how many of the 7 species a gene is seen) and on chicken abundance decile.

Layers come from the cached matrices of script 23, so nothing is requantified.
"""
from __future__ import annotations

import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    DATA,
    RES,
    SPECIES,
    SPECIES_ORDER,
    load_geneset,
    normalize_symbol,
    nuclear_oxphos,
)

CACHE = RES / "reverify_confounds" / "gene_matrices.pkl"
RMD = Path("/path/to/storage/deviatiiarov/chicken/bird_sizes/sizes.Rmd")
OUT = RES / "unified_setscan"
LOGM = np.log10(np.array([SPECIES[sp]["mass_g"] for sp in SPECIES_ORDER], float))
N_NULL = 2000
RNG = np.random.default_rng(20260920)

LAYERS = {"CAGE_fib": "CAGE фибробласты", "CAGE_epi": "CAGE эпибласт", "GEO_RNA": "RNA-seq фибробласты"}


def fit(y: np.ndarray) -> tuple[float, float]:
    y = np.asarray(y, float)
    ok = np.isfinite(y)
    if ok.sum() < 3:
        return np.nan, np.nan
    res = stats.linregress(LOGM[ok], y[ok])
    return float(res.rvalue), float(res.pvalue)


def deck_lists() -> dict[str, list[str]]:
    """mt_genes / cell_cycle_genes as typed in the colleagues' Rmd."""
    out: dict[str, list[str]] = {}
    if not RMD.exists():
        return out
    txt = RMD.read_text(errors="ignore")
    for var, key in (("mt_genes", "их митохондриальный 123"), ("cell_cycle_genes", "их клеточный цикл 116")):
        m = re.search(rf"^{var}\s*<-\s*c\((.*?)\)\s*$", txt, re.S | re.M)
        if not m:
            continue
        toks = re.findall(r"[A-Za-z][A-Za-z0-9\-]*", m.group(1))
        out[key] = sorted({normalize_symbol(t) for t in toks} - {""})
    return out


def core29() -> list[str]:
    t = pd.read_csv(RES / "definitions" / "oxphos_module.tsv", sep="\t")
    return sorted({normalize_symbol(g) for g in t["gene"]} - {""})


def gene_sets() -> dict[str, list[str]]:
    def norm(it) -> list[str]:
        return sorted({normalize_symbol(g) for g in it} - {""})

    sets: dict[str, list[str]] = {
        "его OXPHOS 133": norm(load_geneset("ruslan_oxphos_133")),
        "его термогенез 229": norm(load_geneset("ruslan_thermo_229")),
        "KEGG hsa00190": norm(load_geneset("kegg_00190")),
        "KEGG hsa04714": norm(load_geneset("kegg_04714")),
        "наш ядерный OXPHOS 77": norm(nuclear_oxphos()),
        "наше связное ядро 29": core29(),
        "контроль: гликолиз": norm(load_geneset("GLYCOLYSIS")),
        "контроль: TCA": norm(load_geneset("TCA")),
        "контроль: биогенез митохондрий": norm(load_geneset("MITO_BIOGENESIS")),
    }
    sets.update(deck_lists())
    return sets


def matched_null_p(wide: pd.DataFrame, genes: list[str], observed_r: float) -> tuple[float, float]:
    """Empirical two-sided p for the MEAN metric against breadth+abundance nulls."""
    detected = wide.notna()
    breadth = detected.sum(axis=1)
    chick = wide["Chicken"]
    pool = wide.index[breadth > 0]
    decile = pd.qcut(chick.reindex(pool).rank(method="first"), 10, labels=False, duplicates="drop")
    key = pd.DataFrame({"breadth": breadth.reindex(pool), "decile": decile}, index=pool)
    hit = [g for g in genes if g in wide.index]
    if len(hit) < 5:
        return np.nan, np.nan
    want = key.loc[hit].groupby(["breadth", "decile"]).size()
    buckets = {k: key.index[(key["breadth"] == k[0]) & (key["decile"] == k[1])].to_numpy() for k in want.index}
    null = np.empty(N_NULL)
    for i in range(N_NULL):
        pick: list[str] = []
        for k, need in want.items():
            avail = buckets[k]
            if len(avail) == 0:
                continue
            pick.extend(RNG.choice(avail, size=min(need, len(avail)), replace=False))
        vals = np.log10(wide.loc[pick].mean(axis=0).to_numpy(float))
        null[i] = fit(vals)[0]
    null = null[np.isfinite(null)]
    p = float((np.abs(null) >= abs(observed_r)).mean())
    return p, float(np.median(null))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with CACHE.open("rb") as fh:
        layers = pickle.load(fh)

    wide = {
        lay: pd.DataFrame({sp: layers[lay][sp][0] for sp in SPECIES_ORDER})
        for lay in LAYERS
    }

    sets = gene_sets()
    rows = []
    for name, genes in sets.items():
        for lay, lab in LAYERS.items():
            m = wide[lay]
            hit = m.reindex(genes)
            n_det = hit.notna().sum(axis=0).to_numpy(float)
            mean_det = np.log10(hit.mean(axis=0).to_numpy(float))
            set_sum = np.log10(hit.sum(axis=0).to_numpy(float))
            r_sum, p_sum = fit(set_sum)
            r_mean, p_mean = fit(mean_det)
            r_cnt, p_cnt = fit(n_det)
            p_null, null_med = matched_null_p(m, genes, r_mean)
            rows.append(
                dict(
                    gene_set=name,
                    n_list=len(genes),
                    layer=lab,
                    n_detected=";".join(str(int(v)) for v in n_det),
                    r_sum=round(r_sum, 3),
                    p_sum=round(p_sum, 4),
                    r_mean=round(r_mean, 3),
                    p_mean=round(p_mean, 4),
                    r_count=round(r_cnt, 3),
                    p_count=round(p_cnt, 4),
                    p_mean_matched_null=None if np.isnan(p_null) else round(p_null, 4),
                    null_median_r=None if np.isnan(null_med) else round(null_med, 3),
                )
            )
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "unified_setscan.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    for lab in LAYERS.values():
        print("=" * 120)
        print(lab)
        sub = tab[tab.layer == lab].drop(columns=["layer"])
        print(sub.to_string(index=False))
    print(f"\nwrote {OUT / 'unified_setscan.tsv'}")


if __name__ == "__main__":
    main()
