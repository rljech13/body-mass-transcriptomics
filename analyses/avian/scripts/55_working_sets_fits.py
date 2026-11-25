#!/usr/bin/env python3
"""The working sets and their numbers.

Working set, fixed by the user 20.09.2026: after gene recovery, the genes of the
official KEGG map that are detected in EVERY ONE of the seven birds.

    hsa00190 Oxidative phosphorylation -> 53 genes
    hsa04714 Thermogenesis             -> 71 genes
    GLYCOLYSIS (control)               -> whatever passes the same filter

With the list frozen this way the count term is a constant by construction, so
log10 SUM and log10 MEAN give the same correlation and there is nothing left for
annotation completeness to do.

Reported per set: CAGE fibroblasts, CAGE epiblast, and the GEO RNA-seq of the
same fibroblasts on the common within-table scale. Each with an empirical p
against 2000 nulls matched on detection breadth and on chicken abundance decile.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
s53 = importlib.import_module("53_three_checks")
from oxphos_thread_lib import (  # noqa: E402
    RES,
    SPECIES_ORDER,
    load_geneset,
    load_geo_rnaseq,
    normalize_symbol,
)

OUT = RES / "working_sets"
SETS = {"hsa00190 OXPHOS": "kegg_00190", "hsa04714 Thermogenesis": "kegg_04714",
        "GLYCOLYSIS (контроль)": "GLYCOLYSIS"}


def rna_common_scale() -> pd.DataFrame:
    cols = {}
    for sp in SPECIES_ORDER:
        g = load_geo_rnaseq(sp)
        m = g["gene_mat"].div(g["lib_total"], axis=1) * 1e6
        cols[sp] = m.mean(axis=1)
    return pd.DataFrame(cols)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mats = {sp: s53.measure(sp, promoter_only=False, use_miniprot=True) for sp in SPECIES_ORDER}
    cage = {t: pd.DataFrame({sp: mats[sp][t] for sp in SPECIES_ORDER}) for t in ("fib", "epi")}
    idx = sorted(set(cage["fib"].index) | set(cage["epi"].index))
    ok = None
    for t in ("fib", "epi"):
        m = cage[t].reindex(idx)
        g = m.notna() & (m > 0)
        ok = g if ok is None else (ok & g)
    rna = rna_common_scale()

    rows = []
    for name, key in SETS.items():
        genes = sorted({normalize_symbol(g) for g in load_geneset(key)} - {""})
        det = ok.reindex(genes).fillna(False).infer_objects(copy=False)
        core = sorted(det.index[det.all(axis=1)])
        (OUT / f"working_{key}.txt").write_text("\n".join(core) + "\n")
        print(f"\n{name}: карта {len(genes)}, рабочий набор {len(core)}")
        for layer, w in (("CAGE фибробласты", cage["fib"]), ("CAGE эпибласт", cage["epi"]),
                         ("RNA-seq фибробласты", rna)):
            present = [g for g in core if g in w.index]
            vals = np.log10(w.reindex(present).mean(axis=0).to_numpy(float))
            r, p = s53.fit(vals)
            pn, nm = s53.matched_null_p(w, present, r)
            n_det = w.reindex(present).notna().sum(axis=0).to_numpy(int)
            rows.append(dict(
                набор=name, n_карта=len(genes), n_рабочий=len(core), слой=layer,
                n_измерено=";".join(map(str, n_det)),
                r=round(r, 3), p=round(p, 4),
                p_против_нуля=None if np.isnan(pn) else round(pn, 4),
                медиана_r_нуля=None if np.isnan(nm) else round(nm, 3),
            ))
            print(f"  {layer:22s} r={r:+.3f} p={p:.4f}  против нуля {pn if np.isnan(pn) else round(pn,4)}")
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "working_sets_fits.tsv", sep="\t", index=False)
    pd.set_option("display.width", 240)
    print("\n" + tab.to_string(index=False))
    print(f"\nwrote {OUT / 'working_sets_fits.tsv'}")


if __name__ == "__main__":
    main()
