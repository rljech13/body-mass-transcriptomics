#!/usr/bin/env python3
"""Step one of the story: the two gene sets and their seven-bird cores.

Two sets only, both official KEGG, unedited:
    hsa00190 Oxidative phosphorylation   136 names
    hsa04714 Thermogenesis               232 names   (104 shared with 00190)

For each set, three counts are reported per species and for the panel:
    in map          - names in the KEGG map
    detected        - names with CAGE signal in that bird
    core            - names detected in ALL SEVEN birds  <- what the fits use

Two variants:
    аннотация  - the gene position comes only from the annotation shipped with
                 the peak tables (Augustus for Peacock, Turkey, Emu)
    после      - positions of the genes the annotation lacks are taken from the
                 miniprot loci found in the source assemblies (scripts 48, 50)

Measurement is identical for all seven species in both variants: peaks within
1000 bp of the gene's 5' end on the gene's strand, each peak thresholded at
1 TPM per library, summed, averaged over the tissue's libraries.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

s53 = importlib.import_module("53_three_checks")
from oxphos_thread_lib import RES, SPECIES_ORDER, load_geneset, normalize_symbol  # noqa: E402

OUT = RES / "bird_core_sets"
SETS = {"hsa00190 OXPHOS": "kegg_00190", "hsa04714 Thermogenesis": "kegg_04714"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sets = {k: sorted({normalize_symbol(g) for g in load_geneset(v)} - {""}) for k, v in SETS.items()}

    variants = {}
    for tag, use_miniprot in (("аннотация", False), ("после восстановления", True)):
        s53.USE_MINIPROT = use_miniprot
        mats = {sp: s53.measure(sp, promoter_only=False, use_miniprot=use_miniprot) for sp in SPECIES_ORDER}
        w = {t: pd.DataFrame({sp: mats[sp][t] for sp in SPECIES_ORDER}) for t in ("fib", "epi")}
        idx = sorted(set(w["fib"].index) | set(w["epi"].index))
        ok = None
        for t in ("fib", "epi"):
            m = w[t].reindex(idx)
            g = m.notna() & (m > 0)
            ok = g if ok is None else (ok & g)
        variants[tag] = ok

    rows = []
    for name, genes in sets.items():
        for tag, ok in variants.items():
            det = ok.reindex(genes).fillna(False)
            per = det.sum(axis=0)
            core = sorted(det.index[det.all(axis=1)])
            rows.append(dict(
                набор=name, n_в_карте=len(genes), вариант=tag,
                **{sp: int(per[sp]) for sp in SPECIES_ORDER},
                хотя_бы_у_одного=int(det.any(axis=1).sum()),
                ядро_7_из_7=len(core),
            ))
            if tag.startswith("после"):
                (OUT / f"core7_{name.split()[0]}.txt").write_text("\n".join(core) + "\n")
                print(f"\n=== {name}: ядро 7/7 после восстановления, {len(core)} генов ===")
                print(", ".join(core))
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "bird_core_sets.tsv", sep="\t", index=False)
    pd.set_option("display.width", 240)
    print("\n" + tab.to_string(index=False))
    print(f"\nwrote {OUT / 'bird_core_sets.tsv'}")


if __name__ == "__main__":
    main()
