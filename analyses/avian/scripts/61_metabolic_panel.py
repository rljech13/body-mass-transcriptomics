#!/usr/bin/env python3
"""Does anything else scale with body mass, or only the respiratory chain?

This is the direct test of the colleagues' reading of the data (deck slide 7:
"all cells know end point size of the adult body and set up metabolism speed at
very early stages"). That claim makes two predictions:

  1. a metabolic set-point should move MANY energy modules together, not one;
  2. it should be visible in both cell types, since both belong to the same animal.

Prediction 2 is already refuted on the frozen core: epiblast is flat while
fibroblasts are not. This script tests prediction 1 by running a panel of
official KEGG maps through the identical machinery: the genes of the map detected
in all seven birds, mean of log level against log10 adult mass, empirical p
against 2000 nulls matched on detection breadth and chicken abundance decile.

Gene lists come from the KEGG REST API and are not edited.
"""
from __future__ import annotations

import importlib
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
s53 = importlib.import_module("53_three_checks")
from oxphos_thread_lib import DATA, RES, SPECIES_ORDER, normalize_symbol  # noqa: E402

SCAN = RES / "unbiased_scan"
OUT = RES / "metabolic_panel"
GSDIR = DATA / "genesets"

MAPS = {
    "hsa00190": "Oxidative phosphorylation",
    "hsa00020": "Цитратный цикл (TCA)",
    "hsa00010": "Гликолиз / глюконеогенез",
    "hsa00620": "Метаболизм пирувата",
    "hsa00030": "Пентозофосфатный путь",
    "hsa00071": "Бета-окисление жирных кислот",
    "hsa00061": "Синтез жирных кислот",
    "hsa01212": "Метаболизм жирных кислот (сводная)",
    "hsa00280": "Деградация BCAA",
    "hsa00230": "Метаболизм пуринов",
    "hsa04146": "Пероксисома",
    "hsa03050": "Протеасома",
    "hsa03010": "Рибосома (цитозольная)",
    "hsa03040": "Сплайсосома",
    "hsa04110": "Клеточный цикл",
    "hsa04150": "Путь mTOR",
    "hsa04152": "Путь AMPK",
    "hsa04140": "Аутофагия",
}


def kegg_genes(mid: str) -> set[str]:
    f = GSDIR / f"kegg_{mid[3:]}.txt"
    if f.exists():
        return {normalize_symbol(l.strip()) for l in f.read_text().splitlines() if l.strip()} - {""}
    url = f"https://rest.kegg.jp/get/{mid}"
    raw = urllib.request.urlopen(url, timeout=60).read().decode()
    time.sleep(0.4)
    genes, inblock = set(), False
    for line in raw.splitlines():
        if line.startswith("GENE"):
            inblock = True
            line = "     " + line[4:]
        elif line and not line.startswith(" "):
            inblock = False
        if not inblock:
            continue
        body = line.strip()
        if not body:
            continue
        parts = body.split(None, 1)
        if len(parts) < 2:
            continue
        sym = parts[1].split(";")[0].strip()
        s = normalize_symbol(sym.split(",")[0].strip())
        if s:
            genes.add(s)
    f.write_text("\n".join(sorted(genes)) + "\n")
    return genes


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fib = pd.read_csv(SCAN / "cage_fib.tsv", sep="\t", index_col=0)
    epi = pd.read_csv(SCAN / "cage_epi.tsv", sep="\t", index_col=0)
    idx = sorted(set(fib.index) | set(epi.index))
    fib, epi = fib.reindex(idx), epi.reindex(idx)
    det = (fib.notna() & (fib > 0)).all(axis=1) & (epi.notna() & (epi > 0)).all(axis=1)
    universe = set(np.asarray(idx)[det.to_numpy()])
    print(f"вселенная: {len(universe)} генов у всех семи птиц в обеих тканях\n")

    rows = []
    for mid, label in MAPS.items():
        g = kegg_genes(mid)
        core = sorted(g & universe)
        if len(core) < 8:
            print(f"{label:38s} ядро {len(core):>3} — мало, пропускаю")
            continue
        out = dict(карта=mid, набор=label, n_в_карте=len(g), n_ядро=len(core))
        for tag, w in (("fib", fib), ("epi", epi)):
            r, p = s53.fit(np.log10(w.reindex(core).mean(axis=0).to_numpy(float)))
            pn, nm = s53.matched_null_p(w, core, r)
            out[f"r_{tag}"] = round(r, 3)
            out[f"p_{tag}"] = round(p, 4)
            out[f"null_{tag}"] = None if np.isnan(pn) else round(pn, 4)
        rows.append(out)
        print(f"{label:38s} ядро {len(core):>3}  фибробласты r={out['r_fib']:+.3f} "
              f"нуль {out['null_fib']}   эпибласт r={out['r_epi']:+.3f} нуль {out['null_epi']}")
    t = pd.DataFrame(rows).sort_values("r_fib")
    t.to_csv(OUT / "metabolic_panel.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("\n" + t.to_string(index=False))
    print(f"\nwrote {OUT / 'metabolic_panel.tsv'}")


if __name__ == "__main__":
    main()
