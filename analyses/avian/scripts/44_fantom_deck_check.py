#!/usr/bin/env python3
"""Check the colleagues' FANTOM figure (corr - FANTOM - oxphos.pdf).

Their deck plots log10 of the SUMMED OXPHOS TPM per species against log10 body
mass for three FANTOM5 cultured cell types, on four species only (mouse,
chicken, dog, rat), and reports:

    Aortic      cor. -0.894  p = 0.106
    Hepato      cor. -0.601  p = 0.399
    Mesenchymal cor. -0.769  p = 0.231

This script reuses the per-library table already produced by script 38
(results/oxphos_species_thread/fantom_msc_oxphos/libraries.tsv), rebuilds the
summed metric as n_genes_any x mean_max, and asks three questions:

  1. does their four-species fit reproduce on our assembly;
  2. what happens when human is put back in;
  3. how the sum splits into the gene-count term and the level term.

Nothing is refetched; any disagreement with script 38 is a disagreement in the
reading of the same numbers.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import RES  # noqa: E402

SRC = RES / "fantom_msc_oxphos" / "libraries.tsv"
OUT = RES / "fantom_msc_oxphos"

# script 38 masses. Dog: FANTOM donors are beagle (~10 kg), not AnAge 40 kg.
MASS = {"mouse": 20.5, "rat": 300.0, "chicken": 912.5, "dog": 10000.0, "human": 62035.0}
DECK_SPECIES = ["mouse", "chicken", "dog", "rat"]
DECK_REPORTED = {"smc": (-0.894, 0.106), "hep": (-0.601, 0.399), "msc": (-0.769, 0.231)}
CELL_LABEL = {"smc": "Aortic", "hep": "Hepato", "msc": "Mesenchymal"}


def fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 3:
        return float("nan"), float("nan")
    r, p = stats.pearsonr(x, y)
    return float(r), float(p)


def main() -> None:
    libs = pd.read_csv(SRC, sep="\t")
    # mean_sum (per-gene sum over its peaks, averaged over detected genes) is the
    # only level column present for both the peak tables and the gene tables;
    # mean_max is absent for the gene-level species, so it cannot be compared.
    libs["set_sum"] = libs["n_genes_any"] * libs["mean_sum"]

    per_species = (
        libs.groupby(["cell", "species"])
        .agg(
            n_libs=("col", "size"),
            source=("source", lambda s: ";".join(sorted(set(s)))),
            n_genes=("n_genes_any", "mean"),
            mean_level=("mean_sum", "mean"),
            set_sum=("set_sum", "mean"),
        )
        .reset_index()
    )
    per_species["mass_g"] = per_species["species"].map(MASS)
    per_species["log10_mass"] = np.log10(per_species["mass_g"])
    per_species = per_species.sort_values(["cell", "mass_g"])
    per_species.to_csv(OUT / "deck_check_species.tsv", sep="\t", index=False)

    rows = []
    for cell in ["smc", "hep", "msc"]:
        sub = per_species[per_species["cell"] == cell]
        for label, spp in [
            ("deck_4spp", DECK_SPECIES),
            ("plus_human", DECK_SPECIES + ["human"]),
        ]:
            s = sub[sub["species"].isin(spp)]
            if len(s) < 3:
                continue
            x = s["log10_mass"].to_numpy()
            for term, col in [
                ("sum", "set_sum"),
                ("count", "n_genes"),
                ("level", "mean_level"),
            ]:
                y = np.log10(s[col].to_numpy()) if term == "sum" else s[col].to_numpy()
                r, p = fit(x, y)
                rows.append(
                    {
                        "cell": CELL_LABEL[cell],
                        "subset": label,
                        "n": len(s),
                        "term": term,
                        "r": round(r, 3),
                        "p": round(p, 3),
                        "species": ",".join(s.sort_values("mass_g")["species"]),
                    }
                )
    fits = pd.DataFrame(rows)
    fits.to_csv(OUT / "deck_check_fits.tsv", sep="\t", index=False)

    loo = []
    for cell in ["smc", "hep", "msc"]:
        sub = per_species[
            (per_species["cell"] == cell) & (per_species["species"].isin(DECK_SPECIES))
        ]
        if len(sub) < 4:
            continue
        for drop in DECK_SPECIES:
            s = sub[sub["species"] != drop]
            r, p = fit(s["log10_mass"].to_numpy(), np.log10(s["set_sum"].to_numpy()))
            loo.append(
                {"cell": CELL_LABEL[cell], "dropped": drop, "n": len(s), "r": round(r, 3), "p": round(p, 3)}
            )
    pd.DataFrame(loo).to_csv(OUT / "deck_check_loo.tsv", sep="\t", index=False)

    print("=== per-species, summed OXPHOS metric (n_genes x mean_max) ===")
    cols = ["cell", "species", "mass_g", "n_libs", "source", "n_genes", "mean_level", "set_sum"]
    print(per_species[cols].to_string(index=False, float_format=lambda v: f"{v:,.1f}"))

    print("\n=== fits vs log10 mass; deck reported in brackets ===")
    for cell in ["smc", "hep", "msc"]:
        lab = CELL_LABEL[cell]
        dr, dp = DECK_REPORTED[cell]
        print(f"\n{lab}  deck: r={dr} p={dp} on 4 species")
        print(fits[fits["cell"] == lab].to_string(index=False))

    print("\n=== leave-one-out on their four species, summed metric ===")
    print(pd.DataFrame(loo).to_string(index=False))

    print(f"\nwrote {OUT / 'deck_check_species.tsv'}")
    print(f"wrote {OUT / 'deck_check_fits.tsv'}")
    print(f"wrote {OUT / 'deck_check_loo.tsv'}")


if __name__ == "__main__":
    main()
