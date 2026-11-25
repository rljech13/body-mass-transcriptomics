#!/usr/bin/env python3
"""Targeted gene recovery: which OXPHOS / thermogenesis genes are missing where.

Scope fixed by the user (20.09.2026): only the two official KEGG sets,
hsa00190 Oxidative phosphorylation and hsa04714 Thermogenesis. No whole-genome
re-annotation — we hunt the specific names.

A gene is a RECOVERY TARGET for a species if chicken has a CAGE promoter peak
carrying that symbol and the species does not. Chicken is the reference because
its galGal6 annotation is native and the libraries are the same experiment, so
"chicken has it, this species does not" is the cleanest statement of a gap.

Genes absent from chicken too are reported separately: for those the CAGE panel
has no evidence at all, and 13 of them are mtDNA names that are absent from every
species' peak table by construction.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from oxphos_thread_lib import (  # noqa: E402
    RES,
    RUSLAN,
    SPECIES,
    annotate_table,
    load_geneset,
    normalize_symbol,
    peak_to_symbol_de,
)

EXP = Path("/path/to/user-home/chicken/kleiber_birds_gemcat/data/dropbox_20260918/exp_tab")
OUT = RES / "missing_gene_targets"
NATIVE = ["Chicken", "Duck", "Quail", "Zebra"]
AUGUSTUS = ["Turkey", "Peacock", "Emu"]
MTDNA = {
    "ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB",
    "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6",
}
COMPLEX = [
    ("I (NDUF)", lambda g: g.startswith("NDUF")),
    ("II (SDH)", lambda g: g.startswith("SDH")),
    ("III (UQCR/CYC)", lambda g: g.startswith(("UQCR", "CYC"))),
    ("IV (COX)", lambda g: g.startswith("COX")),
    ("V (ATP5)", lambda g: g.startswith("ATP5")),
    ("V-ATPase (ATP6V/ATP4/ATP12)", lambda g: g.startswith(("ATP6V", "ATP6AP", "ATP4", "ATP12", "TCIRG"))),
    ("мтДНК", lambda g: g in MTDNA),
]


def norm(it) -> set[str]:
    return {normalize_symbol(g) for g in it} - {""}


def symbols(species: str) -> set[str]:
    """Canonical symbol assignment: same as annotate_table in the rest of the project.

    geneId / chicken_genes first, then the DE-table rescue by peak id. Skipping the
    rescue undercounts Zebra, Quail and Duck, whose peak tables carry non-symbol
    identifiers for part of the genes.
    """
    d = pd.read_csv(EXP / f"{species}_fib_epi_exp.txt", sep="\t", low_memory=False)
    id2sym = peak_to_symbol_de(RUSLAN / SPECIES[species]["de"])
    d = annotate_table(d, id2sym)
    d = d[d["class"] == "promoter"]
    return set(d["symbol"]) - {""}


def complex_of(g: str) -> str:
    for name, test in COMPLEX:
        if test(g):
            return name
    return "прочее"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sets = {"hsa00190": norm(load_geneset("kegg_00190")), "hsa04714": norm(load_geneset("kegg_04714"))}
    union = sets["hsa00190"] | sets["hsa04714"]
    found = {sp: symbols(sp) for sp in NATIVE + AUGUSTUS}

    rows = []
    for g in sorted(union):
        rows.append(
            dict(
                gene=g,
                complex=complex_of(g),
                in_00190=g in sets["hsa00190"],
                in_04714=g in sets["hsa04714"],
                **{sp: (g in found[sp]) for sp in NATIVE + AUGUSTUS},
                n_native=sum(g in found[sp] for sp in NATIVE),
            )
        )
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "gene_presence_matrix.tsv", sep="\t", index=False)

    print("Присутствие символа среди пиков класса promoter\n")
    print(f"{'вид':10s} {'hsa00190':>12s} {'hsa04714':>12s}")
    for sp in NATIVE + AUGUSTUS:
        a = len(sets["hsa00190"] & found[sp])
        b = len(sets["hsa04714"] & found[sp])
        mark = "  (Augustus)" if sp in AUGUSTUS else ""
        print(f"{sp:10s} {a:>6d}/{len(sets['hsa00190']):<5d} {b:>6d}/{len(sets['hsa04714']):<5d}{mark}")

    print("\n" + "=" * 100)
    print("ЦЕЛИ ВОССТАНОВЛЕНИЯ: курица ген имеет, этот вид — нет")
    print("=" * 100)
    targets = {}
    for sp in AUGUSTUS:
        t = sorted((found["Chicken"] & union) - found[sp])
        targets[sp] = t
        by_c = pd.Series([complex_of(g) for g in t]).value_counts()
        print(f"\n--- {sp}: {len(t)} генов ---")
        print("  по комплексам: " + ", ".join(f"{k} {v}" for k, v in by_c.items()))
        for name, _ in COMPLEX + [("прочее", None)]:
            sub = [g for g in t if complex_of(g) == name]
            if sub:
                print(f"    {name:30s} {', '.join(sub)}")

    common = sorted(set.intersection(*[set(targets[sp]) for sp in AUGUSTUS]))
    print(f"\nПропали у всех трёх сразу: {len(common)}")
    print("  " + ", ".join(common))

    pd.DataFrame(
        [{"species": sp, "n_targets": len(t), "genes": ";".join(t)} for sp, t in targets.items()]
    ).to_csv(OUT / "recovery_targets.tsv", sep="\t", index=False)
    Path(OUT / "targets_common_all3.txt").write_text("\n".join(common) + "\n")

    absent_ch = sorted(union - found["Chicken"])
    mt = [g for g in absent_ch if g in MTDNA]
    print(f"\nНет и у курицы: {len(absent_ch)} имён, из них мтДНК {len(mt)}")
    print("  не мтДНК: " + ", ".join(g for g in absent_ch if g not in MTDNA))

    print(f"\nwrote {OUT / 'gene_presence_matrix.tsv'}")
    print(f"wrote {OUT / 'recovery_targets.tsv'}")
    print(f"wrote {OUT / 'targets_common_all3.txt'}")


if __name__ == "__main__":
    main()
