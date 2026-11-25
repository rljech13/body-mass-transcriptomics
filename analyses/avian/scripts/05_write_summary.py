#!/usr/bin/env python3
"""Compile short results summary answering Q1–Q8 and B1."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
RES = ROOT / "results"
NOTES = ROOT / "notes"
NOTES.mkdir(exist_ok=True)


def main():
    lines = []
    lines.append("# Kleiber / OXPHOS fibroblast analysis — results summary\n\n")
    lines.append("Project: `chicken/kleiber_birds_gemcat`\n\n")

    lines.append("## Clarification: atlas & torpor\n\n")
    lines.append(
        "Atlas and torpor papers were **not rewritten**. Their data (CAGE atlas ESC/fib TPM, CRE passport) were used as inputs only.\n\n"
    )

    # A2 fits
    fits = RES / "allometry_ols_fits.tsv"
    if fits.exists():
        df = pd.read_csv(fits, sep="\t")
        lines.append("## Q1–Q2 Fibroblast OXPHOS ~ mass (RNAseq, 4 symbol-annotated spp)\n\n")
        lines.append("| module | n | slope | r | p | |err| vs -1/3 | |err| vs -1/4 |\n|---|---:|---:|---:|---:|---:|---:|\n")
        for _, r in df.iterrows():
            lines.append(
                f"| {r['module']} | {r['n']} | {r['slope']:.3f} | {r['r']:.3f} | {r['p']:.3f} | {r['abs_err_m13']:.3f} | {r['abs_err_m14']:.3f} |\n"
            )
        ox = df[df.module == "OXPHOS"]
        if len(ox):
            p = float(ox.iloc[0]["p"])
            slope = float(ox.iloc[0]["slope"])
            if p > 0.05:
                verdict = f"**Null / not significant** (slope={slope:.3f}, p={p:.3f}) — consistent with Jimenez fibroblast cellular MR null."
            else:
                verdict = f"**Significant association** (slope={slope:.3f}, p={p:.3f}); compare to mass-specific targets −1/3 and −1/4 (not organismal 2/3)."
            lines.append(f"\n**OXPHOS verdict:** {verdict}\n\n")
        lines.append(
            "Species with portable gene symbols: ZebraFinch, Quail, Duck, Chicken. "
            "Peacock/Turkey/Emu GEO tables use assembly `g####` IDs — excluded from symbol intersection.\n\n"
        )

    # A3
    for tag in ("m13", "m14"):
        p = RES / f"exploratory_best_summary_{tag}.json"
        if p.exists():
            s = json.loads(p.read_text())
            lines.append(f"## Q3 Exploratory gene subset (target {s['target_slope']:.3f}) — POST-HOC\n\n")
            lines.append(
                f"Best slope={s.get('best_slope')}, |err|={s.get('abs_err'):.3f}, n_genes={s.get('n_genes')}. "
                f"**Not confirmatory.** Genes in `results/exploratory_best_genes_{tag}.txt`.\n\n"
            )

    # A4
    mod = RES / "chicken_esc_fib_module_scores.tsv"
    if mod.exists():
        m = pd.read_csv(mod, sep="\t")
        lines.append("## Q4–Q5 Chicken ESC vs fibroblast (+ GEMCAT)\n\n")
        lines.append(m.to_markdown(index=False) + "\n\n")
    ov = RES / "A4_gene_overlap_allometry_vs_escfib.json"
    if ov.exists():
        o = json.loads(ov.read_text())
        lines.append(
            f"Gene overlap (neg allometry ∩ down in ESC): n={o.get('n_overlap_neg_down')}; "
            f"(neg allometry ∩ up in ESC): n={o.get('n_overlap_neg_up')}.\n\n"
        )
    if (RES / "gemcat_esc_over_fib_top50.tsv").exists():
        lines.append("GEMCAT top metabolites: `results/gemcat_esc_over_fib_top50.tsv`.\n\n")

    q67 = RES / "A5_Q6_Q7_multispecies_unavailable.json"
    if q67.exists():
        lines.append("## Q6–Q7 Multi-species ΔOXPHOS ~ mass\n\n")
        lines.append("**Unavailable** — no paired ESC/fib across bird species. Chicken-only ΔOXPHOS reported in JSON.\n\n")

    # A6
    if (RES / "cage_promoter_species_summary.tsv").exists():
        lines.append("## Q8 CAGE promoter regulation\n\n")
        lines.append(pd.read_csv(RES / "cage_promoter_species_summary.tsv", sep="\t").to_markdown(index=False) + "\n\n")
        lines.append("Figures: `figures/fig_A6_cage_promoter_vs_mass.png`.\n\n")

    # B1
    b1 = RES / "B1_gwas_cre_inventory.json"
    if b1.exists():
        d = json.loads(b1.read_text())
        lines.append("## Paper B scout\n\n")
        lines.append(f"**{d['decision']['recommendation']}**\n\n")
        lines.append("See `notes/B1_gwas_cre_inventory.md`.\n\n")

    lines.append("## Key figures\n\n")
    for fig in sorted((ROOT / "figures").glob("*.png")):
        lines.append(f"- `{fig.relative_to(ROOT)}`\n")

    (NOTES / "RESULTS_SUMMARY.md").write_text("".join(lines))
    print("Wrote", NOTES / "RESULTS_SUMMARY.md")


if __name__ == "__main__":
    main()
