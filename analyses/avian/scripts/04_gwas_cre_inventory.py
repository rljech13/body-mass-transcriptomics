#!/usr/bin/env python3
"""B1: inventory chicken GWAS resources vs own CAGE CRE classes — go/no-go scout."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path("/path/to/user-home/chicken/kleiber_birds_gemcat")
RES = ROOT / "results"
NOTES = ROOT / "notes"
ATLAS_REL = Path("/path/to/user-home/chicken/CAGE_atlases/release")
for d in (RES, NOTES):
    d.mkdir(parents=True, exist_ok=True)


def bed_count(path: Path) -> int | None:
    if not path.exists():
        return None
    if path.suffix == ".gz" or str(path).endswith(".bed.gz") or str(path).endswith(".tsv.gz"):
        r = subprocess.run(["bash", "-lc", f"zcat -f {path} | grep -v '^#' | wc -l"], capture_output=True, text=True)
    else:
        r = subprocess.run(["bash", "-lc", f"grep -v '^#' {path} | wc -l"], capture_output=True, text=True)
    try:
        return int(r.stdout.strip().split()[0])
    except Exception:
        return None


def peek_tsv(path: Path, n: int = 3) -> list[str]:
    if not path.exists():
        return []
    cmd = f"zcat -f {path} | head -n {n}"
    r = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    return r.stdout.splitlines()


def main():
    # --- Local CRE inventory ---
    cre_files = {
        "atlas_peak_passport": ATLAS_REL / "atlas_peak_passport.tsv.gz",
        "atlas_novel_cres_unified": ATLAS_REL / "atlas_novel_cres_unified.tsv.gz",
        "atlas_confident_permissive_bed": ATLAS_REL / "atlas_confident_permissive.bed.gz",
        "atlas_tissue_specific_rescue": ATLAS_REL / "atlas_tissue_specific_rescue.bed.gz",
        "atlas_super_enhancers": ATLAS_REL / "atlas_super_enhancers.bed.gz",
        "eg_links_highconf": ATLAS_REL / "atlas_cage_eg_links_highconf.tsv.gz",
    }
    cre_inventory = []
    for name, path in cre_files.items():
        cre_inventory.append(
            {
                "resource": name,
                "path": str(path),
                "exists": path.exists(),
                "n_rows": bed_count(path) if path.exists() else None,
                "head": peek_tsv(path, 2),
            }
        )
    # class breakdown from passport if available
    passport = ATLAS_REL / "atlas_peak_passport.tsv.gz"
    class_counts = {}
    if passport.exists():
        df = pd.read_csv(passport, sep="\t", nrows=500000)
        # find class-like columns
        class_cols = [c for c in df.columns if "class" in c.lower() or c.lower() in ("pla", "ela", "cre_class", "tss_class", "final_class")]
        for c in class_cols[:5]:
            class_counts[c] = df[c].astype(str).value_counts().head(20).to_dict()
        cre_inventory.append({"resource": "passport_columns", "columns": df.columns.tolist()[:40], "class_counts": class_counts})

    # --- Public chicken GWAS inventory (literature / databases; local files if any) ---
    local_gwas = []
    search_roots = [
        Path("/path/to/user-home/chicken"),
        Path("/path/to/storage"),
        Path("/path/to/storage/deviatiiarov/chicken"),
    ]
    for root in search_roots:
        if not root.exists():
            continue
        r = subprocess.run(
            ["bash", "-lc", f"find {root} -maxdepth 4 \\( -iname '*gwas*' -o -iname '*sumstat*' -o -iname '*qtls*' \\) 2>/dev/null | head -80"],
            capture_output=True,
            text=True,
        )
        for line in r.stdout.splitlines():
            local_gwas.append(line)

    # Curated public resources to chase (not downloaded in this scout)
    public_gwas = [
        {
            "name": "Animal QTLdb — Chicken",
            "url": "https://www.animalgenome.org/cgi-bin/QTLdb/GG/index",
            "traits_of_interest": ["body weight", "growth", "feed conversion / residual feed intake", "disease resistance", "viral"],
            "format": "QTL/SNP coordinates; not always full GWAS sumstats",
            "priority": "high",
        },
        {
            "name": "Chicken GTEx / FarmGTEx GWAS catalog links",
            "url": "https://www.farmgtex.org/",
            "traits_of_interest": ["growth", "metabolic", "immunity"],
            "format": "varies; eQTL + some GWAS",
            "priority": "high",
        },
        {
            "name": "Pan et al. 2023 / related chicken epigenome GWAS enrichment",
            "url": "https://doi.org/10.1038/s41467-023-?",
            "note": "Already did chromHMM CRE × GWAS enrichment — baseline to beat with own CAGE CREs",
            "priority": "must_compare",
        },
        {
            "name": "Published broiler growth / FCR GWAS sumstats (various; e.g. body weight, RFI)",
            "examples": [
                "Body weight / growth rate GWAS in commercial broilers",
                "Residual feed intake / feed efficiency GWAS",
                "Marek's disease / ALV / IBDV resistance GWAS",
            ],
            "priority": "high_if_sumstats_public",
        },
        {
            "name": "GWAS Atlas / EVA / Ensembl variation (chicken)",
            "url": "https://www.ebi.ac.uk/eva/",
            "priority": "medium",
        },
    ]

    # Go/no-go criteria
    decision = {
        "local_CAGE_CRE_ready": any(x.get("exists") for x in cre_inventory if x.get("resource") != "passport_columns"),
        "local_GWAS_sumstats_found": len(local_gwas) > 0,
        "recommendation": None,
        "next_steps": [],
    }
    if decision["local_CAGE_CRE_ready"] and not decision["local_GWAS_sumstats_found"]:
        decision["recommendation"] = "CONDITIONAL GO — CRE layer ready; need to download trait-specific chicken GWAS sumstats (growth / virus / feed efficiency) then test enrichment vs Pan chromHMM."
        decision["next_steps"] = [
            "Download Animal QTLdb chicken SNPs + 1–2 public growth/FCR GWAS tables with coordinates on galGal6/galGal5 (liftOver if needed).",
            "Intersect with passport PLA/ELA / novel STRONG / embryo-biased CREs.",
            "Compare enrichment odds ratios vs Pan chromHMM classes; proceed to paper only if CAGE adds signal.",
        ]
    elif decision["local_CAGE_CRE_ready"] and decision["local_GWAS_sumstats_found"]:
        decision["recommendation"] = "GO — both CRE and some local GWAS-like files found; run formal enrichment next."
        decision["next_steps"] = [
            "Normalize genome builds",
            "Enrichment: trait SNPs in CAGE CRE classes vs matched background",
            "Benchmark against Pan chromHMM",
        ]
    else:
        decision["recommendation"] = "NO-GO for now — CRE or path issues."

    out = {
        "cre_inventory": cre_inventory,
        "local_gwas_paths": local_gwas,
        "public_gwas_targets": public_gwas,
        "decision": decision,
        "traits_priority": ["growth / body weight", "viral disease resistance", "metabolic efficiency / FCR / RFI"],
    }
    with open(RES / "B1_gwas_cre_inventory.json", "w") as f:
        json.dump(out, f, indent=2, default=str)

    # human-readable note
    md = []
    md.append("# B1 GWAS × CAGE CRE inventory (scout)\n")
    md.append(f"**Recommendation:** {decision['recommendation']}\n")
    md.append("\n## Local CRE resources\n")
    for x in cre_inventory:
        if "path" in x:
            md.append(f"- `{x['resource']}`: exists={x.get('exists')} n={x.get('n_rows')}\n")
    md.append("\n## Local GWAS-like files\n")
    if local_gwas:
        for p in local_gwas:
            md.append(f"- `{p}`\n")
    else:
        md.append("- None found under chicken / mnt search roots.\n")
    md.append("\n## Public targets (priority)\n")
    for g in public_gwas:
        md.append(f"- **{g['name']}** (priority={g.get('priority')})\n")
    md.append("\n## Next steps\n")
    for s in decision["next_steps"]:
        md.append(f"1. {s}\n")
    (NOTES / "B1_gwas_cre_inventory.md").write_text("".join(md))
    print("".join(md))
    print("Done B1.")


if __name__ == "__main__":
    main()
