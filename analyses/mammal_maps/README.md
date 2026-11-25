# analyses/mammal_maps

Existing scripts; original names and analytical logic retained. See the [analysis map](../../docs/analysis-map.md) and [input/configuration guide](../../docs/inputs-and-configuration.md).

| File | Source purpose / first documentation line |
|---|---|
| `code/00_gene_sets.py` | Fix Ruslan's gene sets: his mitochondrial list (142, from his workspace) and the KEGG human pathways he used |
| `code/01_sets_mass.py` | Step 1 of Ruslan's approach in mammalian fibroblasts: his gene sets versus adult mass. |
| `code/02_fetch_genomes.sh` | Download genome FASTA + GFF3 for each species in config/genomes.tsv from NCBI Datasets; a species is DONE when unzip -t passes. |
| `code/03_promoters_fimo.py` | Step 2a: native promoter windows and FIMO motif presence per species (plan_step2.json). |
| `code/03b_run_all.sh` | Run 03_promoters_fimo.py for every species once its genome zip is verified; up to 28 in parallel; then step 2b. |
| `code/04_enrichment_mass.py` | Step 2b: per species x set x motif Mantel-Haenszel log odds ratio (GC quintiles, Haldane 0.5), |
| `code/05_gsea_phylo.py` | Step 3: genome-wide phylogenetic GSEA in three mammalian fibroblast panels with cross-panel concordance (plan_step3.json). |
| `code/common_rm.py` | Panel loader shared with 01_sets_mass.py (copied verbatim) plus the fixed per-panel background used in step 1. |
| `code/describe_correlations.py` | Requested descriptive interpretation; retains phylogenetic results separately. |
| `code/export_go.R` | Original analysis script; inspect configuration and read/write calls. |
| `code/run_all_maps.py` | Whole measured-map linear sums, no trait-driven membership selection. |

Inputs and outputs retain the original research-run hierarchy. No new automatic runner was added. Machine-specific prefixes were replaced by explicit placeholders; configure these before executing an original stage.
