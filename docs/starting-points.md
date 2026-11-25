# Entry points, inputs and outputs

These commands refer to existing scripts. They are a reading/launch guide after restoring the original inputs and configuring paths, not a newly added runner. Do not interpret the complete chronological script list as one linear pipeline.

## Birds

Scripts remain in `analyses/avian/scripts/`; imports of numbered modules and `oxphos_thread_lib` remain in the same directory. Their shared `RES` is the original `results/oxphos_species_thread` directory.

| Existing stage | Required context | Main saved output |
|---|---|---|
| `57_working_sets_round2.py --window 1000` | Correct species annotations, both gene-recovery rounds, CAGE matrices, KEGG maps | `working_sets_round2/`: complete-case working sets and fits |
| `74_nineteen_via_cage.py` | Fixed19genes, chicken-protein placements, lifted chicken promoter peaks, species CAGE and neighboring-gene constraints | `nineteen_via_cage/chicken_promoters.tsv`, `identity_cards.tsv`, `per_gene_species.tsv`, `M_nineteen_tpm.tsv` |
| `78_holding_cores.py` | Measured common genes, retained maps, WGCNA modules and fixed species masses | `holding_cores/map_cores_vs_mass.tsv`, `map_species_means.tsv`, `module_cores_vs_mass.tsv`, `module_genes.tsv` |
| `94_nine_set_motifs.py` and `94_monalisa_nine.R` | Nine locked gene sets, matched region backgrounds and precomputed JASPAR/HOCOMOCO hits | `story_cre_motifs/nine_sets/`: set/species motif analyses |
| `98_epi_own_promoters.py` | Epiblast's own promoter/background assignments | Epiblast-specific promoter analysis in the original results hierarchy |

The earlier scripts are retained because they document how measurements and diagnoses changed. The corrected19-to16gene branch does not erase the preceding window-based analyses. The scope of each comparison is recorded in its source docstring and outputs.

## Mammalian WGCNA

Scripts remain in `analyses/mammal_wgcna/code/`. Preserve the `code/` level: many scripts resolve their run directory from `__file__.parents[1]`.

1. `10_prepare_full.py`: source panel matrices and metadata → full linear/log-expression inputs and `config/full_inventory.tsv`.
2. `11_wgcna_full.R`: prepared expression inputs and network settings → full-gene WGCNA results/eigengenes/membership.
3. `12_statistics_full.py`: frozen networks, corrected masses and covariance/admission registry → module association and model-comparison tables.
4. `15_wgcna_full_LOO.R`, `16_tree_sensitivity.py`, `19_Hallmark_and_gene_jackknife.py`: separate saved sensitivity/annotation analyses.
5. `module_scores.R`: consumes `data/saved_eigengenes.tsv`, membership, species and historical tests in its original drawing-workspace context; writes descriptive module scores and within-study/all-three-study BH tables. Its single positional argument is that drawing workspace. This is separate from network discovery.

Other existing numbered scripts retain the earlier restricted-gene and native-motif branches. They are described in the branch inventory rather than silently mixed into the full-gene WGCNA sequence.

## Mammalian gene maps

Scripts remain in `analyses/mammal_maps/code/`.

- `00_gene_sets.py` and `01_sets_mass.py`: stored maps and panel-specific measured expression → per-map sums, coverage and mass association results.
- `run_all_maps.py`: existing extended screen, with explicit `--results-root`, `--genesets-root` and `--out` arguments. These point to the original measured panels, fixed map catalogs and a separate output directory. `export_go.R` is part of the original catalog preparation.
- Promoter, enrichment and phylogenetic scripts form a separate downstream branch. Retain their original gene/region background and correction families.

## Mammalian pySCENIC

Scripts remain in `analyses/mammal_scenic/code/`; the run folder contains `inputs/`, `reference/`, `results/`, `config/` and `audit/`.

| Existing stage | Input | Output |
|---|---|---|
| `01_prepare.py` | Ma2016 count/metadata workbooks, fixed technical pairs, one-to-one orthology, HGNC mappings and reference annotation | Collapsed/filtered expression, CPM and mappings in `inputs/`; saved plan/configuration |
| `02_grn.py --seed 20261023` | Prepared expression and TF list | `results/20261023/adjacencies.tsv.gz` and GRN summary |
| `03_build_db.py` | Reference sequence/windows, motif matrices and custom ranking tools | Custom cisTarget resources under `reference/` |
| `04_prune.py --seed 20261023` | Adjacencies, expression, rankings and motif annotation | Pruned regulons, nonself target inventory and `AUCell.tsv` |
| `05_analyze.py` | Saved primary/repeated runs, metadata, selected modules and corrected masses | Species-level AUCell, regulon–mass associations, module overlaps and sensitivity diagnostics |
| `06_verify.R`, `08_stability.py` | Saved analysis tables and repeated seeds | Independent numerical checks and stability summaries |

The existing seed-sensitivity runs used20261024and20261025. The original shell entry points remain alongside these scripts; their environments and working directories must be restored before use.

## Figures

`figures/figure_bank/code/`, `figures/regulon_scatter/code/` and `figures/target_heatmap/code/` retain the original `code/` plus `data/` layout. The compact saved inputs identify exactly what was displayed. Full inferential testing families and upstream model files are not replaced by the displayed subset.
