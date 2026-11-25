# analyses/avian

Existing scripts; original names and analytical logic retained. See the [analysis map](../../docs/analysis-map.md) and [input/configuration guide](../../docs/inputs-and-configuration.md).

| File | Source purpose / first documentation line |
|---|---|
| `scripts/01_ingest_and_allometry.py` | A1–A3: ingest GSE213253, OXPHOS module allometry, exploratory gene-subset search. |
| `scripts/02_chicken_esc_fib_gemcat.py` | A4–A5: chicken ESC vs fibroblast metabolic modules + GEMCAT on fold-change. |
| `scripts/02b_gemcat_ies_primary_recon_baseline.py` | A5: GEMCAT on ESC/fib — iES1300 primary, Recon3D baseline (recomputed, not Ruslan's outputs). |
| `scripts/03_cage_promoter_regulation.py` | A6: CAGE promoter usage for metabolic genes across species and chicken ESC/fib. |
| `scripts/04_gwas_cre_inventory.py` | B1: inventory chicken GWAS resources vs own CAGE CRE classes — go/no-go scout. |
| `scripts/05_write_summary.py` | Compile short results summary answering Q1–Q8 and B1. |
| `scripts/06_recover3_salmon_allometry.py` | Recompute fibroblast module allometry with n=7 after salmon→chicken quant for Emu/Turkey/Peacock. |
| `scripts/07_gemcat_multispp_ies_vs_recon.py` | Same design as Ruslan (per-species epi/fib FC → GEMCAT → mito summary ~ mass), |
| `scripts/08_paperA_genes_gemcat_oxphos_cage.py` | Paper-A differentiation package (arms 1–4) + corpus-ready tables. |
| `scripts/08b_cage_oxphos_n7.py` | CAGE OXPHOS/mito promoter metrics for all 7 fibroblast species. |
| `scripts/09_alt_scaling_not_kleiber.py` | Paper A — alternative scaling (NOT Kleiber). |
| `scripts/10_pgls_simplicity_S.py` | PGLS + simplicity index S for OXPHOS CAGE promoters (n=7). |
| `scripts/11_salmon_all7_chicken.sh` | Quantify GSE213253 fibroblast RNA-seq (7 spp × 1F+1M) against chicken cDNA symbols. |
| `scripts/11b_aggregate_salmon_all7.py` | Aggregate salmon→chicken quants (7 spp) into gene TPM matrix; OXPHOS coverage vs Ruslan 3055. |
| `scripts/12_promoter_honest_reanalysis.py` | Honest promoter reanalysis (not a claim of 'simplification'). |
| `scripts/13_dropbox_oxphos_promoter_architecture.py` | OXPHOS promoter architecture & regulation from Dropbox CAGE fib+epi tables (n=7). |
| `scripts/14_oxphos_gene_drivers_vs_mass.py` | Which OXPHOS genes/complexes drive fib TPM↓ with mass; multi-promoter context. |
| `scripts/15_ortholog_expand_and_confounders.py` | (a) Expand OXPHOS mapping / fair nuclear-only recompute; (b) confounder residuals. |
| `scripts/16_oxphos_promoter_shape_motifs_bidir.py` | OXPHOS promoter shape: width, GC/motifs (chicken), bidirectionality (n=7). |
| `scripts/17_prep_cage_remap_genomes.sh` | Prepare better genomes (+MT) for CAGE remap of GSE213253 fib libraries. |
| `scripts/18_cage_remap_mt_pipeline.sh` | Remap GSE213253 fib CAGE onto genomes_with_mt (chrM included). |
| `scripts/19_oxphos_species_thread_phase2.py` | Phase 2: fib+epi × set/module/regulon — Y1/Y2/Y3, matched, residual, coherence, coupling. |
| `scripts/19a_oxphos_module_consistency_chicken.py` | Phase 1: CEF coexpression module + regulon-lite (prior TFs × FIMO × coupling). |
| `scripts/20_oxphos_cis_promoter_orthologs_n7.py` | Phase 2 cis: GC / CpG o/e / prior-motif density in OXPHOS TSS windows across n=7. |
| `scripts/21_rnaseq_geo_native_recompute.py` | Recompute fibroblast RNA-seq OXPHOS ~ mass on official GSE213253 tables. |
| `scripts/22_gene_correspondence_and_cage_recompute.py` | Build 7-species gene correspondence and recompute CAGE fib+epi. |
| `scripts/23_reverify_oxphos_mass_confounds.py` | Re-verify the OXPHOS-vs-mass decline against explanations the corpus did not test. |
| `scripts/24_random_geneset_null_and_coverage.py` | Is the OXPHOS-vs-mass decline specific to OXPHOS, or does any gene set do it? |
| `scripts/25_matched_null_oxphos_vs_mass.py` | Matched nulls for the OXPHOS-vs-mass decline. |
| `scripts/26_what_explains_the_cage_decline.py` | What, if anything, explains the surviving CAGE fibroblast OXPHOS decline? |
| `scripts/27_abundance_gradient_and_rank_test.py` | Is the CAGE fibroblast decline a property of OXPHOS or of highly expressed genes? |
| `scripts/28_dropout_mechanism_and_fixed_null.py` | The mechanism behind the CAGE fibroblast OXPHOS decline, and why the null missed it. |
| `scripts/29_colleagues_ortholog_matrix_check.py` | Re-run the OXPHOS-vs-mass question on the colleagues' 1-to-1 ortholog matrix. |
| `scripts/30_reproduce_colleagues_heatmap.py` | Reproduce the colleagues' per-gene trait-correlation heatmap and check its statistics. |
| `scripts/31_reproduce_colleagues_setsum_figure.py` | Reproduce the colleagues' 'corr epi fib' figures and find what drives them. |
| `scripts/32_gene_length_discriminator.py` | Is gene length the thing that decides which gene set 'works' against body mass? |
| `scripts/33_same_metric_fixed_list.py` | The deck's own metric and the deck's own gene lists, on a gene list that cannot move. |
| `scripts/34_pvalue_audit.py` | Independent audit of every p-value the OXPHOS-vs-mass answer rests on. |
| `scripts/35_rna_layer_independence.py` | Is the RNA-seq layer an independent confirmation of the CAGE result? |
| `scripts/36_gse167580_fib_oxphos_module.py` | Fibroblast OXPHOS module fraction vs mass: GSE167580 + GSE213253 birds. |
| `scripts/37_wagner_sf_oxphos.py` | Wagner / Ba skin-fibroblast panel (PRJNA564062 + PRJNA674672 SF only). |
| `scripts/38_fantom_msc_oxphos.py` | FANTOM5 cultured mesenchymal CAGE: OXPHOS mean TPM vs mass. |
| `scripts/39_ruslan_metric_fib_panels.py` | Ruslan metric on every local cultured-fibroblast panel. |
| `scripts/40_ruslan_style_scatters.py` | Ruslan-style log–log scatters: theme_light, lm, species labels, 'cor. ; p ='. |
| `scripts/41_ruslan_kegg_thermogenesis.py` | Ruslan metric on his OXPHOS list and on KEGG thermogenesis (hsa04714). |
| `scripts/42_file_panels_ruslan.py` | Ruslan metric on Ma 2016 and GSE105160 mock fibroblasts. |
| `scripts/43_wagner_kegg_panels.py` | Wagner / Ba skin fibroblasts on official KEGG 00190 and 04714. |
| `scripts/44_fantom_deck_check.py` | Check the colleagues' FANTOM figure (corr - FANTOM - oxphos.pdf). |
| `scripts/45_unified_setscan_birds.py` | One metric, every gene set, all three bird layers. |
| `scripts/46_fixed_core_gradient.py` | Fixed-core test on the colleagues' 133-gene OXPHOS list. |
| `scripts/47_missing_genes_targets.py` | Targeted gene recovery: which OXPHOS / thermogenesis genes are missing where. |
| `scripts/48_fetch_source_assemblies.sh` | Fetch the exact assemblies the colleagues' CAGE peaks were called on. |
| `scripts/50_recover_missing_genes.py` | Find the unannotated copies of the missing OXPHOS / thermogenesis genes. |
| `scripts/51_rerun_with_recovered.py` | Recount everything after recovering the unannotated genes. |
| `scripts/52_recover_rna_layer.py` | Carry the gene recovery over to the RNA-seq layer. |
| `scripts/53_three_checks.py` | The three checks the recovered result needs before it can be stated. |
| `scripts/54_bird_core_sets.py` | Step one of the story: the two gene sets and their seven-bird cores. |
| `scripts/55_working_sets_fits.py` | The working sets and their numbers. |
| `scripts/56_recovery_round2.py` | Second recovery round, driven by thermogenesis. |
| `scripts/57_working_sets_round2.py` | Working sets and fits after BOTH recovery rounds. |
| `scripts/58_full_reannotation.sh` | Full protein-level re-annotation of the three Augustus-annotated assemblies. |
| `scripts/59_emu_liftover.py` | Lift the emu CAGE peaks from ZJU1.0 onto bDroNov1 and use the real annotation. |
| `scripts/60_unbiased_mass_scan.py` | Which genes track body mass, without asking a pathway map first. |
| `scripts/61_metabolic_panel.py` | Does anything else scale with body mass, or only the respiratory chain? |
| `scripts/62_build_library_matrices.py` | Gene x library matrices for WGCNA, plus the widened gene universe. |
| `scripts/63_variance_partition.py` | How much of the variance is species, how much tissue, how much left over. |
| `scripts/64_wgcna_run.R` | Standard WGCNA on the bird CAGE matrices built by script 62. |
| `scripts/64b_consensus_fix.R` | Consensus modules over the seven species as seven sets, inside fibroblasts. |
| `scripts/65_modules_vs_mass.py` | Modules against body mass, and where the 55 OXPHOS genes landed. |
| `scripts/66_module_function.py` | What the WGCNA modules are, functionally. |
| `scripts/67_modules_within_species.py` | Refining the modules: co-expression within species, not between them. |
| `scripts/68_rebuild_matrices_asymmetric.py` | Gene x library matrices rebuilt on a symmetric footing. |
| `scripts/69_fetch_native_assemblies.sh` | Assemblies and proteomes for the four species that never got a re-annotation. |
| `scripts/70_common_reannotation.sh` | One proteome projected onto all seven genomes, same query everywhere. |
| `scripts/71_wgcna_v2.R` | WGCNA on the rebuilt matrices of script 68 (common_only, window 5000 up / 1000 down). |
| `scripts/71b_consensus_v2.R` | Consensus on the v2 fibroblast matrix. Same fix as 64b: multiData() instead of |
| `scripts/72_modules_v2.py` | Modules vs mass on the rebuilt panel, plus the centred-within-species network. |
| `scripts/73_phylo_oxphos.py` | Phylogenetic calibration of the OXPHOS working-set decline, n = 7. |
| `scripts/74_nineteen_via_cage.py` | The 19 genes that carried the old ±1000 decline, measured via CAGE. |
| `scripts/75_honest_gene_mass_scan.py` | Per-gene correlation with adult mass under the homologous-CAGE rule. |
| `scripts/76_set_scan_honest.py` | Which gene SETS, as a sum, track body mass on the honest CAGE matrix. |
| `scripts/77_kegg_zeros.py` | KEGG set scan in Ruslan's metric: official list, missing gene×species = 0. |
| `scripts/78_holding_cores.py` | Holding cores inside kept KEGG maps and inside each WGCNA module. |
| `scripts/79_story_ruslan_scatters.py` | Ruslan-style log–log scatters for the five locked story sets. |
| `scripts/80_story_epi_and_traits.py` | Story sets in epiblast, and the same sums vs Ruslan's other three traits. |
| `scripts/81_story_cre_motifs.py` | Story-set CRE: proximal vs distal homologous CAGE, and motif enrichment. |
| `scripts/82_create_fib_atlas.sh` | Fibroblast CREate for the seven birds. |
| `scripts/82_monalisa_binned.R` | monaLisa::calcBinnedMotifEnrR on one species. |
| `scripts/82_motif_method_compare.py` | Three lightweight motif tests next to the frozen JASPAR FIMO. |
| `scripts/83_ruslan_bidir_enhancers.py` | Transcribed enhancers the way Ruslan called them. |
| `scripts/84_create_peaks_to_story_genes.py` | Assign fibroblast CREate regions to the story gene sets. |
| `scripts/85_ela_distal_fimo.py` | Motifs on distal ELA of the story sets, the layer HEY2 came from. |
| `scripts/86_ela_distance_to_pla.py` | Distance of proximal ELA to the promoter PLA of the same gene. |
| `scripts/87_cross_tissue_size_set.py` | Find a size-linked gene set on honest promoters, confirmed across tissues. |
| `scripts/88_mammal_honest_matrices.py` | Honest gene x species matrices for three mammal fibroblast panels. |
| `scripts/88_size_group_search.py` | Gene groups vs body mass, beyond KEGG. |
| `scripts/89_mammal_set_scan.py` | A priori gene-set sums vs log10 mass on the three honest mammal matrices. |
| `scripts/90_mammal_camera_gsea.py` | CAMERA and rank GSEA on the honest mammal matrices. |
| `scripts/91_mammal_pgls.py` | PGLS for mammal gene sets that passed the abundance-null BH. |
| `scripts/92_mammal_size_report.py` | Report for the mammal size scan. Does not edit the frozen bird nine. |
| `scripts/93_panel_background.py` | Where a set's mass slope sits against its own panel. |
| `scripts/94_monalisa_nine.R` | monaLisa::calcBinnedMotifEnrR on one species. |
| `scripts/94_nine_set_motifs.py` | Motif enrichment for the nine locked bird sets. |
| `scripts/95_nine_set_slopes.py` | Nine locked sets, fibroblast and epiblast, Ruslan-style log-log. |
| `scripts/96_mammal_pass_slopes.py` | One Ruslan-style grid of mammal sets that passed abundance-null BH. |
| `scripts/97_mammal_core_slopes.py` | Cores of the seven mammal sets, cut the way the bird slices were cut. |
| `scripts/98_epi_own_promoters.py` | Epiblast promoters with the same 2 kb rule, without the fibroblast filter. |
| `scripts/bg_full_pipeline.sh` | Redownload if missing or suspiciously small (<200MB; full gene_info is usually larger) |
| `scripts/oxphos_thread_lib.py` | Shared helpers for OXPHOS species-thread pipeline. |
| `scripts/run_all.sh` | Original analysis script; inspect configuration and read/write calls. |

Inputs and outputs retain the original research-run hierarchy. No new automatic runner was added. Machine-specific prefixes were replaced by explicit placeholders; configure these before executing an original stage.
