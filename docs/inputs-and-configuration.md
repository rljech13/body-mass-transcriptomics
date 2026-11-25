# Inputs and configuration

## Existing inputs retained here

`data/frozen/poster_inputs.json` contains the corrected seven-species bird points, saved18fits, fixed group memberships and ERG/FOXP4 activity points. `figure_inputs.json` contains saved splicing points, module statistics, relative avian profiles, motif-count cells and target-expression values. The other files retain their source names and source provenance.

The input tables document already selected results. They are not raw counts for all upstream analyses. In particular, the complete multiple-testing families for the map, regulon and motif screens are not reducible to the few displayed rows.

## Configure an analysis branch

The original scripts assumed a research directory with `data/`, `ref/`, `results/` and neighboring run directories. That directory is not recreated by this repository. Paths marked `/path/to/user-home`, `/path/to/storage` or `/path/to/workstation` require an explicit replacement with the corresponding original input location. Preserve the relative run hierarchy used in the selected branch and keep its helper scripts together.

Before running a stage, inspect its top-level constants, command-line arguments and read/write calls. Match the original genomes, annotation versions, orthology, sample registry, normalization, species masses, foreground/background and reference databases. Do not substitute current references and describe the result as an exact historical reproduction.

## Larger source data

- Avian fibroblasts: GEO **GSE213253**.
- Avian epiblast CAGE: study-author source; no public accession in the saved metadata.
- Mammalian skin fibroblasts: Ma et al. (2016), eLife **10.7554/eLife.19130**, supplementary count and metadata workbooks.
- Additional mammalian comparison panels used in historical scripts: source-specific inputs and admission/covariance registries identified by the respective scripts and manifests.
- Motif inference: recorded JASPAR collections, background gene/region sets, genome/annotation resources, FIMO outputs and custom cisTarget rankings.

Historical environment information is retained under `environments/` where it was recorded. There is no invented universal requirements file: CAGE tools, R/Bioconductor WGCNA/motif packages and pySCENIC used different environments. The assembly syntax checks do not mean these environments have been reinstalled or all stages rerun.
