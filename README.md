# Body-mass transcriptomics: analysis code

Existing analysis scripts for **Adult body size is associated with distinct transcriptional signatures in avian epiblast and primary fibroblasts**.

Ilya Rubinshteyn, Ruslan Deviatiiarov, Oleg Gusev and Guojun Sheng.

This repository organizes the code used during the project. It does not introduce a replacement analysis, a new statistical implementation or a new reproduction framework. Original script names, numbering and analytical logic are retained. Failed development variants and corpus-management scripts are excluded.

## Start here

1. Read the [analysis map](docs/analysis-map.md) to follow the actual sequence of work.
2. Open the README inside the relevant analysis branch and the [entry-point guide](docs/starting-points.md).
3. Check the [inputs and configuration](docs/inputs-and-configuration.md) and restore the original input versions before running that branch.
4. Consult `provenance/source_manifest.json` for the original source checksum and the repository-copy checksum.

| Directory | Analysis branch |
|---|---|
| `analyses/avian/scripts/` | Original bird CAGE, reassignment, gene-group, network and promoter scripts, including exploratory and diagnostic steps |
| `analyses/mammal_wgcna/code/` | Original mammalian expression preparation, WGCNA, mass tests, model/leave-one-species checks and figures |
| `analyses/mammal_maps/code/` | Original fixed gene-set sums, promoter analyses and phylogenetic comparisons |
| `analyses/mammal_scenic/code/` | Original pySCENIC preparation, network inference, rankings, pruning, AUCell and seed sensitivity |
| `figures/` | Existing panel-drawing and target-heatmap scripts |
| `data/frozen/` | Selected existing numerical inputs, memberships and frozen scores used for the poster |
| `environments/` | Saved historical environment information, where available |
| `docs/` | Order of work, inputs, configuration and result boundaries |
| `provenance/` | File identities and explicit changes made during organization |

The scripts within each branch remain together because they import shared helpers and numbered modules. Moving each numbered script into its own stage folder would break these dependencies. The stage READMEs provide the ordering while preserving the original executable context.

## What was changed during repository assembly

Only organization and documentation were added. Machine-specific path prefixes in source files were replaced by explicit `/path/to/…` placeholders. Source and repository checksums make these edits visible. Configure these paths for your own copy of the original inputs; they are not working defaults. No algorithm, threshold, seed, model or correction family was changed. No full analyses were rerun during assembly.

## Scientific scope

The project follows gene-level measurements through fixed or refined gene groups to associations with adult mass and candidate shared regulation. Birds and mammals are complementary parts. Pathway selection and its later tests used the same data; the displayed nine-group correction does not account for the preceding selection. OXPHOS conclusions depend on the specified filter. Splicing significance depends on the model, and pySCENIC significance varied across runs. Motifs and inferred targets are candidate evidence, not established causal size regulators.

See [result boundaries](docs/result-boundaries.md) for details. A source checksum establishes which code was preserved; it does not by itself establish independent biological replication or availability of every upstream dependency.

## Sources and reuse

Avian fibroblasts: **GSE213253**. Mammalian fibroblasts: **Ma et al. (2016), eLife 5:e19130**, https://doi.org/10.7554/eLife.19130. Avian epiblast CAGE was supplied by the study authors; the saved project metadata do not provide a public accession. Raw sequencing, genomes, large motif/ranking databases, model environments and intermediate poster versions are not redistributed here.

A reuse license has not yet been selected by the project authors. Third-party inputs retain their own terms.
