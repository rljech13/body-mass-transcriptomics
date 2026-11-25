# Reproduction: scope, commands and remaining gaps

## Three levels

| Level | Starting point | What the repository provides | Status |
|---|---|---|---|
| Saved-result replay | Checked-in measurements, memberships, eigengenes and fit tables | Portable plotting and the original descriptive module-statistics procedure | Bounded replay; see the execution report |
| Numerical component verification | Fixed bird points, saved sample expression and original fit functions | 18 bird associations through two existing implementations; 123 module tests; 300 target-heatmap cells and primary target rankings | Numerical checks, not independent biological replication |
| Complete study reconstruction | Raw CAGE/RNA data, annotation, orthology, reference sequences and ranking databases | Existing analysis procedures, selected historical configurations, source hashes and input map | Not independently reproduced end to end; external inputs and historical environment gaps remain |

`python -m bodymass --verify-data` checks packaged bytes against `provenance/data_manifest.json`. Hash agreement establishes integrity, not scientific correctness. The source checkout is the supported layout; an isolated binary wheel is not provided.

## Python replay without private data roots

Use Python 3.11 and `environments/replay-requirements.txt` in a new virtual environment. These pins define the bounded replay environment. They are not advertised as the environment for CAGE processing or full pySCENIC.

```bash
python -m bodymass --verify-data
python -m unittest discover -s tests -v
python -m bodymass figures all-maps --out results/all_maps
python -m bodymass figures refine-motifs --out results/motifs
```

The first figure command renders the existing three whole-map figures from `data/snapshots/all_maps/`. It reads the complete saved map-test table, saved plotted species scores, leave-one-species diagnostics and descriptive summary. It does not rebuild gene expression, select new hits or recompute the original catalogue-wide test family. The motif command renders both saved motif families separately. Output directories are created; checked-in inputs are not overwritten.

The bird-fit tests recompute Pearson r and P from the saved corrected group scores, through two original fit implementations. They do not validate the upstream CAGE assignments, the historical gene selection or independent significance after that selection.

## R replay without the original workstation

Required R packages: `jsonlite` for module statistics; additionally `ggplot2`, `ggrepel` for regulon plots, and `patchwork` for the target heatmap. `svglite` is optional: if absent, the existing plot is exported through the base R SVG device. This changes the graphics backend, not the data or calculations. The historical WGCNA R session is included under `environments/historical/`; the current replay execution report records the actually used R and package versions separately. A historical session report is not an installable lock.

```bash
python -m bodymass networks summary --out results/network_summary
python -m bodymass figures regulons --out results/regulons
python -m bodymass figures targets --out results/targets
```

`Rscript` must be on PATH. Alternatively supply a small configuration such as `{"rscript":"/your/Rscript"}` with `--config`. No research-workspace roots are needed for these three commands.

`networks summary` runs the original procedure on all 123 saved modules from Ma2016, GSE105160 and Wagner. It checks each ordinary-correlation P against its historical OLS value, computes the original within-panel and pooled BH families, and preserves historical phylogenetic P/q in separate columns. It does not discover networks again. The historical broader correction and the descriptive within-panel correction answer different questions.

The regulon plot checks r against the saved fit table before rendering. Its P/q come from the existing primary-run tests; rendering two selected regulons does not redefine their original 187-regulon correction family.

The target heatmap checks the top ten nonself targets against the included primary target-edge table, then verifies all 300 species-mean/z-score cells from the included 520 sample-level values. This checks the compact snapshot, not re-quantification of all RNA or a new GRN fit. The original wide-expression/workstation checks remain available with `--verify-originals` and a configured `paths.workstation`.

## Full analyses: exact external input contexts still matter

| Branch | Required inputs beyond the checkout | Existing configurations / outputs |
|---|---|---|
| Avian measurement and groups | Original CAGE/promoter tables; gene-location/orthology evidence; native genomes/annotations; measured matrices; historical frozen working-set tables | Configured avian root; output subdirectories under the original `results/oxphos_species_thread`; source mapping in `provenance/consolidation.json` |
| Mammalian WGCNA discovery | Panel expression/features/sample tables, phylogenetic covariance and original references | `config/historical/wgcna/`; original run's `inputs/`, `results/wgcna/` and reference tree contexts |
| Map and motif searches | Normalized matrices, canonical gene-set versions, coverage/orthology, native promoter/reference resources and FIMO | `config/historical/maps/`; `groups` workflow and saved all-map test/diagnostic tables |
| pySCENIC | Frozen expression and metadata, TF/ortholog catalogues, primary reference DB/annotation and sequence resources; relevant WGCNA memberships for overlaps | `config/historical/scenic/`; seeds 20261023/20261024/20261025; full seed-stability tables included |

Historical run configurations are evidence, not machine-local configurations ready to execute unchanged. Before running a full stage, place the required original files under the configured run/input layout and check the original source hashes. Do not treat an absent file as permission to replace it with another release or silently download today's version.

The epiblast accession is not specified in the selected records; an external reader cannot independently acquire those raw inputs from this repository alone. The exact historical pySCENIC environment is not present in the selected method, and the mammalian-map environment record is incomplete. Earlier empty environment objects have been replaced by actual available evidence and explicit missing-record statuses. Neither a container nor a complete historical environment lock is claimed.

The latest [execution report](../reports/replay.json) distinguishes a fresh checkout/environment replay from the earlier refactor tests and from full model training. Image-file hashes are recorded for traceability; rendering can differ with fonts/graphics libraries, so numeric agreement is the scientific check.
