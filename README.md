# Body-mass transcriptomics

Existing analyses for **Adult body size is associated with distinct transcriptional signatures in avian epiblast and primary fibroblasts**.

Ilya Rubinshteyn, Ruslan Deviatiiarov, Oleg Gusev and Guojun Sheng.

## Five scientific workflows, one entry point

```bash
python -m bodymass --list
```

| Workflow | Module | Work performed |
|---|---|---|
| `avian` | `bodymass/avian.py` | CAGE measurement, gene reassignment, library matrices and avian module summaries |
| `groups` | `bodymass/groups.py` | Gene-group refinement, nine-group tests, map-level associations and promoter analyses |
| `networks` | `bodymass/networks.py`, `R/networks.R` | Mammalian expression preparation, WGCNA, statistics and sensitivity analyses |
| `scenic` | `bodymass/scenic.py`, `R/regulation.R` | Preparation, GRNBoost2, rankings, motif pruning, AUCell and seed sensitivity |
| `figures` | `bodymass/figures.py`, `R/figures.R` | Existing scientific panels, regulon plots and target-expression preparation |

The calculations are in these modules. They do not dispatch to the old numbered Python files. Shared dependencies are bound in package-local namespaces so similarly named helpers from different analyses cannot overwrite one another. Loading the package or listing stages does not read datasets, create output folders or start an analysis.

## Configure the original inputs

```bash
cp config/example.json config/local.json
# Edit the roots in config/local.json to point to your existing data and tools.
python -m bodymass scenic grn --config config/local.json --plan --seed 20261023
```

`--plan` only describes the selected operation. Remove it to execute that operation using its original parameters and environment:

```bash
python -m bodymass avian cage --config config/local.json
python -m bodymass groups nine-tests --config config/local.json
python -m bodymass networks prepare --config config/local.json
python -m bodymass networks fit --config config/local.json
python -m bodymass networks statistics --config config/local.json
python -m bodymass scenic grn --config config/local.json --seed 20261023
python -m bodymass scenic prune --config config/local.json --seed 20261023
```

Arguments specific to a stage are passed to its original argument parser. The primary and sensitivity seeds remain 20261023, 20261024 and 20261025. Both original SCENIC pruning variants remain explicitly named; they are not silently substituted for each other.

Use the corresponding original scientific environment: CAGE tools, R/Bioconductor and pySCENIC have different dependencies. Saved environment descriptions are in `environments/`. Listing and planning use only the Python standard library. The package can also be installed with `python -m pip install --no-deps -e .` to expose the `bodymass` command.

## What the refactor changes

The existing analysis bodies were moved into coherent modules, with explicit execution scope and centralized paths. No new scientific method, statistical formula, membership filter, correction family or seed was introduced. The old script layout is preserved in Git commit `d230e56ab942305ec88fbfc04cc1ec8a7a6d67e2`, rather than duplicated in the current tree. That revision also retains earlier exploratory and superseded analyses.

`provenance/consolidation.json` maps each active operation back to its exact source file and hash. `provenance/original_sources.json` records the earlier source export and path sanitization. Existing numerical data and saved environment files are retained unchanged.

## Checks

```bash
python -m unittest discover -s tests -v
```

The checks compare the moved original bird-fit functions with all 18 frozen associations, exercise the GRN worker interface without training a model, and verify that planning does not start analyses. Refactoring checks are not independent biological replication. Full CAGE reassignment, WGCNA discovery and pySCENIC runs were not repeated during repository cleanup.

See [workflow and input guide](docs/workflows.md). Scientific limitations, source-data versions and the original configuration context remain relevant.

## Source data and reuse

Avian fibroblasts: GSE213253. Mammalian fibroblasts: Ma et al. (2016), eLife 5:e19130, https://doi.org/10.7554/eLife.19130. The avian epiblast measurements were supplied by the study authors; the saved metadata do not specify a public accession. Raw reads, genomes and large ranking databases remain external. A reuse license has not yet been selected by the authors.
