# analyses/mammal_scenic

Existing scripts; original names and analytical logic retained. See the [analysis map](../../docs/analysis-map.md) and [input/configuration guide](../../docs/inputs-and-configuration.md).

| File | Source purpose / first documentation line |
|---|---|
| `code/01_prepare.py` | Normalize unambiguous historic symbols to HGNC accepted names once, before matching. |
| `code/02_grn.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/03_build_db.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/04_prune.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/04_prune_existing.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/04_run_existing.sh` | Original analysis script; inspect configuration and read/write calls. |
| `code/05_analyze.py` | All statistical families specified before motif pruning, separate from seed diagnostics. |
| `code/06_verify.R` | Original analysis script; inspect configuration and read/write calls. |
| `code/07_draw.R` | Display the two leading primary candidates per module plus the preexisting E2F2 candidate. |
| `code/08_stability.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/09_export.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/run_analysis.sh` | Original analysis script; inspect configuration and read/write calls. |
| `code/run_databases.sh` | Original analysis script; inspect configuration and read/write calls. |
| `code/run_networks.sh` | Original analysis script; inspect configuration and read/write calls. |

Inputs and outputs retain the original research-run hierarchy. No new automatic runner was added. Machine-specific prefixes were replaced by explicit placeholders; configure these before executing an original stage.
