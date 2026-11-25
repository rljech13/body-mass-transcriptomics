# analyses/mammal_wgcna

Existing scripts; original names and analytical logic retained. See the [analysis map](../../docs/analysis-map.md) and [input/configuration guide](../../docs/inputs-and-configuration.md).

| File | Source purpose / first documentation line |
|---|---|
| `code/01_prepare.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/02_wgcna.R` | Original analysis script; inspect configuration and read/write calls. |
| `code/03_native_pwm.py` | Exact initial TF list, kept explicit to make compatibility with plan schema auditable. |
| `code/04_statistics.py` | Learn only in discovery panels. q correction includes ALL candidate edges in all six panels. |
| `code/05_human_pwm.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/06_wgcna_LOO.R` | Original analysis script; inspect configuration and read/write calls. |
| `code/07_pwm_statistics.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/08_additional.py` | An independent whitening check against base R lm is prepared without importing this fitter into R. |
| `code/09_independent_GLS.R` | Original analysis script; inspect configuration and read/write calls. |
| `code/10_prepare_full.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/11_wgcna_full.R` | Original analysis script; inspect configuration and read/write calls. |
| `code/12_statistics_full.py` | Learn only in discovery panels. q correction includes ALL candidate edges in all six panels. |
| `code/13_own5_full_original.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/14_full_pipeline.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/15_wgcna_full_LOO.R` | Original analysis script; inspect configuration and read/write calls. |
| `code/16_tree_sensitivity.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/17_human_own5_union.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/18_human_pipeline.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/19_Hallmark_and_gene_jackknife.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/20_QA.py` | Exhaustive tiny conditional null, independent of sequence enrichment data. |
| `code/21_figures.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/common.py` | Original analysis script; inspect configuration and read/write calls. |
| `code/module_scores.R` | Original analysis script; inspect configuration and read/write calls. |
| `code/native_bridge_lib.py` | Original analysis script; inspect configuration and read/write calls. |

Inputs and outputs retain the original research-run hierarchy. No new automatic runner was added. Machine-specific prefixes were replaced by explicit placeholders; configure these before executing an original stage.
