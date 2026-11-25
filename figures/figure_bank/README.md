# figures/figure_bank

Existing scripts; original names and analytical logic retained. See the [analysis map](../../docs/analysis-map.md) and [input/configuration guide](../../docs/inputs-and-configuration.md).

| File | Source purpose / first documentation line |
|---|---|
| `code/draw_scientific_panels.py` | Render fixed, existing P-KLB measurements on the research host. |
| `code/refine_motif_panels.py` | Re-render only the two motif tables exported by draw_scientific_panels.py. |

Inputs and outputs retain the original research-run hierarchy. No new automatic runner was added. Machine-specific prefixes were replaced by explicit placeholders; configure these before executing an original stage.
