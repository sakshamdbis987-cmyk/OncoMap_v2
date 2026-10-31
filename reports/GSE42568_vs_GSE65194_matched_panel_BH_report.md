# Matched-panel BH significance comparison and GEO2R assignment sheets

## Matched-panel significance comparison

The comparison matches **exact probe IDs only** across the original GSE42568 GEO2R table and the GSE65194 tumour-versus-healthy Welch table. All 38 shared probes were retained. No effect estimate, contrast definition, or direction comparison was recalculated or changed; the existing values are carried forward unchanged.

Benjamini?Hochberg (BH) correction was recomputed separately for exactly these 38 probes in each cohort from the raw two-sided p-values (`P.Value` for GSE42568 and `p_value_two_sided` for GSE65194). The existing adjusted p-values are preserved in separate `*_adjusted_p_value_original` columns, alongside the new `*_adjusted_p_value_matched_38_BH` values. Significance is defined as adjusted p < 0.05.

| Cohort | Significant using original adjusted p-values | Significant using matched-panel BH | Significant both | No longer significant | Newly significant |
|---|---:|---:|---:|---:|---:|
| GSE42568 | 18 | 18 | 18 | 0 | 0 |
| GSE65194 | 27 | 27 | 27 | 0 | 0 |

Because GSE42568's original adjusted p-values cover its genome-wide GEO2R testing family, while GSE65194's original values cover its targeted 38-probe family, the original significance columns use different correction families. The matched-panel values provide a like-for-like panel-size correction, but do not harmonize cohort design, preprocessing, or statistical method.

## GEO2R sample assignment sheets

The assignments were generated from the existing sample manifest. Tumour samples were grouped by `tumor_id`; paired patients contribute one array per sheet (repA in the repA sheet and repB in the repB sheet), while each singleton patient contributes its only array to both sheets. Healthy controls are included individually. All CellLine arrays are excluded. The GEO2R group is `Tumour` or `Healthy`; the `original_subtype` field retains the manifest's `sample_group` label (including `Healthy` for controls).

| Sheet | Tumour GSMs | Unique tumour patients | Healthy GSMs | Cell-line GSMs | Total GSMs | Unique GSMs |
|---|---:|---:|---:|---:|---:|---:|
| repA | 130 | 130 | 11 | 0 | 141 | 141 |
| repB | 130 | 130 | 11 | 0 | 141 | 141 |

For both sheets, all 130 tumour patients are represented exactly once, all 11 healthy GSMs are included, and none of the 14 cell-line GSMs is included. Each selected GSM is present in the manifest. The `selection_reason` values identify `repA`, `repB`, `singleton`, or `healthy control` as applicable.

These are assignment files only. GEO2R was not run, and no existing analysis output or source file was modified.
