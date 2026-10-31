# GSE65194 targeted differential-expression report

## Scope and caution

This is a **targeted analysis of 38 validated GPL570 probe sets** representing
eight genes. It is not genome-wide differential-expression analysis. Results
are probe-level; multiple probes for one gene are not averaged into an invented
gene-expression value.

The analysis is exploratory and does not establish causation, diagnostic
performance, or independent validation. A statistically significant p-value
does not by itself validate a biological conclusion.

No input matrix, manifest, compressed GEO source, or prior report was modified.
No samples were removed. The patient-level matrix already represents the
documented aggregation of replicate pairs; this script did not re-average
replicates.

## Inputs and validation

- Patient-level matrix: 54,673 probes × 141 observations.
- Tumour patients: 130; healthy controls: 11.
- Tumour subtype counts: TNBC 41, Her2 30, Luminal A 29, Luminal B 30.
- All 141 sample identifiers matched exactly between matrix columns and
  manifest rows, with no duplicates.
- All 38 target probes were present.
- No missing, non-numeric, or infinite expression values were found.
- SciPy version: 1.18.1.

## Methods

For each probe, two-sided Welch independent-samples t-tests were run using
`scipy.stats.ttest_ind(equal_var=False)`. The expression values are the
processed log2-scale values as stored in the derived patient-level matrix;
there was no additional transformation or normalization.

### Tumour versus healthy

All 130 tumour patients were compared with all 11 healthy controls. The
reported mean log2 difference is tumour mean minus healthy-control mean.
Positive values indicate higher mean expression in tumour; negative values
indicate lower mean expression in tumour. Benjamini-Hochberg (BH) correction
was applied across the complete family of **38 probe tests**. Significance is
defined as BH-adjusted p < 0.05.

### Pairwise subtype comparisons

All six two-sided contrasts were tested: TNBC−Her2, TNBC−Luminal A,
TNBC−Luminal B, Her2−Luminal A, Her2−Luminal B, and Luminal A−Luminal B.
The reported difference is first subtype mean minus second subtype mean.
Benjamini-Hochberg correction was applied once across the complete family of
**228 tests (6 contrasts × 38 probes)**, not separately within each contrast.
Significance is defined as adjusted p < 0.05.

The recorded subtypes are analyzed as listed; this analysis does not infer
ER-positive status from Luminal A or Luminal B labels.

## Tumour versus healthy results

- Probes passing BH-adjusted p < 0.05: **27 of 38**.
- Positive tumour-minus-healthy differences: **19**.
- Negative tumour-minus-healthy differences: **19**.

Significant probe results:

| probe_id | gene_symbol | n_tumour | n_healthy | mean_log2_difference_tumour_minus_healthy | p_value_two_sided | bh_adjusted_p_value |
| --- | --- | --- | --- | --- | --- | --- |
| 211851_x_at | BRCA1 | 130 | 11 | 0.156438 | 0.00231525 | 0.00439897 |
| 208368_s_at | BRCA2 | 130 | 11 | 0.752426 | 0.000510473 | 0.00121237 |
| 214727_at | BRCA2 | 130 | 11 | 0.205522 | 0.00738645 | 0.0122037 |
| 201746_at | TP53 | 130 | 11 | -0.825305 | 0.000127065 | 0.000536498 |
| 211300_s_at | TP53 | 130 | 11 | 1.22475 | 1.46979e-08 | 1.86173e-07 |
| 210930_s_at | ERBB2 | 130 | 11 | 2.09064 | 1.38919e-10 | 4.3154e-09 |
| 216836_s_at | ERBB2 | 130 | 11 | 1.88627 | 2.68098e-07 | 2.03754e-06 |
| 234354_x_at | ERBB2 | 130 | 11 | 1.21077 | 2.90309e-08 | 2.75794e-07 |
| 204369_at | PIK3CA | 130 | 11 | 2.93926 | 0.000250664 | 0.000924865 |
| 231854_at | PIK3CA | 130 | 11 | -0.392799 | 0.012845 | 0.0196281 |
| 235980_at | PIK3CA | 130 | 11 | -1.3486 | 7.38827e-05 | 0.000401078 |
| 211233_x_at | ESR1 | 130 | 11 | 0.144858 | 0.000702759 | 0.0014836 |
| 211234_x_at | ESR1 | 130 | 11 | 0.152416 | 0.000267724 | 0.000924865 |
| 211235_s_at | ESR1 | 130 | 11 | 0.344181 | 0.000491962 | 0.00121237 |
| 211627_x_at | ESR1 | 130 | 11 | -0.879169 | 0.000326471 | 0.00103382 |
| 215552_s_at | ESR1 | 130 | 11 | 0.674832 | 0.000966891 | 0.00193378 |
| 204054_at | PTEN | 130 | 11 | -0.306884 | 0.00398781 | 0.00688804 |
| 211711_s_at | PTEN | 130 | 11 | -0.95064 | 0.000476217 | 0.00121237 |
| 217492_s_at | PTEN | 130 | 11 | 0.711807 | 0.0298744 | 0.0420454 |
| 225363_at | PTEN | 130 | 11 | -1.00955 | 6.28491e-07 | 3.98045e-06 |
| 228006_at | PTEN | 130 | 11 | -1.59645 | 2.27126e-10 | 4.3154e-09 |
| 242622_x_at | PTEN | 130 | 11 | -0.615822 | 0.0024858 | 0.00449811 |
| 201983_s_at | EGFR | 130 | 11 | -1.89773 | 0.000398197 | 0.00116396 |
| 201984_s_at | EGFR | 130 | 11 | -0.993656 | 0.0129132 | 0.0196281 |
| 211550_at | EGFR | 130 | 11 | -0.151177 | 0.000668825 | 0.0014836 |
| 211551_at | EGFR | 130 | 11 | -0.0262287 | 0.0255528 | 0.0373464 |
| 224999_at | EGFR | 130 | 11 | -1.95798 | 0.000112155 | 0.000532736 |

Complete probe-level results are in
`GSE65194_tumor_vs_healthy_welch.tsv`.

## Pairwise subtype results

| Contrast (first minus second) | First n | Second n | BH-significant probes in contrast |
|---|---:|---:|---:|
| TNBC - Her2 | 41 | 30 | 12 |
| TNBC - Luminal A | 41 | 29 | 19 |
| TNBC - Luminal B | 41 | 30 | 22 |
| Her2 - Luminal A | 30 | 29 | 17 |
| Her2 - Luminal B | 30 | 30 | 15 |
| Luminal A - Luminal B | 29 | 30 | 10 |

There are **95 of 228** subtype-probe tests with
BH-adjusted p < 0.05 across the full family. Complete results, including
non-significant probes, are in
`GSE65194_subtype_pairwise_welch.tsv`.

## Probe directions by gene

The following direction summary is based on tumour-minus-healthy mean
differences. The sign counts refer to separate probes, not a combined gene
measurement.

| gene_symbol | probe_count | positive_probe_count_tumour_vs_healthy | negative_probe_count_tumour_vs_healthy | tumour_vs_healthy_direction_consistency | tumour_vs_healthy_bh_significant_probes |
| --- | --- | --- | --- | --- | --- |
| BRCA1 | 2 | 2 | 0 | all_positive | 211851_x_at |
| BRCA2 | 2 | 2 | 0 | all_positive | 208368_s_at;214727_at |
| TP53 | 2 | 1 | 1 | mixed | 201746_at;211300_s_at |
| ERBB2 | 3 | 3 | 0 | all_positive | 210930_s_at;216836_s_at;234354_x_at |
| PIK3CA | 3 | 1 | 2 | mixed | 204369_at;231854_at;235980_at |
| ESR1 | 9 | 5 | 4 | mixed | 211233_x_at;211234_x_at;211235_s_at;211627_x_at;215552_s_at |
| PTEN | 8 | 3 | 5 | mixed | 204054_at;211711_s_at;217492_s_at;225363_at;228006_at;242622_x_at |
| EGFR | 9 | 2 | 7 | mixed | 201983_s_at;201984_s_at;211550_at;211551_at;224999_at |

`217492_s_at` retains the annotation `PTENP1 /// PTEN`; because the annotation
is ambiguous, it should not be attributed uniquely to PTEN. All individual
probes remain separately visible in both result tables.

## ESR1 and ERBB2 across the recorded subtypes

For ESR1 and ERBB2, the gene summary records each probe's mean expression
within TNBC, Her2, Luminal A, and Luminal B separately. These values are
probe-specific; no gene-level average is calculated. Adjusted-significant
pairwise subtype probes are listed below.

| gene_symbol | probe_count | subtype_pairwise_bh_significant_probes | ESR1_ERBB2_probe_means_by_subtype |
| --- | --- | --- | --- |
| ERBB2 | 3 | TNBC_vs_Her2:210930_s_at;TNBC_vs_Her2:216836_s_at;TNBC_vs_Her2:234354_x_at;TNBC_vs_Luminal A:216836_s_at;TNBC_vs_Luminal B:210930_s_at;TNBC_vs_Luminal B:216836_s_at;TNBC_vs_Luminal B:234354_x_at;Her2_vs_Luminal A:210930_s_at;Her2_vs_Luminal A:216836_s_at;Her2_vs_Luminal A:234354_x_at;Her2_vs_Luminal B:210930_s_at;Her2_vs_Luminal B:216836_s_at;Her2_vs_Luminal B:234354_x_at;Luminal A_vs_Luminal B:210930_s_at;Luminal A_vs_Luminal B:216836_s_at;Luminal A_vs_Luminal B:234354_x_at | 210930_s_at[TNBC=2.596544, Her2=8.835914, Luminal A=2.923014, Luminal B=5.990229]; 216836_s_at[TNBC=9.292487, Her2=14.482989, Luminal A=10.335252, Luminal B=13.560604]; 234354_x_at[TNBC=2.464550, Her2=6.004396, Luminal A=2.463770, Luminal B=4.352527] |
| ESR1 | 9 | TNBC_vs_Luminal A:205225_at;TNBC_vs_Luminal A:211233_x_at;TNBC_vs_Luminal A:211234_x_at;TNBC_vs_Luminal A:211235_s_at;TNBC_vs_Luminal A:215552_s_at;TNBC_vs_Luminal A:217190_x_at;TNBC_vs_Luminal B:205225_at;TNBC_vs_Luminal B:211234_x_at;TNBC_vs_Luminal B:211235_s_at;TNBC_vs_Luminal B:215551_at;TNBC_vs_Luminal B:215552_s_at;TNBC_vs_Luminal B:217190_x_at;Her2_vs_Luminal A:205225_at;Her2_vs_Luminal A:211233_x_at;Her2_vs_Luminal A:211234_x_at;Her2_vs_Luminal A:211235_s_at;Her2_vs_Luminal A:215552_s_at;Her2_vs_Luminal B:205225_at;Her2_vs_Luminal B:211234_x_at;Her2_vs_Luminal B:211235_s_at;Her2_vs_Luminal B:215551_at;Her2_vs_Luminal B:215552_s_at | 205225_at[TNBC=3.251582, Her2=3.598402, Luminal A=11.462228, Luminal B=11.207202]; 211233_x_at[TNBC=1.879967, Her2=1.880961, Luminal A=2.333760, Luminal B=2.095701]; 211234_x_at[TNBC=2.041318, Her2=2.043773, Luminal A=2.510899, Luminal B=2.334434]; 211235_s_at[TNBC=2.374801, Her2=2.378560, Luminal A=3.563506, Luminal B=3.119239]; 211627_x_at[TNBC=3.479000, Her2=3.480437, Luminal A=3.534301, Luminal B=3.498226]; 215551_at[TNBC=2.354307, Her2=2.348282, Luminal A=2.592509, Luminal B=2.675824]; 215552_s_at[TNBC=2.496554, Her2=2.512877, Luminal A=4.741416, Luminal B=4.402672]; 217163_at[TNBC=2.041454, Her2=2.042292, Luminal A=2.036555, Luminal B=2.039051]; 217190_x_at[TNBC=2.298791, Her2=2.337544, Luminal A=2.544221, Luminal B=2.405082] |

Use the full subtype result table for each individual probe's six contrasts,
mean difference, and family-wide adjusted p-value. A subtype mean pattern is
descriptive and does not by itself establish receptor status or biological
mechanism.

## Limitations

- This tests only the 38 selected probes; it is not genome-wide screening.
- Multiple probes may measure different transcript regions or have ambiguous
  annotations. Probe-level direction disagreement is retained rather than
  concealed by averaging.
- The healthy-control sample size is 11, compared with 130 tumour patients;
  group sizes are reported for transparency.
- The upstream audit notes that the exact CDF and explicit log2 transformation
  history are not fully resolved. The analysis uses the processed values as
  supplied and does not renormalize them.
- GEO's processing text mentions averaging technical replicates even though
  separate GSM columns are present in its downloadable matrix. The patient-
  level matrix is the project's explicit documented aggregation; this targeted
  analysis does not resolve GEO's original processing provenance.
- Luminal labels were not converted to ER-positive labels. Clinical receptor
  status is not inferred by this subtype analysis.
- Statistical significance is not proof of biological validation, causality,
  or diagnostic performance.

## Outputs

- `GSE65194_tumor_vs_healthy_welch.tsv`
- `GSE65194_subtype_pairwise_welch.tsv`
- `GSE65194_gene_direction_summary.tsv`
- `GSE65194_differential_expression_report.md`
