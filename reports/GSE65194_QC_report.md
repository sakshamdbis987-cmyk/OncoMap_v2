# GSE65194 preliminary data-loading and QC report

## Summary

**Result:** PASSED the preliminary structural and descriptive QC checks.

These are preliminary file-structure, sample-mapping, expression-summary, and
technical-replicate checks. Passing them does **not** biologically validate the
dataset, establish clinical validity, or confirm that the GEO preprocessing
description is correct.

No normalization, log transformation, sample removal, replicate averaging,
differential-expression analysis, or subtype expression comparison was performed.
Expression values were analyzed as stored in the GEO series matrix.

## Software

- Python: 3.15.2
- pandas: 3.1.0
- NumPy: 2.5.4
- matplotlib: 3.11.2
- SciPy: not used

## Matrix and manifest validation

| Item | Result |
|---|---:|
| Matrix dimensions (probes × arrays) | 54673 × 178 |
| Unique probes | 54673 |
| Unique matrix GSMs | 178 |
| Manifest records | 178 |
| Unique manifest GSMs | 178 |
| Matrix/manifest exact GSM match | PASS |
| Probe index | ID_REF |
| Matrix minimum | -1.222691 |
| Matrix maximum | 19.999146 |
| Source SHA-256 | f92d2109fe0d64ed056d9673e8733b2af25173230f3632e370201b6a4054caba |

The manifest columns found were: `gsm_id, sample_title, tumor_id, replicate_label, sample_group, array_batch, hybridation_batch, source_name, cell_line, is_cell_line, platform, clinical_record_gsm_id, patient_clinical_record_matched, direct_gsm_clinical_match`.
The `sample_group` values were preserved as stored. `tumor_id` was used as the
patient identifier; healthy and cell-line labels were kept distinct.

## Sample composition

- TNBC: 55
- Her2: 39
- Luminal A: 29
- Luminal B: 30
- Healthy: 11
- CellLine: 14

Expected-count check:

- All sample-group counts match the expected counts.

- Tumour arrays: 153
- Unique tumour patients: 130 (expected 130)
- Patients with one array: 107 (expected 107)
- Patients with two arrays: 23 (expected 23)
- Patients with more than two arrays: 0

The matrix/manifest sample-order match was validated by GSM accession. The
manifest was aligned to the expression matrix column order without changing
either input file.

## Expression distribution findings

| Item | Result |
|---|---:|
| Probe rows | 54673 |
| Unique probe identifiers | 54673 |
| Array columns | 178 |
| Unique GSM identifiers | 178 |
| Missing expression values | 0 |
| Non-numeric expression values | 0 |
| Minimum expression value | -1.222691 |
| Maximum expression value | 19.999146 |
| Median of array medians | 3.217397 |
| Lowest array median | 3.142912 |
| Highest array median | 3.705825 |
| Median-array IQR | 0.044323 |

Missing-value check: **0** missing values were found
in the expression matrix. Non-numeric value check: **0**
non-numeric values were found.

The following arrays are candidates for visual/follow-up inspection because their sample median is outside the 1.5×IQR screening range across all arrays: GSM1588986 (TNBC, TUM013_repA, median 3.3118); GSM1589121 (Luminal B, TUM123, median 3.3236); GSM1589129 (CellLine, 184B5, median 3.3082); GSM1589130 (Healthy, H3, median 3.6046); GSM1589132 (Healthy, H13, median 3.7058); GSM1589136 (Healthy, H10, median 3.3431); GSM1589139 (Healthy, H18, median 3.6158); GSM1589142 (Healthy, H12, median 3.4636); GSM1589144 (Healthy, H11, median 3.4032); GSM1589145 (Healthy, H15, median 3.4132); GSM1589147 (CellLine, BT-549, median 3.3301); GSM1589148 (Healthy, H4, median 3.5425); GSM1589150 (Healthy, H14, median 3.3543); GSM1589151 (Healthy, H5, median 3.3533); GSM1589152 (CellLine, MDA-MB-157, median 3.3132). This is a review flag only, not an automatic outlier classification; no samples were removed.

The distribution and sample-median figures use distinct colors and a legend for
the preserved manifest sample groups. They are descriptive plots and do not
perform subtype comparisons.

## Technical replicate findings

Exactly 23 tumour replicate pairs were evaluated across all
54673 available probes per pair. Correlation is Pearson's
correlation; difference metrics use array 1 minus array 2 for the signed
difference. No replicate values were averaged.

| tumor_id | gsm_id_array_1 | replicate_label_array_1 | gsm_id_array_2 | replicate_label_array_2 | probe_count | pearson_correlation | median_absolute_difference | mean_signed_difference_array_1_minus_array_2 | rmse |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TUM001 | GSM1588970 | repA | GSM1588971 | repB | 54673 | 0.993865 | 0.071268 | -0.009532 | 0.328061 |
| TUM002 | GSM1588972 | repA | GSM1588973 | repB | 54673 | 0.996806 | 0.070404 | 0.011000 | 0.235214 |
| TUM006 | GSM1588977 | repA | GSM1588978 | repB | 54673 | 0.996350 | 0.066041 | -0.001595 | 0.252300 |
| TUM010 | GSM1588982 | repA | GSM1588983 | repB | 54673 | 0.994025 | 0.078458 | -0.012053 | 0.318345 |
| TUM013 | GSM1588986 | repA | GSM1588987 | repB | 54673 | 0.993960 | 0.094571 | -0.003379 | 0.321757 |
| TUM015 | GSM1588989 | repA | GSM1588990 | repB | 54673 | 0.995748 | 0.069750 | 0.033427 | 0.270774 |
| TUM016 | GSM1588991 | repA | GSM1588992 | repB | 54673 | 0.996801 | 0.058803 | -0.004193 | 0.231353 |
| TUM017 | GSM1588993 | repA | GSM1588994 | repB | 54673 | 0.996711 | 0.059755 | -0.024424 | 0.241683 |
| TUM018 | GSM1588995 | repA | GSM1588996 | repB | 54673 | 0.995430 | 0.064186 | 0.004627 | 0.282360 |
| TUM020 | GSM1588998 | repA | GSM1588999 | repB | 54673 | 0.996953 | 0.066498 | 0.008970 | 0.233415 |
| TUM021 | GSM1589000 | repA | GSM1589001 | repB | 54673 | 0.996546 | 0.059648 | -0.008881 | 0.244197 |
| TUM026 | GSM1589006 | repA | GSM1589007 | repB | 54673 | 0.996407 | 0.060661 | 0.000956 | 0.243836 |
| TUM028 | GSM1589009 | repA | GSM1589010 | repB | 54673 | 0.995959 | 0.088768 | 0.009195 | 0.267514 |
| TUM032 | GSM1589014 | repA | GSM1589015 | repB | 54673 | 0.994833 | 0.067782 | -0.016355 | 0.297272 |
| TUM035 | GSM1589018 | repA | GSM1589019 | repB | 54673 | 0.996726 | 0.061548 | 0.019391 | 0.238030 |
| TUM038 | GSM1589022 | repA | GSM1589023 | repB | 54673 | 0.996111 | 0.057786 | 0.003491 | 0.259651 |
| TUM043 | GSM1589028 | repA | GSM1589029 | repB | 54673 | 0.995537 | 0.063648 | 0.007971 | 0.275086 |
| TUM045 | GSM1589031 | repA | GSM1589032 | repB | 54673 | 0.996351 | 0.064504 | -0.031137 | 0.255930 |
| TUM048 | GSM1589035 | repA | GSM1589036 | repB | 54673 | 0.996679 | 0.064171 | -0.015344 | 0.238370 |
| TUM051 | GSM1589039 | repA | GSM1589040 | repB | 54673 | 0.996226 | 0.067779 | 0.025459 | 0.258833 |
| TUM058 | GSM1589047 | repA | GSM1589048 | repB | 54673 | 0.996578 | 0.059839 | -0.020386 | 0.244022 |
| TUM064 | GSM1589054 | repA | GSM1589055 | repB | 54673 | 0.994059 | 0.076476 | -0.005425 | 0.318183 |
| TUM067 | GSM1589058 | repA | GSM1589059 | repB | 54673 | 0.995728 | 0.063606 | 0.002839 | 0.271062 |

High technical-replicate correlation does not by itself prove that the GEO
processing description's statement about replicate averaging is correct. The
matrix still contains separate GSM columns for these arrays.

## Warnings and unresolved questions

- The GEO record describes GC-RMA processing and mentions averaging technical
  replicates, while this matrix retains separate replicate GSM columns. These
  QC calculations do not resolve that discrepancy.
- The exact Affymetrix CDF used for processing remains unresolved in the GEO
  description, which mentions either a standard AffyCDF or a Brainarray custom
  definition.
- Sample-median screening is a heuristic only. Investigate any candidates with
  relevant technical and biological context before making decisions.
- No differential-expression tests, log fold changes, multiple-testing
  corrections, sample exclusions, or biological conclusions were produced.

## Input integrity

The compressed matrix SHA-256 before and after this run was identical:
`f92d2109fe0d64ed056d9673e8733b2af25173230f3632e370201b6a4054caba`. The GEO matrix and sample manifest were read-only inputs.
All generated files were written under the project `reports` directory.
