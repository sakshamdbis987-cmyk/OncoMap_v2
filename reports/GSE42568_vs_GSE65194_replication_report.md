# GSE42568 vs GSE65194 cross-cohort replication audit

## Result source and contrast

**Usable GSE42568 result table found:** the existing probe-level GSE42568 result table used for this comparison. The source full table is not distributed in this repository.

This is the existing probe-level GEO2R/limma result table. Its actual columns are `ID`, `adj.P.Val`, `P.Value`, `t`, `B`, `logFC`, `Gene.symbol`, and `Gene.title`; it contains 54,675 rows and 54,675 unique probe IDs. The source annotation is carried into this report without replacing it with a gene-symbol-only match.

The source analysis documentation describes the contrast as **cancer tissue minus normal breast tissue**, with **104 cancer samples and 17 normal samples**. It states the GEO2R logFC direction is positive for higher cancer expression and negative for lower cancer expression. GSE42568 adjusted p-values are the source table's `adj.P.Val` values, with the documentation describing Benjamini-Hochberg adjustment across the GEO2R probe tests.

GSE65194 values are taken from `GSE65194_tumor_vs_healthy_welch.tsv`: mean tumour minus healthy-control log2 expression difference and its BH-adjusted p-value. Sample sizes are 130 tumour patients and 11 healthy controls.

## Exact-probe matching and summary

Probe rows were joined using the **exact probe ID** only. The target panel yielded **38 exact shared probes**. Gene labels are included for readability and do not determine the match. The GSE42568 `Gene.symbol` annotation is retained exactly; `217492_s_at` remains `PTENP1///PTEN`.

- Concordant direction: **24** probes.
- Discordant direction: **14** probes.
- Neutral direction comparison (zero difference in at least one cohort): **0** probes.
- Significant in both cohorts (adjusted p < 0.05 in each): **15** probes.
- Significant only in GSE42568: **3** probes.
- Significant only in GSE65194: **12** probes.
- Significant in neither: **8** probes.

### Counts by gene

| Gene | Matching probes | Concordant | Discordant | Significant in both | GSE42568 only | GSE65194 only | Neither |
|---|---:|---:|---:|---:|---:|---:|---:|
| BRCA1 | 2 | 2 | 0 | 0 | 0 | 1 | 1 |
| BRCA2 | 2 | 2 | 0 | 1 | 0 | 1 | 0 |
| TP53 | 2 | 1 | 1 | 1 | 0 | 1 | 0 |
| ERBB2 | 3 | 3 | 0 | 1 | 0 | 2 | 0 |
| PIK3CA | 3 | 2 | 1 | 3 | 0 | 0 | 0 |
| ESR1 | 9 | 5 | 4 | 1 | 0 | 4 | 4 |
| PTEN | 8 | 5 | 3 | 5 | 1 | 1 | 1 |
| EGFR | 9 | 4 | 5 | 3 | 2 | 2 | 2 |

### Probe-level comparison

Direction is compared with both contrasts oriented cancer/tumour minus normal/healthy. Significance means adjusted p < 0.05 in that cohort's existing analysis.

| Probe ID | Gene | GSE42568 annotation | GSE42568 logFC | GSE42568 adjusted p | GSE65194 tumour-minus-healthy mean difference | GSE65194 adjusted p | Direction | GSE42568 significant | GSE65194 significant |
|---|---|---|---:|---:|---:|---:|---|---|---|
| 204531_s_at | BRCA1 | BRCA1 | 0.357 | 0.17 | 0.4305 | 0.0526298 | concordant | no | no |
| 211851_x_at | BRCA1 | BRCA1 | 0.00335 | 0.979 | 0.156438 | 0.00439897 | concordant | no | yes |
| 208368_s_at | BRCA2 | BRCA2 | 0.536 | 0.00151 | 0.752426 | 0.00121237 | concordant | yes | yes |
| 214727_at | BRCA2 | BRCA2 | 0.019 | 0.654 | 0.205522 | 0.0122037 | concordant | no | yes |
| 201746_at | TP53 | TP53 | 0.671 | 0.0167 | -0.825305 | 0.000536498 | discordant | yes | yes |
| 211300_s_at | TP53 | TP53 | 0.0371 | 0.688 | 1.22475 | 1.86173e-07 | concordant | no | yes |
| 210930_s_at | ERBB2 | ERBB2 | 0.0683 | 0.652 | 2.09064 | 4.3154e-09 | concordant | no | yes |
| 216836_s_at | ERBB2 | ERBB2 | 3.12 | 3.52e-08 | 1.88627 | 2.03754e-06 | concordant | yes | yes |
| 234354_x_at | ERBB2 | ERBB2 | 0.499 | 0.162 | 1.21077 | 2.75794e-07 | concordant | no | yes |
| 204369_at | PIK3CA | PIK3CA | -1.47 | 1.89e-08 | 2.93926 | 0.000924865 | discordant | yes | yes |
| 231854_at | PIK3CA | PIK3CA | -1.29 | 3.7e-16 | -0.392799 | 0.0196281 | concordant | yes | yes |
| 235980_at | PIK3CA | PIK3CA | -1.08 | 9.62e-09 | -1.3486 | 0.000401078 | concordant | yes | yes |
| 205225_at | ESR1 | ESR1 | 1.75 | 0.105 | -0.92931 | 0.186141 | discordant | no | no |
| 211233_x_at | ESR1 | ESR1 | 0.0555 | 0.648 | 0.144858 | 0.0014836 | concordant | no | yes |
| 211234_x_at | ESR1 | ESR1 | 0.0612 | 0.429 | 0.152416 | 0.000924865 | concordant | no | yes |
| 211235_s_at | ESR1 | ESR1 | 0.114 | 0.414 | 0.344181 | 0.00121237 | concordant | no | yes |
| 211627_x_at | ESR1 | ESR1 | 0.103 | 0.722 | -0.879169 | 0.00103382 | discordant | no | yes |
| 215551_at | ESR1 | ESR1 | 0.0348 | 0.831 | -0.0252843 | 0.733609 | discordant | no | no |
| 215552_s_at | ESR1 | ESR1 | 0.924 | 0.0419 | 0.674832 | 0.00193378 | concordant | yes | yes |
| 217163_at | ESR1 | ESR1 | 0.00127 | 0.966 | -0.00711654 | 0.169517 | discordant | no | no |
| 217190_x_at | ESR1 | ESR1 | 0.17 | 0.18 | 0.0531932 | 0.0879195 | concordant | no | no |
| 204053_x_at | PTEN | PTEN | -2.03 | 9.78e-13 | 0.379938 | 0.186141 | discordant | yes | no |
| 204054_at | PTEN | PTEN | -1.12 | 3.73e-08 | -0.306884 | 0.00688804 | concordant | yes | yes |
| 211711_s_at | PTEN | PTEN | -2.37 | 4.41e-15 | -0.95064 | 0.00121237 | concordant | yes | yes |
| 217492_s_at | PTEN | PTENP1///PTEN | -2.02 | 1.66e-13 | 0.711807 | 0.0420454 | discordant | yes | yes |
| 225363_at | PTEN | PTEN | -1.21 | 9.55e-05 | -1.00955 | 3.98045e-06 | concordant | yes | yes |
| 228006_at | PTEN | PTEN | -1.98 | 4.34e-12 | -1.59645 | 4.3154e-09 | concordant | yes | yes |
| 233254_x_at | PTEN | PTEN | 0.241 | 0.375 | 0.00653324 | 0.975686 | concordant | no | no |
| 242622_x_at | PTEN | PTEN | 0.348 | 0.0853 | -0.615822 | 0.00449811 | discordant | no | yes |
| 1565483_at | EGFR | EGFR | 0.116 | 0.826 | -0.0380935 | 0.464233 | discordant | no | no |
| 1565484_x_at | EGFR | EGFR | 0.163 | 0.547 | -0.218496 | 0.172927 | discordant | no | no |
| 201983_s_at | EGFR | EGFR | -2.51 | 3.89e-05 | -1.89773 | 0.00116396 | concordant | yes | yes |
| 201984_s_at | EGFR | EGFR | -2.09 | 9.35e-08 | -0.993656 | 0.0196281 | concordant | yes | yes |
| 210984_x_at | EGFR | EGFR | -0.357 | 0.000592 | 0.0960042 | 0.210111 | discordant | yes | no |
| 211550_at | EGFR | EGFR | 0.00226 | 0.928 | -0.151177 | 0.0014836 | discordant | no | yes |
| 211551_at | EGFR | EGFR | -0.111 | 0.179 | -0.0262287 | 0.0373464 | concordant | no | yes |
| 211607_x_at | EGFR | EGFR | -0.226 | 0.0445 | 0.0778203 | 0.278278 | discordant | yes | no |
| 224999_at | EGFR | EGFR | -2.39 | 7.14e-08 | -1.95798 | 0.000532736 | concordant | yes | yes |

## Interpretation limits

- **Different sample sizes:** GSE42568 uses 104 cancer and 17 normal samples; GSE65194 uses 130 tumour patients and 11 healthy controls. Different group sizes affect precision and power.
- **Different cohort composition:** group definitions and patient composition may differ. The comparison follows the documented cancer/tumour versus normal/healthy contrasts; it does not establish perfectly equivalent populations.
- **Different preprocessing:** both datasets report GC-RMA-related processing, but GSE65194 also documents batch/hybridization correction. The v1 audit does not document the same batch correction for GSE42568. The exact GSE65194 CDF and explicit transformation history are unresolved, so scales and preprocessing should not be assumed identical.
- **Technical replicates:** GSE65194's tested tumour observations are patient-level values created by the project by arithmetic-averaging its validated replicate pairs. GSE42568's GEO2R table is the pre-existing source result; it was not reconstructed or reprocessed here.
- **Different statistical methods:** GSE42568 results were produced with GEO2R/limma; GSE65194 uses Welch's independent two-sample t-test. The effect estimates and p-values are therefore not method-identical.
- **Different multiple-testing families:** GSE42568 uses the supplied `adj.P.Val` from its GEO2R analysis across the broader probe table; GSE65194 adjusted its 38 selected probes as a target-panel family. Significance counts are useful descriptions but are not directly equivalent tests of replication strength.
- **Probe annotations:** exact probe IDs were required. Annotation strings, including the ambiguous PTENP1/PTEN entry, are preserved; probe-level differences should not be collapsed into gene-level claims.
- **Exploratory evidence:** concordant direction or adjusted significance is not proof of biological validation, causality, or diagnostic performance.

## Provenance and scope

The GSE42568 table was supplied as a pre-existing result; its full source table and supporting report are not included in this repository. No result was reconstructed or invented, no matrix was pooled across cohorts, and no new tests were run for this audit. The existing GSE65194 targeted Welch results were compared as supplied.
