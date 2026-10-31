# GSE65194 local limma sensitivity analysis

## Validation and reproducibility

- R version: 4.6.2 (x86_64-w64-mingw32).
- limma version: 3.68.5.
- Input GEO matrix: 54,673 probes and 178 source GSM columns; the expression table was read between the GEO series-matrix table markers.
- Assignment validation: each sheet contains 141 unique GSMs (130 Tumour and 11 Healthy), representing 130 unique tumour patients and 11 uniquely identified healthy controls.
- repA and repB differ only in the selected arrays for 23 paired tumour patients; 107 singleton tumour selections and all 11 healthy controls are unchanged.
- Group sizes per fit: Tumour = 130; Healthy = 11.
- Contrast: Tumour - Healthy. Positive logFC means higher expression in Tumour.
- The matrix contains processed log2 expression values; no further log transformation, normalization, batch correction, sample removal, or probe filtering was applied.
- repA and repB are alternative technical-replicate selections, not independent cohorts.

## Methods

Both full-probe analyses used `limma::lmFit()`, `limma::contrasts.fit()`, `limma::eBayes()`, and `limma::topTable()`, with BH adjustment across all 54,673 probes. The full results retain all probes.

The comparison with the prior Welch analysis uses exact probe IDs. Welch BH-adjusted p-values were calculated across its 38 target probes, while limma BH-adjusted p-values were calculated across all 54,673 probes; significance differences should be interpreted with these different correction families in mind. Significance is adjusted p < 0.05.

This is a local limma analysis of the downloaded GEO matrix using the planned sample assignments. It does not claim exact reproduction of the GEO2R website.

## Limma versus Welch: probe-level comparison

- Directional agreement with Welch: repA 38/38; repB 38/38.
- Adjusted-significance comparison (both / Welch only / limma only / neither): repA 18 / 9 / 1 / 10; repB 18 / 9 / 1 / 10.

| Probe ID | Gene | Welch tumour - healthy | limma repA logFC | limma repB logFC | Welch BH p | repA BH p | repB BH p | Direction agrees repA | Direction agrees repB | Welch significant | repA significant | repB significant |
|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|---|
| 204531_s_at | BRCA1 | 0.4305 | 0.4316 | 0.4294 | 0.05263 | 0.2072 | 0.2097 | Yes | Yes | No | No | No |
| 211851_x_at | BRCA1 | 0.1564 | 0.1565 | 0.1564 | 0.004399 | 0.3319 | 0.3622 | Yes | Yes | Yes | No | No |
| 208368_s_at | BRCA2 | 0.7524 | 0.751 | 0.7538 | 0.001212 | 0.01268 | 0.01364 | Yes | Yes | Yes | Yes | Yes |
| 214727_at | BRCA2 | 0.2055 | 0.2026 | 0.2085 | 0.0122 | 0.3693 | 0.3662 | Yes | Yes | Yes | No | No |
| 201746_at | TP53 | -0.8253 | -0.8251 | -0.8255 | 0.0005365 | 0.0375 | 0.03774 | Yes | Yes | Yes | Yes | Yes |
| 211300_s_at | TP53 | 1.225 | 1.225 | 1.225 | 1.862e-07 | 0.004015 | 0.003946 | Yes | Yes | Yes | Yes | Yes |
| 210930_s_at | ERBB2 | 2.091 | 2.086 | 2.096 | 4.315e-09 | 0.02761 | 0.03026 | Yes | Yes | Yes | Yes | Yes |
| 216836_s_at | ERBB2 | 1.886 | 1.886 | 1.886 | 2.038e-06 | 0.02097 | 0.02153 | Yes | Yes | Yes | Yes | Yes |
| 234354_x_at | ERBB2 | 1.211 | 1.208 | 1.213 | 2.758e-07 | 0.0594 | 0.06157 | Yes | Yes | Yes | No | No |
| 204369_at | PIK3CA | 2.939 | 2.938 | 2.941 | 0.0009249 | 1.746e-16 | 9.683e-17 | Yes | Yes | Yes | Yes | Yes |
| 231854_at | PIK3CA | -0.3928 | -0.4007 | -0.3849 | 0.01963 | 1.674e-06 | 0.0003339 | Yes | Yes | Yes | Yes | Yes |
| 235980_at | PIK3CA | -1.349 | -1.352 | -1.345 | 0.0004011 | 5.244e-12 | 4.851e-12 | Yes | Yes | Yes | Yes | Yes |
| 205225_at | ESR1 | -0.9293 | -0.9277 | -0.931 | 0.1861 | 0.5713 | 0.5721 | Yes | Yes | No | No | No |
| 211233_x_at | ESR1 | 0.1449 | 0.1449 | 0.1449 | 0.001484 | 0.3967 | 0.3988 | Yes | Yes | Yes | No | No |
| 211234_x_at | ESR1 | 0.1524 | 0.1524 | 0.1524 | 0.0009249 | 0.3466 | 0.3483 | Yes | Yes | Yes | No | No |
| 211235_s_at | ESR1 | 0.3442 | 0.3442 | 0.3442 | 0.001212 | 0.2781 |  0.28 | Yes | Yes | Yes | No | No |
| 211627_x_at | ESR1 | -0.8792 | -0.8789 | -0.8794 | 0.001034 | 2.033e-20 | 4.753e-21 | Yes | Yes | Yes | Yes | Yes |
| 215551_at | ESR1 | -0.02528 | -0.02528 | -0.02529 | 0.7336 | 0.8886 | 0.8904 | Yes | Yes | No | No | No |
| 215552_s_at | ESR1 | 0.6748 | 0.6748 | 0.6748 | 0.001934 | 0.1986 | 0.1997 | Yes | Yes | Yes | No | No |
| 217163_at | ESR1 | -0.007117 | -0.007078 | -0.007155 | 0.1695 | 0.1432 | 0.136 | Yes | Yes | No | No | No |
| 217190_x_at | ESR1 | 0.05319 | 0.0532 | 0.05319 | 0.08792 | 0.6151 | 0.6179 | Yes | Yes | No | No | No |
| 204053_x_at | PTEN | 0.3799 | 0.3802 | 0.3797 | 0.1861 | 0.2002 | 0.1917 | Yes | Yes | No | No | No |
| 204054_at | PTEN | -0.3069 | -0.3066 | -0.3072 | 0.006888 | 0.2326 | 0.2517 | Yes | Yes | Yes | No | No |
| 211711_s_at | PTEN | -0.9506 | -0.95 | -0.9512 | 0.001212 | 0.002491 | 0.003098 | Yes | Yes | Yes | Yes | Yes |
| 217492_s_at | PTEN | 0.7118 | 0.712 | 0.7116 | 0.04205 | 0.01157 | 0.009763 | Yes | Yes | Yes | Yes | Yes |
| 225363_at | PTEN | -1.01 | -1.009 | -1.01 | 3.98e-06 | 7.813e-05 | 5.771e-05 | Yes | Yes | Yes | Yes | Yes |
| 228006_at | PTEN | -1.596 | -1.596 | -1.597 | 4.315e-09 | 1.535e-05 | 1.651e-05 | Yes | Yes | Yes | Yes | Yes |
| 233254_x_at | PTEN | 0.006533 | 0.007647 | 0.005419 | 0.9757 | 0.9786 | 0.984 | Yes | Yes | No | No | No |
| 242622_x_at | PTEN | -0.6158 | -0.6106 | -0.621 | 0.004498 | 1.706e-06 | 5.22e-06 | Yes | Yes | Yes | Yes | Yes |
| 1565483_at | EGFR | -0.03809 | -0.04219 | -0.034 | 0.4642 | 0.5231 | 0.665 | Yes | Yes | No | No | No |
| 1565484_x_at | EGFR | -0.2185 | -0.2219 | -0.2151 | 0.1729 | 7.405e-07 | 0.0002265 | Yes | Yes | No | Yes | Yes |
| 201983_s_at | EGFR | -1.898 | -1.897 | -1.898 | 0.001164 | 0.007189 | 0.007076 | Yes | Yes | Yes | Yes | Yes |
| 201984_s_at | EGFR | -0.9937 | -0.9931 | -0.9942 | 0.01963 | 0.06502 | 0.06821 | Yes | Yes | Yes | No | No |
| 210984_x_at | EGFR | 0.096 | 0.09637 | 0.09564 | 0.2101 | 0.6798 | 0.6882 | Yes | Yes | No | No | No |
| 211550_at | EGFR | -0.1512 | -0.1515 | -0.1508 | 0.001484 | 8.242e-08 | 5.627e-09 | Yes | Yes | Yes | Yes | Yes |
| 211551_at | EGFR | -0.02623 | -0.02622 | -0.02624 | 0.03735 | 3.033e-05 | 2.318e-05 | Yes | Yes | Yes | Yes | Yes |
| 211607_x_at | EGFR | 0.07782 | 0.07821 | 0.07743 | 0.2783 | 0.7228 | 0.7315 | Yes | Yes | No | No | No |
| 224999_at | EGFR | -1.958 | -1.958 | -1.958 | 0.0005327 | 0.0003464 | 0.0002528 | Yes | Yes | Yes | Yes | Yes |
