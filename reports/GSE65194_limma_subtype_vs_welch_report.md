# GSE65194 patient-level limma subtype analysis

## Validation and reproducibility

- R version: 4.6.2 (x86_64-w64-mingw32).
- limma version: 3.68.5.
- Expression input: existing patient-level averaged matrix; no re-averaging, log transformation, normalization, or batch correction was performed.
- The patient-level matrix contains 54,673 probes × 141 observations. The analysis selected the 130 tumour patients only; 11 healthy controls were excluded. The 14 cell-line arrays are not in the patient-level matrix.
- Tumour subtype counts: TNBC 41, Her2 30, Luminal A 29, Luminal B 30.
- All 38 target probes were validated against the compatibility audit and were present in the matrix. Probe measurements remain separate; no gene-level averaging was performed.
- Design group order: TNBC, Her2, Luminal A, Luminal B.
- Positive logFC means higher expression in the first-named subtype.

## Model and multiple testing

The model used `limma::lmFit()`, `limma::makeContrasts()`, `limma::contrasts.fit()`, and `limma::eBayes()`. Six pairwise contrasts were extracted for each of the 38 target probes. Benjamini-Hochberg adjustment was applied once across all 228 limma raw p-values; the Welch adjusted p-values use the existing complete 228-test family.

## Comparison with Welch

- Across all 228 probe-by-contrast results, limma has 100 significant probes and Welch has 95.
- Significant in both methods: 87; limma only: 13; Welch only: 8; neither: 120.
- Directional agreement: 228/228 probe-by-contrast results.

| Contrast | Limma significant | Welch significant | Both | Limma only | Welch only | Neither | Direction agrees |
|---|---:|---:|---:|---:|---:|---:|---:|
| TNBC - Her2 | 18 | 12 | 12 | 6 | 0 | 20 | 38/38 |
| TNBC - Luminal A | 24 | 19 | 19 | 5 | 0 | 14 | 38/38 |
| TNBC - Luminal B | 21 | 22 | 20 | 1 | 2 | 15 | 38/38 |
| Her2 - Luminal A | 17 | 17 | 16 | 1 | 1 | 20 | 38/38 |
| Her2 - Luminal B | 13 | 15 | 13 | 0 | 2 | 23 | 38/38 |
| Luminal A - Luminal B | 7 | 10 | 7 | 0 | 3 | 28 | 38/38 |

All comparisons below are probe-level and joined using exact probe IDs and contrast names. LogFC difference is limma minus Welch. Adjusted-p differences compare the same 228-test Welch family against limma BH correction across all 228 results.

| Contrast | Probe ID | Annotation | limma logFC | Welch difference | logFC difference | limma raw p | Welch raw p | limma BH p | Welch BH p | Direction agrees | Significance |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| TNBC - Her2 | 204531_s_at | BRCA1 | 0.12271 | 0.12271 | 9.2287e-15 | 0.58064 | 0.60886 | 0.72739 | 0.70442 | Yes | Neither |
| TNBC - Her2 | 211851_x_at | BRCA1 | 0.23328 | 0.23328 | -2.911e-13 | 0.032423 | 0.041652 | 0.070403 | 0.084041 | Yes | Neither |
| TNBC - Her2 | 208368_s_at | BRCA2 | 0.23073 | 0.23073 | 1.7389e-13 | 0.21836 | 0.27572 | 0.34815 | 0.4077 | Yes | Neither |
| TNBC - Her2 | 214727_at | BRCA2 | 0.043489 | 0.043489 | -1.6639e-14 | 0.76644 | 0.80895 | 0.85351 | 0.87369 | Yes | Neither |
| TNBC - Her2 | 201746_at | TP53 | 0.68254 | 0.68254 | 8.5931e-14 | 0.015094 | 0.04861 | 0.036641 | 0.094807 | Yes | Limma only |
| TNBC - Her2 | 211300_s_at | TP53 | 0.74829 | 0.74829 | 1.3123e-13 | 0.011586 | 0.030673 | 0.030018 | 0.0676 | Yes | Limma only |
| TNBC - Her2 | 210930_s_at | ERBB2 | -6.2394 | -6.2394 | 2.5997e-12 | 5.4328e-40 | 4.6781e-33 | 6.1934e-38 | 5.3331e-31 | Yes | Both |
| TNBC - Her2 | 216836_s_at | ERBB2 | -5.1905 | -5.1905 | 4.964e-12 | 6.9866e-41 | 4.6779e-33 | 1.593e-38 | 5.3331e-31 | Yes | Both |
| TNBC - Her2 | 234354_x_at | ERBB2 | -3.5398 | -3.5398 | -1.5423e-12 | 1.1092e-22 | 9.4926e-13 | 2.5289e-21 | 1.6649e-11 | Yes | Both |
| TNBC - Her2 | 204369_at | PIK3CA | 0.61677 | 0.61677 | -3.8303e-13 | 0.0020326 | 0.0023864 | 0.0068152 | 0.0080016 | Yes | Both |
| TNBC - Her2 | 231854_at | PIK3CA | 0.16965 | 0.16965 | 1.1349e-13 | 0.0041566 | 0.017637 | 0.011846 | 0.043574 | Yes | Both |
| TNBC - Her2 | 235980_at | PIK3CA | 0.31524 | 0.31524 | -1.6648e-13 | 0.011406 | 0.024852 | 0.029893 | 0.057819 | Yes | Limma only |
| TNBC - Her2 | 205225_at | ESR1 | -0.34682 | -0.34682 | 1.8635e-13 | 0.41074 | 0.46764 | 0.57102 | 0.58263 | Yes | Neither |
| TNBC - Her2 | 211233_x_at | ESR1 | -0.00099352 | -0.00099352 | 7.223e-16 | 0.99252 | 0.57898 | 0.99791 | 0.68168 | Yes | Neither |
| TNBC - Her2 | 211234_x_at | ESR1 | -0.0024558 | -0.0024558 | -2.8475e-15 | 0.98019 | 0.28059 | 0.99791 | 0.40888 | Yes | Neither |
| TNBC - Her2 | 211235_s_at | ESR1 | -0.0037588 | -0.0037588 | 6.2016e-16 | 0.98343 | 0.68846 | 0.99791 | 0.7657 | Yes | Neither |
| TNBC - Her2 | 211627_x_at | ESR1 | -0.0014374 | -0.0014374 | -2.0359e-15 | 0.97599 | 0.95846 | 0.99791 |  0.989 | Yes | Neither |
| TNBC - Her2 | 215551_at | ESR1 | 0.0060254 | 0.0060254 | 5.0116e-15 | 0.95176 | 0.34133 | 0.99791 | 0.46049 | Yes | Neither |
| TNBC - Her2 | 215552_s_at | ESR1 | -0.016323 | -0.016323 | 4.1023e-14 | 0.95102 | 0.38907 | 0.99791 | 0.50981 | Yes | Neither |
| TNBC - Her2 | 217163_at | ESR1 | -0.00083726 | -0.00083726 | -3.7058e-16 | 0.77043 | 0.68034 | 0.85351 | 0.7639 | Yes | Neither |
| TNBC - Her2 | 217190_x_at | ESR1 | -0.038753 | -0.038753 | 2.4765e-14 | 0.54296 | 0.14725 | 0.69159 | 0.2398 | Yes | Neither |
| TNBC - Her2 | 204053_x_at | PTEN | -0.74019 | -0.74019 | 9.7811e-14 | 9.8311e-05 | 0.00066979 | 0.0004666 | 0.0029368 | Yes | Both |
| TNBC - Her2 | 204054_at | PTEN | -0.34722 | -0.34722 | -2.3093e-14 | 0.049914 | 0.085052 | 0.098107 | 0.14917 | Yes | Neither |
| TNBC - Her2 | 211711_s_at | PTEN | -0.72128 | -0.72128 | 1.3145e-13 | 0.0011318 | 0.0028247 | 0.0041622 | 0.009071 | Yes | Both |
| TNBC - Her2 | 217492_s_at | PTENP1 /// PTEN | -0.68563 | -0.68563 | -3.0553e-13 | 0.00025725 | 0.0010742 | 0.0010113 | 0.004151 | Yes | Both |
| TNBC - Her2 | 225363_at | PTEN | -0.70185 | -0.70185 | -1.2146e-13 | 1.2424e-06 | 2.7954e-05 | 7.8686e-06 | 0.00016342 | Yes | Both |
| TNBC - Her2 | 228006_at | PTEN | -0.84016 | -0.84016 | -1.3645e-13 | 0.0001196 | 0.00091854 | 0.00053363 | 0.0037398 | Yes | Both |
| TNBC - Her2 | 233254_x_at | PTEN | -0.55615 | -0.55615 | 4.5375e-13 | 6.3404e-05 | 0.00032768 | 0.00031427 | 0.0015565 | Yes | Both |
| TNBC - Her2 | 242622_x_at | PTEN | -0.22621 | -0.22621 | -1.8674e-13 | 0.0083725 | 0.022859 | 0.022999 | 0.053731 | Yes | Limma only |
| TNBC - Her2 | 1565483_at | EGFR | -0.021178 | -0.021178 | 3.0222e-14 | 0.62286 | 0.32586 | 0.76764 | 0.45028 | Yes | Neither |
| TNBC - Her2 | 1565484_x_at | EGFR | -0.0050872 | -0.0050872 | 3.2621e-15 | 0.76879 | 0.82685 | 0.85351 | 0.87685 | Yes | Neither |
| TNBC - Her2 | 201983_s_at | EGFR | 1.0618 | 1.0618 | 4.9405e-13 | 0.0096918 | 0.019846 | 0.025997 | 0.047629 | Yes | Both |
| TNBC - Her2 | 201984_s_at | EGFR | 0.73227 | 0.73227 | 2.1649e-14 | 0.015769 | 0.077642 | 0.037846 | 0.13723 | Yes | Limma only |
| TNBC - Her2 | 210984_x_at | EGFR | -0.14424 | -0.14424 | 1.7103e-13 | 0.32133 | 0.5121 | 0.46964 | 0.62105 | Yes | Neither |
| TNBC - Her2 | 211550_at | EGFR | 0.0060883 | 0.0060883 | -4.146e-16 | 0.72452 | 0.56509 | 0.84713 | 0.67798 | Yes | Neither |
| TNBC - Her2 | 211551_at | EGFR | -0.00053831 | -0.00053831 | 2.6591e-15 | 0.88043 | 0.81741 | 0.95589 | 0.87369 | Yes | Neither |
| TNBC - Her2 | 211607_x_at | EGFR | -0.1125 | -0.1125 | 2.4578e-14 | 0.40817 | 0.57875 | 0.57094 | 0.68168 | Yes | Neither |
| TNBC - Her2 | 224999_at | EGFR | 0.81338 | 0.81338 | -3.0698e-13 | 0.015106 | 0.043017 | 0.036641 | 0.086035 | Yes | Limma only |
| TNBC - Luminal A | 204531_s_at | BRCA1 | 0.81581 | 0.81581 | -1.088e-13 | 0.00038726 | 0.0011502 | 0.0014965 | 0.0043708 | Yes | Both |
| TNBC - Luminal A | 211851_x_at | BRCA1 | 0.31917 | 0.31917 | -2.9282e-13 | 0.0040183 | 0.0042181 | 0.011597 | 0.012022 | Yes | Both |
| TNBC - Luminal A | 208368_s_at | BRCA2 | 1.2358 | 1.2358 | -9.6723e-13 | 1.2361e-09 | 2.5492e-10 | 1.0438e-08 | 3.229e-09 | Yes | Both |
| TNBC - Luminal A | 214727_at | BRCA2 | 0.44841 | 0.44841 | -2.8139e-13 | 0.002882 | 0.001176 | 0.0092549 | 0.0043956 | Yes | Both |
| TNBC - Luminal A | 201746_at | TP53 | 0.20638 | 0.20638 | 4.7878e-13 | 0.46211 | 0.41131 | 0.61849 | 0.53119 | Yes | Neither |
| TNBC - Luminal A | 211300_s_at | TP53 | 0.64834 | 0.64834 | -4.6518e-13 | 0.029781 | 0.033768 | 0.065924 | 0.071288 | Yes | Neither |
| TNBC - Luminal A | 210930_s_at | ERBB2 | -0.32647 | -0.32647 | 2.3331e-13 | 0.31396 | 0.070068 | 0.46182 | 0.12883 | Yes | Neither |
| TNBC - Luminal A | 216836_s_at | ERBB2 | -1.0428 | -1.0428 | -1.3645e-12 | 0.0001217 | 4.7104e-05 | 0.00053363 | 0.00025571 | Yes | Both |
| TNBC - Luminal A | 234354_x_at | ERBB2 | 0.00078 | 0.00078 | -7.8063e-17 | 0.99791 | 0.99095 | 0.99791 | 0.99554 | Yes | Neither |
| TNBC - Luminal A | 204369_at | PIK3CA | 0.35811 | 0.35811 | 3.6321e-13 | 0.07245 | 0.031948 | 0.13652 | 0.068719 | Yes | Neither |
| TNBC - Luminal A | 231854_at | PIK3CA | 0.092696 | 0.092696 | -3.024e-14 | 0.11677 | 0.12817 | 0.20169 | 0.2133 | Yes | Neither |
| TNBC - Luminal A | 235980_at | PIK3CA | -0.039681 | -0.039681 | 3.8448e-14 | 0.74949 | 0.73333 | 0.85351 |    0.8 | Yes | Neither |
| TNBC - Luminal A | 205225_at | ESR1 | -8.2106 | -8.2106 | 4.8193e-12 | 1.2207e-39 | 4.8663e-30 | 9.2771e-38 | 3.6984e-28 | Yes | Both |
| TNBC - Luminal A | 211233_x_at | ESR1 | -0.45379 | -0.45379 | -2.2332e-13 | 4.1001e-05 | 0.0037575 | 0.00020774 | 0.011273 | Yes | Both |
| TNBC - Luminal A | 211234_x_at | ESR1 | -0.46958 | -0.46958 | 4.8572e-13 | 6.3992e-06 | 0.0016123 | 3.6475e-05 | 0.005835 | Yes | Both |
| TNBC - Luminal A | 211235_s_at | ESR1 | -1.1887 | -1.1887 | 2.1605e-12 | 1.537e-09 | 2.5815e-05 | 1.2084e-08 | 0.00015907 | Yes | Both |
| TNBC - Luminal A | 211627_x_at | ESR1 | -0.055301 | -0.055301 | 1.8367e-14 | 0.25281 | 0.40126 | 0.38684 | 0.52279 | Yes | Neither |
| TNBC - Luminal A | 215551_at | ESR1 | -0.2382 | -0.2382 | 2.9304e-13 | 0.019146 | 0.072515 | 0.043652 | 0.13122 | Yes | Limma only |
| TNBC - Luminal A | 215552_s_at | ESR1 | -2.2449 | -2.2449 | -3.1211e-12 | 8.5359e-14 | 1.8406e-07 | 1.1448e-12 | 1.6101e-06 | Yes | Both |
| TNBC - Luminal A | 217163_at | ESR1 | 0.004899 | 0.004899 | -4.6673e-15 | 0.092655 | 0.06568 |  0.169 | 0.12217 | Yes | Neither |
| TNBC - Luminal A | 217190_x_at | ESR1 | -0.24543 | -0.24543 | -4.7151e-13 | 0.00020401 | 0.011549 | 0.0008457 | 0.029922 | Yes | Both |
| TNBC - Luminal A | 204053_x_at | PTEN | -0.44857 | -0.44857 | -3.8725e-13 | 0.017214 | 0.0285 | 0.040883 | 0.06498 | Yes | Limma only |
| TNBC - Luminal A | 204054_at | PTEN | -0.42092 | -0.42092 | 2.1738e-13 | 0.018986 | 0.025232 | 0.043652 | 0.05811 | Yes | Limma only |
| TNBC - Luminal A | 211711_s_at | PTEN | -0.70909 | -0.70909 | 4.6574e-13 | 0.0015113 | 0.0031381 | 0.0052843 | 0.0098012 | Yes | Both |
| TNBC - Luminal A | 217492_s_at | PTENP1 /// PTEN | -0.45975 | -0.45975 | 1.3162e-13 | 0.013795 | 0.022299 | 0.034947 | 0.05296 | Yes | Limma only |
| TNBC - Luminal A | 225363_at | PTEN | -1.2859 | -1.2859 | -4.0898e-12 | 7.3774e-16 | 3.9958e-12 | 1.2015e-14 | 6.0736e-11 | Yes | Both |
| TNBC - Luminal A | 228006_at | PTEN | -1.7514 | -1.7514 | 4.613e-12 | 2.3544e-13 | 2.8163e-10 | 2.8253e-12 | 3.3795e-09 | Yes | Both |
| TNBC - Luminal A | 233254_x_at | PTEN | -0.77588 | -0.77588 | 1.2179e-13 | 7.3982e-08 | 1.4986e-06 | 5.1115e-07 | 1.1022e-05 | Yes | Both |
| TNBC - Luminal A | 242622_x_at | PTEN | -0.18466 | -0.18466 | 3.8808e-13 | 0.032248 | 0.036239 | 0.070403 | 0.074436 | Yes | Neither |
| TNBC - Luminal A | 1565483_at | EGFR | -0.075298 | -0.075298 | -3.0212e-14 | 0.085054 | 0.2688 | 0.15639 | 0.40321 | Yes | Neither |
| TNBC - Luminal A | 1565484_x_at | EGFR | 0.00036894 | 0.00036894 | -2.4899e-16 | 0.98316 | 0.97796 | 0.99791 | 0.99554 | Yes | Neither |
| TNBC - Luminal A | 201983_s_at | EGFR | 1.8761 | 1.8761 | 1.8714e-12 | 1.0298e-05 | 1.8685e-06 | 5.7266e-05 | 1.291e-05 | Yes | Both |
| TNBC - Luminal A | 201984_s_at | EGFR | 2.1202 | 2.1202 | -4.4227e-12 | 1.2239e-10 | 3.615e-12 | 1.1162e-09 | 5.8873e-11 | Yes | Both |
| TNBC - Luminal A | 210984_x_at | EGFR | 0.066801 | 0.066801 | 2.5924e-14 | 0.64876 |  0.365 | 0.78908 | 0.48667 | Yes | Neither |
| TNBC - Luminal A | 211550_at | EGFR | 0.044894 | 0.044894 | 3.9413e-14 | 0.011052 | 0.033364 | 0.029299 | 0.071093 | Yes | Limma only |
| TNBC - Luminal A | 211551_at | EGFR | 0.011024 | 0.011024 | -1.9956e-14 | 0.0027306 | 0.0062106 | 0.008894 | 0.016857 | Yes | Both |
| TNBC - Luminal A | 211607_x_at | EGFR | 0.071473 | 0.071473 | 2.6992e-14 | 0.60256 | 0.36406 | 0.75073 | 0.48667 | Yes | Neither |
| TNBC - Luminal A | 224999_at | EGFR |   1.26 |   1.26 | -2.0066e-12 | 0.00024166 | 4.417e-05 | 0.00096662 | 0.00024563 | Yes | Both |
| TNBC - Luminal B | 204531_s_at | BRCA1 | 0.10723 | 0.10723 | 1.0977e-14 | 0.62922 | 0.64602 | 0.77131 | 0.7328 | Yes | Neither |
| TNBC - Luminal B | 211851_x_at | BRCA1 | 0.22914 | 0.22914 | 3.7401e-13 | 0.03556 | 0.087771 | 0.075071 | 0.15276 | Yes | Neither |
| TNBC - Luminal B | 208368_s_at | BRCA2 | 0.63741 | 0.63741 | 1.722e-13 | 0.00085019 | 0.0016838 | 0.0031778 | 0.005929 | Yes | Both |
| TNBC - Luminal B | 214727_at | BRCA2 | 0.28628 | 0.28628 | 3.5033e-13 | 0.052241 | 0.063764 | 0.1018 | 0.12015 | Yes | Neither |
| TNBC - Luminal B | 201746_at | TP53 | 0.58863 | 0.58863 | 8.8041e-14 | 0.035555 | 0.031824 | 0.075071 | 0.068719 | Yes | Neither |
| TNBC - Luminal B | 211300_s_at | TP53 | 1.1735 | 1.1735 | -1.8694e-12 | 0.00010028 | 0.00013647 | 0.0004666 | 0.00072363 | Yes | Both |
| TNBC - Luminal B | 210930_s_at | ERBB2 | -3.3937 | -3.3937 | 2.6015e-12 | 3.2342e-19 | 3.1791e-09 | 5.6723e-18 | 3.6242e-08 | Yes | Both |
| TNBC - Luminal B | 216836_s_at | ERBB2 | -4.2681 | -4.2681 | 1.6271e-12 | 4.1661e-33 | 2.029e-18 | 1.1873e-31 | 5.1401e-17 | Yes | Both |
| TNBC - Luminal B | 234354_x_at | ERBB2 | -1.888 | -1.888 | 1.7892e-12 | 2.5982e-09 | 9.5818e-06 | 1.9747e-08 | 6.3993e-05 | Yes | Both |
| TNBC - Luminal B | 204369_at | PIK3CA | 0.63381 | 0.63381 | -3.8081e-13 | 0.0015381 | 0.0058595 | 0.0052843 | 0.016096 | Yes | Both |
| TNBC - Luminal B | 231854_at | PIK3CA | 0.03776 | 0.03776 | 4.7622e-14 | 0.51707 | 0.5492 | 0.66231 | 0.66253 | Yes | Neither |
| TNBC - Luminal B | 235980_at | PIK3CA | 0.00075169 | 0.00075169 | 6.0553e-16 | 0.99512 | 0.99554 | 0.99791 | 0.99554 | Yes | Neither |
| TNBC - Luminal B | 205225_at | ESR1 | -7.9556 | -7.9556 | -8.118e-13 | 9.4608e-39 | 3.87e-29 | 5.3926e-37 | 2.2059e-27 | Yes | Both |
| TNBC - Luminal B | 211233_x_at | ESR1 | -0.21573 | -0.21573 | -2.9232e-13 | 0.043351 | 0.030126 | 0.088251 | 0.0676 | Yes | Neither |
| TNBC - Luminal B | 211234_x_at | ESR1 | -0.29312 | -0.29312 | 2.4414e-13 | 0.0035709 | 0.0023024 | 0.010713 | 0.0079537 | Yes | Both |
| TNBC - Luminal B | 211235_s_at | ESR1 | -0.74444 | -0.74444 | -3.898e-13 | 6.7476e-05 | 0.00017722 | 0.00032733 | 0.00089792 | Yes | Both |
| TNBC - Luminal B | 211627_x_at | ESR1 | -0.019226 | -0.019226 | 1.3829e-14 | 0.68737 | 0.61173 | 0.81625 | 0.70442 | Yes | Neither |
| TNBC - Luminal B | 215551_at | ESR1 | -0.32152 | -0.32152 | -2.9249e-13 | 0.0015529 | 0.0027969 | 0.0052843 | 0.009071 | Yes | Both |
| TNBC - Luminal B | 215552_s_at | ESR1 | -1.9061 | -1.9061 | 1.7077e-12 | 5.0425e-11 | 1.5667e-07 | 4.9987e-10 | 1.4884e-06 | Yes | Both |
| TNBC - Luminal B | 217163_at | ESR1 | 0.0024029 | 0.0024029 | 5.2692e-16 | 0.40289 | 0.33147 | 0.56703 | 0.45255 | Yes | Neither |
| TNBC - Luminal B | 217190_x_at | ESR1 | -0.10629 | -0.10629 | -2.7611e-13 | 0.096779 | 0.018338 | 0.17512 | 0.044481 | Yes | Welch only |
| TNBC - Luminal B | 204053_x_at | PTEN | -0.60665 | -0.60665 | 4.3365e-13 | 0.0012688 | 0.005081 | 0.0045273 | 0.014302 | Yes | Both |
| TNBC - Luminal B | 204054_at | PTEN | -0.28698 | -0.28698 | 3.1047e-13 | 0.10429 | 0.13866 | 0.18432 | 0.22744 | Yes | Neither |
| TNBC - Luminal B | 211711_s_at | PTEN | -0.53783 | -0.53783 | 1.3145e-13 | 0.014287 | 0.030322 | 0.035797 | 0.0676 | Yes | Limma only |
| TNBC - Luminal B | 217492_s_at | PTENP1 /// PTEN | -0.54638 | -0.54638 | 3.6171e-13 | 0.0032785 | 0.010542 | 0.01024 | 0.027628 | Yes | Both |
| TNBC - Luminal B | 225363_at | PTEN | -1.0192 | -1.0192 | 1.5457e-12 | 1.6895e-11 | 4.5815e-09 | 1.7509e-10 | 4.9742e-08 | Yes | Both |
| TNBC - Luminal B | 228006_at | PTEN | -1.3833 | -1.3833 | 2.1982e-12 | 1.3925e-09 | 1.15e-07 | 1.1339e-08 | 1.14e-06 | Yes | Both |
| TNBC - Luminal B | 233254_x_at | PTEN | -0.64314 | -0.64314 | 1.2301e-13 | 4.6726e-06 | 1.5295e-05 | 2.8036e-05 | 9.687e-05 | Yes | Both |
| TNBC - Luminal B | 242622_x_at | PTEN | -0.20935 | -0.20935 | 1.4647e-13 | 0.014492 | 0.013127 | 0.035916 | 0.033628 | Yes | Both |
| TNBC - Luminal B | 1565483_at | EGFR | -0.021677 | -0.021677 | 3.0798e-14 | 0.61471 | 0.016149 | 0.7617 | 0.040461 | Yes | Welch only |
| TNBC - Luminal B | 1565484_x_at | EGFR | 0.0001134 | 0.0001134 | 2.0173e-16 | 0.99477 | 0.99185 | 0.99791 | 0.99554 | Yes | Neither |
| TNBC - Luminal B | 201983_s_at | EGFR |  3.423 |  3.423 | 4.9249e-13 | 5.338e-14 | 1.0886e-11 | 7.6067e-13 | 1.5513e-10 | Yes | Both |
| TNBC - Luminal B | 201984_s_at | EGFR | 2.2629 | 2.2629 | 2.3577e-12 | 7.0345e-12 | 4.0693e-13 | 7.6375e-11 | 8.4346e-12 | Yes | Both |
| TNBC - Luminal B | 210984_x_at | EGFR | 0.06482 | 0.06482 | -2.9643e-14 | 0.65534 | 0.38171 | 0.79003 | 0.50599 | Yes | Neither |
| TNBC - Luminal B | 211550_at | EGFR | 0.0123 | 0.0123 | -2.3362e-14 | 0.47679 | 0.49625 | 0.62703 | 0.60831 | Yes | Neither |
| TNBC - Luminal B | 211551_at | EGFR | 0.0031418 | 0.0031418 | 5.1677e-15 | 0.38067 | 0.33077 | 0.54245 | 0.45255 | Yes | Neither |
| TNBC - Luminal B | 211607_x_at | EGFR | 0.053387 | 0.053387 | -4.1744e-14 | 0.69439 | 0.50155 | 0.82032 | 0.61152 | Yes | Neither |
| TNBC - Luminal B | 224999_at | EGFR | 2.3262 | 2.3262 | 2.3599e-12 | 1.0527e-10 | 2.5461e-11 | 1.0001e-09 | 3.4147e-10 | Yes | Both |
| Her2 - Luminal A | 204531_s_at | BRCA1 | 0.6931 | 0.6931 | -1.1802e-13 | 0.0045851 | 0.00062035 | 0.012906 | 0.0027733 | Yes | Both |
| Her2 - Luminal A | 211851_x_at | BRCA1 | 0.085888 | 0.085888 | -1.7486e-15 | 0.46387 | 0.092779 | 0.61849 | 0.16026 | Yes | Neither |
| Her2 - Luminal A | 208368_s_at | BRCA2 |  1.005 |  1.005 | -4.1411e-12 | 2.1116e-06 | 1.7128e-06 | 1.3012e-05 | 1.2204e-05 | Yes | Both |
| Her2 - Luminal A | 214727_at | BRCA2 | 0.40492 | 0.40492 | 3.525e-14 | 0.011734 | 0.0036895 | 0.030061 | 0.011216 | Yes | Both |
| Her2 - Luminal A | 201746_at | TP53 | -0.47616 | -0.47616 | 3.9285e-13 | 0.1153 | 0.097767 | 0.20068 | 0.1676 | Yes | Neither |
| Her2 - Luminal A | 211300_s_at | TP53 | -0.099953 | -0.099953 | 3.5805e-15 | 0.75275 | 0.72215 | 0.85351 | 0.79159 | Yes | Neither |
| Her2 - Luminal A | 210930_s_at | ERBB2 | 5.9129 | 5.9129 | -4.3663e-12 | 1.2646e-34 | 1.5428e-28 | 4.8054e-33 | 7.0353e-27 | Yes | Both |
| Her2 - Luminal A | 216836_s_at | ERBB2 | 4.1477 | 4.1477 | 3.6717e-12 | 3.7318e-29 | 3.4866e-28 | 9.4539e-28 | 1.3249e-26 | Yes | Both |
| Her2 - Luminal A | 234354_x_at | ERBB2 | 3.5406 | 3.5406 | -2.6459e-12 | 2.1309e-20 | 6.7015e-13 | 4.0486e-19 | 1.2733e-11 | Yes | Both |
| Her2 - Luminal A | 204369_at | PIK3CA | -0.25867 | -0.25867 | -2.5374e-13 | 0.22507 | 0.18297 | 0.35391 | 0.2877 | Yes | Neither |
| Her2 - Luminal A | 231854_at | PIK3CA | -0.076955 | -0.076955 | -4.3715e-14 | 0.22413 | 0.06039 | 0.35391 | 0.11474 | Yes | Neither |
| Her2 - Luminal A | 235980_at | PIK3CA | -0.35492 | -0.35492 | 2.0495e-13 | 0.0086551 | 0.0012397 | 0.023492 | 0.004559 | Yes | Both |
| Her2 - Luminal A | 205225_at | ESR1 | -7.8638 | -7.8638 | -4.3672e-12 | 4.3787e-35 | 1.2621e-24 | 1.9967e-33 | 4.111e-23 | Yes | Both |
| Her2 - Luminal A | 211233_x_at | ESR1 | -0.4528 | -0.4528 | 6.9e-14 | 0.00012813 | 0.0038237 | 0.00055119 | 0.011322 | Yes | Both |
| Her2 - Luminal A | 211234_x_at | ESR1 | -0.46713 | -0.46713 | -9.1427e-14 | 2.6052e-05 | 0.0016903 | 0.00014142 | 0.005929 | Yes | Both |
| Her2 - Luminal A | 211235_s_at | ESR1 | -1.1849 | -1.1849 | -3.4501e-12 | 1.4968e-08 | 2.7002e-05 | 1.1009e-07 | 0.00016201 | Yes | Both |
| Her2 - Luminal A | 211627_x_at | ESR1 | -0.053864 | -0.053864 | -2.9594e-14 | 0.29914 | 0.41363 | 0.44578 | 0.53119 | Yes | Neither |
| Her2 - Luminal A | 215551_at | ESR1 | -0.24423 | -0.24423 | -8.1962e-14 | 0.025087 | 0.065909 | 0.056076 | 0.12217 | Yes | Neither |
| Her2 - Luminal A | 215552_s_at | ESR1 | -2.2285 | -2.2285 | 1.8381e-12 | 2.559e-12 | 2.089e-07 | 2.9172e-11 | 1.701e-06 | Yes | Both |
| Her2 - Luminal A | 217163_at | ESR1 | 0.0057363 | 0.0057363 | 4.7028e-15 | 0.066863 | 0.038898 | 0.12704 | 0.079186 | Yes | Neither |
| Her2 - Luminal A | 217190_x_at | ESR1 | -0.20668 | -0.20668 | -1.9629e-13 | 0.0032377 | 0.035825 | 0.01024 | 0.074256 | Yes | Limma only |
| Her2 - Luminal A | 204053_x_at | PTEN | 0.29161 | 0.29161 | -4.8506e-13 | 0.14618 | 0.0082732 | 0.24689 | 0.021934 | Yes | Welch only |
| Her2 - Luminal A | 204054_at | PTEN | -0.073704 | -0.073704 | 4.0482e-14 | 0.6989 | 0.56795 | 0.82139 | 0.67798 | Yes | Neither |
| Her2 - Luminal A | 211711_s_at | PTEN | 0.01219 | 0.01219 | 3.4313e-14 | 0.95865 | 0.93634 | 0.99791 | 0.97929 | Yes | Neither |
| Her2 - Luminal A | 217492_s_at | PTENP1 /// PTEN | 0.22587 | 0.22587 | 4.3715e-13 | 0.25513 | 0.048651 | 0.38779 | 0.094807 | Yes | Neither |
| Her2 - Luminal A | 225363_at | PTEN | -0.58409 | -0.58409 | 3.1752e-14 | 0.00014917 | 7.0689e-07 | 0.00062981 | 5.3723e-06 | Yes | Both |
| Her2 - Luminal A | 228006_at | PTEN | -0.91124 | -0.91124 | -2.5047e-13 | 0.00011846 | 2.8321e-07 | 0.00053363 | 2.2266e-06 | Yes | Both |
| Her2 - Luminal A | 233254_x_at | PTEN | -0.21972 | -0.21972 | -3.3196e-13 | 0.13397 | 0.075989 | 0.22795 | 0.13536 | Yes | Neither |
| Her2 - Luminal A | 242622_x_at | PTEN | 0.041552 | 0.041552 | -2.5181e-14 | 0.65065 | 0.62366 | 0.78908 | 0.71455 | Yes | Neither |
| Her2 - Luminal A | 1565483_at | EGFR | -0.05412 | -0.05412 | 3.9566e-14 | 0.24728 | 0.44185 | 0.38354 | 0.55353 | Yes | Neither |
| Her2 - Luminal A | 1565484_x_at | EGFR | 0.0054561 | 0.0054561 | 3.4885e-15 | 0.77115 | 0.82004 | 0.85351 | 0.87369 | Yes | Neither |
| Her2 - Luminal A | 201983_s_at | EGFR | 0.81435 | 0.81435 | 3.7736e-13 | 0.065422 | 0.10936 | 0.12535 | 0.18607 | Yes | Neither |
| Her2 - Luminal A | 201984_s_at | EGFR | 1.3879 | 1.3879 | -3.4446e-12 | 3.6665e-05 | 0.00044841 | 0.00019441 | 0.0020865 | Yes | Both |
| Her2 - Luminal A | 210984_x_at | EGFR | 0.21104 | 0.21104 | -3.4511e-13 | 0.18136 | 0.31962 | 0.30028 | 0.44707 | Yes | Neither |
| Her2 - Luminal A | 211550_at | EGFR | 0.038805 | 0.038805 | 2.983e-14 | 0.039818 | 0.054643 | 0.081787 | 0.10558 | Yes | Neither |
| Her2 - Luminal A | 211551_at | EGFR | 0.011563 | 0.011563 | -1.5614e-14 | 0.0033831 | 0.0034863 | 0.010423 | 0.010742 | Yes | Both |
| Her2 - Luminal A | 211607_x_at | EGFR | 0.18398 | 0.18398 | 1.0242e-13 | 0.21286 | 0.33473 | 0.34419 | 0.45427 | Yes | Neither |
| Her2 - Luminal A | 224999_at | EGFR | 0.44662 | 0.44662 | 3.0026e-13 | 0.21436 | 0.27529 | 0.34419 | 0.4077 | Yes | Neither |
| Her2 - Luminal B | 204531_s_at | BRCA1 | -0.015483 | -0.015483 | 1.7625e-15 | 0.94825 | 0.93236 | 0.99791 | 0.97929 | Yes | Neither |
| Her2 - Luminal B | 211851_x_at | BRCA1 | -0.0041327 | -0.0041327 | -4.8893e-15 | 0.97161 | 0.96509 | 0.99791 | 0.99117 | Yes | Neither |
| Her2 - Luminal B | 208368_s_at | BRCA2 | 0.40669 | 0.40669 | -1.7208e-15 | 0.04456 | 0.0581 | 0.089119 | 0.11132 | Yes | Neither |
| Her2 - Luminal B | 214727_at | BRCA2 | 0.24279 | 0.24279 | -3.3304e-13 | 0.1245 | 0.11598 | 0.21343 | 0.19587 | Yes | Neither |
| Her2 - Luminal B | 201746_at | TP53 | -0.093911 | -0.093911 | 2.1094e-15 | 0.75296 | 0.75394 | 0.85351 | 0.81856 | Yes | Neither |
| Her2 - Luminal B | 211300_s_at | TP53 | 0.42517 | 0.42517 | -6.6613e-16 | 0.17801 | 0.11924 | 0.29843 | 0.19991 | Yes | Neither |
| Her2 - Luminal B | 210930_s_at | ERBB2 | 2.8457 | 2.8457 | 1.3323e-15 | 1.4765e-13 | 1.9067e-07 | 1.8702e-12 | 1.6101e-06 | Yes | Both |
| Her2 - Luminal B | 216836_s_at | ERBB2 | 0.92238 | 0.92238 | -3.3618e-13 | 0.0012708 | 0.0029033 | 0.0045273 | 0.0091937 | Yes | Both |
| Her2 - Luminal B | 234354_x_at | ERBB2 | 1.6519 | 1.6519 | 3.3316e-12 | 7.1175e-07 | 0.00072068 | 4.6365e-06 | 0.0030499 | Yes | Both |
| Her2 - Luminal B | 204369_at | PIK3CA | 0.017039 | 0.017039 | 2.3037e-15 | 0.93557 | 0.94413 | 0.99791 | 0.98293 | Yes | Neither |
| Her2 - Luminal B | 231854_at | PIK3CA | -0.13189 | -0.13189 | 3.3415e-13 | 0.036682 | 0.004209 | 0.076729 | 0.012022 | Yes | Welch only |
| Her2 - Luminal B | 235980_at | PIK3CA | -0.31449 | -0.31449 | 3.3007e-13 | 0.018641 | 0.01333 | 0.043368 | 0.03377 | Yes | Both |
| Her2 - Luminal B | 205225_at | ESR1 | -7.6088 | -7.6088 | 1.7764e-15 | 3.8067e-34 | 4.8413e-24 | 1.2399e-32 | 1.3798e-22 | Yes | Both |
| Her2 - Luminal B | 211233_x_at | ESR1 | -0.21474 | -0.21474 | -5.5511e-17 | 0.061014 | 0.030835 | 0.11789 | 0.0676 | Yes | Neither |
| Her2 - Luminal B | 211234_x_at | ESR1 | -0.29066 | -0.29066 | -3.3301e-13 | 0.007035 | 0.0024744 | 0.019561 | 0.0081763 | Yes | Both |
| Her2 - Luminal B | 211235_s_at | ESR1 | -0.74068 | -0.74068 | -4.4409e-16 | 0.00021103 | 0.0001889 | 0.00085921 | 0.00093628 | Yes | Both |
| Her2 - Luminal B | 211627_x_at | ESR1 | -0.017789 | -0.017789 | -3.4132e-14 | 0.72897 | 0.63889 | 0.84799 | 0.72834 | Yes | Neither |
| Her2 - Luminal B | 215551_at | ESR1 | -0.32754 | -0.32754 | 3.3246e-13 | 0.0026488 | 0.0023824 | 0.0087526 | 0.0080016 | Yes | Both |
| Her2 - Luminal B | 215552_s_at | ESR1 | -1.8898 | -1.8898 | -3.3333e-12 | 8.7036e-10 | 1.825e-07 | 7.6324e-09 | 1.6101e-06 | Yes | Both |
| Her2 - Luminal B | 217163_at | ESR1 | 0.0032401 | 0.0032401 | -1.0278e-16 | 0.29431 | 0.21135 | 0.44146 | 0.32724 | Yes | Neither |
| Her2 - Luminal B | 217190_x_at | ESR1 | -0.067538 | -0.067538 | -8.743e-16 | 0.32444 | 0.18061 | 0.47117 | 0.28597 | Yes | Neither |
| Her2 - Luminal B | 204053_x_at | PTEN | 0.13354 | 0.13354 | 3.3581e-13 | 0.50073 | 0.28155 |  0.645 | 0.40888 | Yes | Neither |
| Her2 - Luminal B | 204054_at | PTEN | 0.060242 | 0.060242 | 3.357e-14 | 0.7498 | 0.66845 | 0.85351 | 0.75449 | Yes | Neither |
| Her2 - Luminal B | 211711_s_at | PTEN | 0.18345 | 0.18345 | -2.7756e-17 | 0.43189 | 0.28705 | 0.58965 | 0.41423 | Yes | Neither |
| Her2 - Luminal B | 217492_s_at | PTENP1 /// PTEN | 0.13925 | 0.13925 | -3.3273e-13 | 0.47853 | 0.29439 | 0.62703 | 0.42215 | Yes | Neither |
| Her2 - Luminal B | 225363_at | PTEN | -0.31733 | -0.31733 | -3.329e-13 | 0.034034 | 0.0040821 | 0.073205 | 0.011932 | Yes | Welch only |
| Her2 - Luminal B | 228006_at | PTEN | -0.5431 | -0.5431 | 3.3462e-13 | 0.018411 | 0.00086768 | 0.043274 | 0.0035969 | Yes | Both |
| Her2 - Luminal B | 233254_x_at | PTEN | -0.086984 | -0.086984 | -3.0725e-14 | 0.5481 | 0.42962 | 0.69426 | 0.54419 | Yes | Neither |
| Her2 - Luminal B | 242622_x_at | PTEN | 0.016867 | 0.016867 | 3.3227e-14 | 0.85286 | 0.83326 | 0.93039 | 0.87955 | Yes | Neither |
| Her2 - Luminal B | 1565483_at | EGFR | -0.00049863 | -0.00049863 | 5.7625e-16 | 0.9914 | 0.98146 | 0.99791 | 0.99554 | Yes | Neither |
| Her2 - Luminal B | 1565484_x_at | EGFR | 0.0052006 | 0.0052006 | -4.0601e-15 | 0.77976 | 0.81978 | 0.85887 | 0.87369 | Yes | Neither |
| Her2 - Luminal B | 201983_s_at | EGFR | 2.3612 | 2.3612 | -1.3323e-15 | 2.7066e-07 | 3.2766e-05 | 1.815e-06 | 0.00018677 | Yes | Both |
| Her2 - Luminal B | 201984_s_at | EGFR | 1.5307 | 1.5307 | 3.336e-12 | 5.1998e-06 | 0.00014275 | 3.0399e-05 | 0.00073972 | Yes | Both |
| Her2 - Luminal B | 210984_x_at | EGFR | 0.20906 | 0.20906 | -6.6613e-16 | 0.18175 | 0.3244 | 0.30028 | 0.45028 | Yes | Neither |
| Her2 - Luminal B | 211550_at | EGFR | 0.0062119 | 0.0062119 | -2.9473e-15 | 0.73792 | 0.71587 | 0.85351 | 0.78849 | Yes | Neither |
| Her2 - Luminal B | 211551_at | EGFR | 0.0036801 | 0.0036801 | -4.9136e-16 | 0.33945 | 0.23548 | 0.48984 | 0.35793 | Yes | Neither |
| Her2 - Luminal B | 211607_x_at | EGFR | 0.16589 | 0.16589 | 3.3368e-13 |  0.257 | 0.38456 | 0.38806 | 0.50682 | Yes | Neither |
| Her2 - Luminal B | 224999_at | EGFR | 1.5129 | 1.5129 | -3.3333e-12 | 3.9008e-05 | 0.00049747 | 0.00020213 | 0.0022685 | Yes | Both |
| Luminal A - Luminal B | 204531_s_at | BRCA1 | -0.70858 | -0.70858 | 1.1979e-13 | 0.0037782 | 0.00029064 | 0.011044 | 0.0014099 | Yes | Both |
| Luminal A - Luminal B | 211851_x_at | BRCA1 | -0.090021 | -0.090021 | -3.314e-14 | 0.44269 | 0.30651 | 0.59724 | 0.43406 | Yes | Neither |
| Luminal A - Luminal B | 208368_s_at | BRCA2 | -0.59835 | -0.59835 | 1.3933e-13 | 0.0036758 | 0.00072235 | 0.010884 | 0.0030499 | Yes | Both |
| Luminal A - Luminal B | 214727_at | BRCA2 | -0.16213 | -0.16213 | -3.6829e-13 | 0.30786 | 0.0712 | 0.45579 | 0.12987 | Yes | Neither |
| Luminal A - Luminal B | 201746_at | TP53 | 0.38225 | 0.38225 | -3.9074e-13 | 0.20536 | 0.048108 | 0.33445 | 0.094807 | Yes | Neither |
| Luminal A - Luminal B | 211300_s_at | TP53 | 0.52513 | 0.52513 | -4.0423e-13 | 0.099672 | 0.017773 | 0.17894 | 0.043574 | Yes | Welch only |
| Luminal A - Luminal B | 210930_s_at | ERBB2 | -3.0672 | -3.0672 | 4.3681e-12 | 6.4897e-15 | 3.014e-08 | 9.8644e-14 | 3.1236e-07 | Yes | Both |
| Luminal A - Luminal B | 216836_s_at | ERBB2 | -3.2254 | -3.2254 | 2.9918e-12 | 3.1759e-21 | 1.0436e-13 | 6.5828e-20 | 2.3794e-12 | Yes | Both |
| Luminal A - Luminal B | 234354_x_at | ERBB2 | -1.8888 | -1.8888 | -4.0226e-12 | 2.8619e-08 | 9.8235e-06 | 2.0391e-07 | 6.3993e-05 | Yes | Both |
| Luminal A - Luminal B | 204369_at | PIK3CA | 0.27571 | 0.27571 | 2.5607e-13 | 0.19617 | 0.21242 | 0.32177 | 0.32724 | Yes | Neither |
| Luminal A - Luminal B | 231854_at | PIK3CA | -0.054936 | -0.054936 | -2.2142e-14 | 0.38482 | 0.035402 | 0.54497 | 0.074052 | Yes | Neither |
| Luminal A - Luminal B | 235980_at | PIK3CA | 0.040432 | 0.040432 | 2.5153e-14 | 0.76177 | 0.68349 | 0.85351 | 0.7639 | Yes | Neither |
| Luminal A - Luminal B | 205225_at | ESR1 | 0.25503 | 0.25503 | 3.6932e-13 | 0.57653 | 0.48214 | 0.72624 | 0.59744 | Yes | Neither |
| Luminal A - Luminal B | 211233_x_at | ESR1 | 0.23806 | 0.23806 | -6.9028e-14 | 0.039761 | 0.17246 | 0.081787 | 0.27534 | Yes | Neither |
| Luminal A - Luminal B | 211234_x_at | ESR1 | 0.17647 | 0.17647 | -2.4161e-13 | 0.10158 | 0.27717 | 0.18094 | 0.4077 | Yes | Neither |
| Luminal A - Luminal B | 211235_s_at | ESR1 | 0.44427 | 0.44427 | 4.4964e-13 | 0.024935 | 0.1357 | 0.056076 | 0.2242 | Yes | Neither |
| Luminal A - Luminal B | 211627_x_at | ESR1 | 0.036075 | 0.036075 | -4.538e-15 | 0.48631 | 0.60878 | 0.62999 | 0.70442 | Yes | Neither |
| Luminal A - Luminal B | 215551_at | ESR1 | -0.083315 | -0.083315 | 1.4461e-14 | 0.44076 | 0.60712 | 0.59724 | 0.70442 | Yes | Neither |
| Luminal A - Luminal B | 215552_s_at | ESR1 | 0.33874 | 0.33874 | -1.7136e-13 | 0.24085 | 0.43331 | 0.37612 | 0.54583 | Yes | Neither |
| Luminal A - Luminal B | 217163_at | ESR1 | -0.0024961 | -0.0024961 | -4.8056e-15 | 0.42269 | 0.4147 | 0.58408 | 0.53119 | Yes | Neither |
| Luminal A - Luminal B | 217190_x_at | ESR1 | 0.13914 | 0.13914 | 1.954e-13 | 0.04543 | 0.17269 | 0.09007 | 0.27534 | Yes | Neither |
| Luminal A - Luminal B | 204053_x_at | PTEN | -0.15808 | -0.15808 | -1.7911e-13 | 0.4295 | 0.16404 | 0.58965 | 0.26526 | Yes | Neither |
| Luminal A - Luminal B | 204054_at | PTEN | 0.13395 | 0.13395 | 9.3092e-14 | 0.48238 | 0.25575 | 0.62847 | 0.38617 | Yes | Neither |
| Luminal A - Luminal B | 211711_s_at | PTEN | 0.17126 | 0.17126 | -3.3434e-13 | 0.46685 | 0.3147 | 0.61885 | 0.44291 | Yes | Neither |
| Luminal A - Luminal B | 217492_s_at | PTENP1 /// PTEN | -0.086628 | -0.086628 | 3.0115e-14 | 0.66182 | 0.49603 | 0.79003 | 0.60831 | Yes | Neither |
| Luminal A - Luminal B | 225363_at | PTEN | 0.26676 | 0.26676 | -3.6465e-13 | 0.076476 | 0.0056034 | 0.14292 | 0.01558 | Yes | Welch only |
| Luminal A - Luminal B | 228006_at | PTEN | 0.36814 | 0.36814 | -4.1489e-13 | 0.11098 | 0.0068807 | 0.19464 | 0.018456 | Yes | Welch only |
| Luminal A - Luminal B | 233254_x_at | PTEN | 0.13274 | 0.13274 | 1.2212e-15 | 0.36393 | 0.22984 | 0.52186 | 0.3517 | Yes | Neither |
| Luminal A - Luminal B | 242622_x_at | PTEN | -0.024685 | -0.024685 | -4.1592e-14 | 0.78785 | 0.71055 | 0.8636 | 0.78643 | Yes | Neither |
| Luminal A - Luminal B | 1565483_at | EGFR | 0.053621 | 0.053621 | -3.899e-14 | 0.25164 | 0.42866 | 0.38684 | 0.54419 | Yes | Neither |
| Luminal A - Luminal B | 1565484_x_at | EGFR | -0.00025554 | -0.00025554 | 4.507e-16 | 0.98913 | 0.9838 | 0.99791 | 0.99554 | Yes | Neither |
| Luminal A - Luminal B | 201983_s_at | EGFR | 1.5468 | 1.5468 | -1.3789e-12 | 0.00057966 | 0.0010506 | 0.0022027 | 0.0041298 | Yes | Both |
| Luminal A - Luminal B | 201984_s_at | EGFR | 0.14274 | 0.14274 | -2.1941e-13 | 0.66061 | 0.30582 | 0.79003 | 0.43406 | Yes | Neither |
| Luminal A - Luminal B | 210984_x_at | EGFR | -0.0019809 | -0.0019809 | 4.4301e-15 | 0.98996 | 0.95864 | 0.99791 |  0.989 | Yes | Neither |
| Luminal A - Luminal B | 211550_at | EGFR | -0.032594 | -0.032594 | 3.722e-14 | 0.083487 | 0.18548 | 0.15476 | 0.28966 | Yes | Neither |
| Luminal A - Luminal B | 211551_at | EGFR | -0.0078825 | -0.0078825 | -4.8763e-15 | 0.043802 | 0.073555 | 0.088379 | 0.13205 | Yes | Neither |
| Luminal A - Luminal B | 211607_x_at | EGFR | -0.018086 | -0.018086 | 3.1256e-14 | 0.90224 | 0.58003 | 0.97493 | 0.68168 | Yes | Neither |
| Luminal A - Luminal B | 224999_at | EGFR | 1.0662 | 1.0662 | 4.3665e-12 | 0.0034651 | 0.00097166 | 0.010534 | 0.0038867 | Yes | Both |
