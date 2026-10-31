# OncoMap v2 — Breast Cancer Expression and Molecular Subtype Explorer

**Final exploratory research report**<br>
**Data source:** GEO series GSE42568 and GSE65194<br>
**Scope:** Probe-level computational analysis of a prespecified eight-gene, 38-probe panel

## 1. Abstract

This exploratory analysis compares probe-level expression for eight breast-cancer-related target genes (38 GPL570 probes) across GSE42568 and GSE65194, and examines tumour–healthy and subtype differences in GSE65194. The source cohort results include 104 cancer and 17 normal samples for GSE42568 and 130 tumour patients plus 11 healthy controls in GSE65194. Across exact-probe cohort matches, 24/38 directions were concordant and 14/38 discordant. In the matched-38 BH comparison, 15 probes were significant in both cohorts. Patient-level subtype comparisons showed 228/228 directional agreement between Welch and limma; 100 and 95 tests, respectively, met adjusted p < 0.05. These results are computational and exploratory, not clinical or diagnostic validation.

## 2. Introduction

Breast tumours show heterogeneous transcriptional profiles. Comparing expression between tumour and non-tumour breast tissue, and across clinically annotated molecular subtype groups, can help describe which probe-level patterns are present in public expression datasets. This project focuses on eight target genes—BRCA1, BRCA2, TP53, ERBB2, PIK3CA, ESR1, PTEN and EGFR—selected in the existing OncoMap target panel. Multiple microarray probes can represent a gene and may behave differently; therefore, all results remain at the exact-probe level rather than averaging probes into a synthetic gene-expression value.

An independent-cohort comparison is useful because a pattern that recurs in another dataset is less likely to depend entirely on one cohort's sample composition or analysis choices. It is not, by itself, biological validation: preprocessing, recruitment, sample sizes and statistical methods differ between cohorts, and the panel is targeted rather than genome-wide.

## 3. Research objectives

1. Describe how the 38 target probes differ between tumour and healthy breast tissue in GSE65194.
2. Describe probe-level differences across TNBC, Her2, Luminal A and Luminal B patient groups.
3. Assess whether tumour-versus-normal expression directions are consistent for exact shared probes in the independent GSE42568 cohort.

## 4. Dataset description

### GSE42568

The existing GSE42568 result source is a probe-level GEO2R/limma table from the GPL570 platform. The documented comparison is cancer tissue minus normal breast tissue, with **104 cancer and 17 normal samples**. Positive logFC indicates higher expression in cancer. Its source GEO2R adjusted p-values come from a broader probe-testing family. This report uses the existing result and does not reconstruct the cohort analysis.

### GSE65194

The downloaded GEO series matrix contains **54,673 probes × 178 GSM arrays** on GPL570. Array composition in the audit is 153 tumour arrays representing 130 patients, 11 healthy breast controls and 14 cell-line arrays. The analysed patient-level cohort has **130 tumour patients** and 11 controls; subtype counts are TNBC 41, Her2 30, Luminal A 29, Luminal B 30. The cell-line arrays are excluded from primary tumour/healthy and patient-level subtype comparisons.

Both series identify GPL570 and report GC-RMA-related processing. GSE65194 also describes fixed-effect correction for array and hybridization batch; equivalent batch correction is not documented for GSE42568. GSE65194's exact CDF choice (standard AffyCDF or Brainarray) and explicit transformation history remain unresolved. The project analyses existing processed values as log2-scale expression and does not repeat transformation or normalization. The GEO records are available at [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568) and [GSE65194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65194).

## 5. Data preparation and quality control

The source-array matrix and manifest were checked for unique probe/GSM identifiers, numeric expression values and exact metadata matching. Preliminary QC found 54,673 unique probe rows and 178 unique GSM columns, no missing or non-numeric measurements, and no sample-median screen that was used to remove arrays. The screen was descriptive only; flagged arrays were retained.

The source manifest identifies 153 tumour arrays for 130 tumour patients: **23 technical replicate pairs** and 107 singleton arrays. Pair-level correlations across 54,673 probes were high, but correlation does not establish sample provenance or prove prior averaging. The deposited series matrix retains separate GSM columns for both members of all 23 pairs even though the GEO processing annotation states, “Samples with technical replicates were subsequently averaged.” The available records do not resolve where or whether that documented averaging occurred. Accordingly, the patient-level matrix is explicitly a **project-derived dataset**: for each validated pair, the two deposited processed values were arithmetically averaged probe-by-probe; single-array patients and healthy controls were retained unchanged. The original 178-array matrix was not rewritten. Replicate-source uncertainty remains a limitation.

The derived matrix has 54,673 probes × 141 observations. Patient-level tumour/healthy and subtype analyses use one tumour observation per patient; the 11 controls are used only for tumour-versus-healthy comparisons. The 14 cell-line arrays are not included in the derived matrix.

## 6. Statistical methodology

### Welch tumour-versus-healthy analysis

The existing patient-level Welch results use two-sided independent-samples Welch t-tests, without assuming equal group variances. Effect is the mean tumour minus healthy log2 expression; positive values indicate higher tumour expression. Benjamini–Hochberg (BH) adjustment was applied over the **38 target-probe tests**.

### limma tumour-versus-healthy sensitivity

Two local limma fits select one array per paired patient using alternative assignments: repA or repB, while retaining the same 107 singleton tumour samples and 11 healthy controls. These are **sensitivity analyses over alternative technical-replicate selections, not independent cohorts**. Both use the contrast Tumour − Healthy; positive logFC means higher tumour expression. Full limma tables retain all **54,673 probes**, and BH adjustment for each fit spans that full set. Hence adjusted-p significance for the 38 targets need not match Welch's 38-test family. The local fits use the downloaded matrix and planned sample assignments; they are not claimed to reproduce the GEO2R website exactly.

### Patient-level subtype analyses

The subtype analyses use the existing averaged patient-level matrix and 130 tumour patients only. For each of 6 contrasts—TNBC − Her2, TNBC − Luminal A, TNBC − Luminal B, Her2 − Luminal A, Her2 − Luminal B and Luminal A − Luminal B—Welch tests and empirical-Bayes limma estimates are available for the same 38 probes. Positive effect means higher expression in the first-named subtype. BH adjustment is across the complete **228-test family (38 probes × 6 contrasts)** for each method, not separately within each contrast. Limma uses a four-group design, `lmFit`, `makeContrasts`, `contrasts.fit` and `eBayes`.

### Independent-cohort matching and interpretation of effects

Cross-cohort comparisons use exact GPL570 probe IDs; gene labels are descriptive, not the join key. GSE42568 logFC and GSE65194 tumour-minus-healthy mean differences have the same stated orientation. Original adjusted p-values reflect the original method-specific families. A separate matched-panel file recalculates BH within the same **38 exact probes separately in each cohort**, which permits a like-for-like panel-size correction but does not harmonize preprocessing or study design. No cohort matrices were pooled.

## 7. Results

### Independent-cohort replication

Among the 38 exact shared probes, **24** had concordant and **14** had discordant tumour/cancer-minus-normal/healthy directions. Under the original cohort-specific adjustments, 15 were significant in both, 3 only in GSE42568, 12 only in GSE65194, and 8 in neither. These are descriptive replication counts, not a pooled meta-analysis.

| Probe | Gene | GSE42568 annotation | GSE42568 logFC | GSE42568 original BH p | GSE65194 tumour−healthy difference | GSE65194 BH p | Direction | Significant in both |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 204531_s_at | BRCA1 | BRCA1 | 0.357 | 0.17 | 0.4305 | 0.05263 | concordant | No |
| 211851_x_at | BRCA1 | BRCA1 | 0.003 | 0.979 | 0.1564 | 0.004399 | concordant | No |
| 208368_s_at | BRCA2 | BRCA2 | 0.536 | 0.00151 | 0.7524 | 0.001212 | concordant | Yes |
| 214727_at | BRCA2 | BRCA2 | 0.019 | 0.654 | 0.2055 | 0.0122 | concordant | No |
| 201746_at | TP53 | TP53 | 0.671 | 0.0167 | -0.8253 | 0.0005365 | discordant | Yes |
| 211300_s_at | TP53 | TP53 | 0.037 | 0.688 | 1.2248 | 1.86e-07 | concordant | No |
| 210930_s_at | ERBB2 | ERBB2 | 0.068 | 0.652 | 2.0906 | 4.32e-09 | concordant | No |
| 216836_s_at | ERBB2 | ERBB2 | 3.120 | 3.52e-08 | 1.8863 | 2.04e-06 | concordant | Yes |
| 234354_x_at | ERBB2 | ERBB2 | 0.499 | 0.162 | 1.2108 | 2.76e-07 | concordant | No |
| 204369_at | PIK3CA | PIK3CA | -1.470 | 1.89e-08 | 2.9393 | 0.0009249 | discordant | Yes |
| 231854_at | PIK3CA | PIK3CA | -1.290 | 3.70e-16 | -0.3928 | 0.01963 | concordant | Yes |
| 235980_at | PIK3CA | PIK3CA | -1.080 | 9.62e-09 | -1.3486 | 0.0004011 | concordant | Yes |
| 205225_at | ESR1 | ESR1 | 1.750 | 0.105 | -0.9293 | 0.1861 | discordant | No |
| 211233_x_at | ESR1 | ESR1 | 0.056 | 0.648 | 0.1449 | 0.001484 | concordant | No |
| 211234_x_at | ESR1 | ESR1 | 0.061 | 0.429 | 0.1524 | 0.0009249 | concordant | No |
| 211235_s_at | ESR1 | ESR1 | 0.114 | 0.414 | 0.3442 | 0.001212 | concordant | No |
| 211627_x_at | ESR1 | ESR1 | 0.103 | 0.722 | -0.8792 | 0.001034 | discordant | No |
| 215551_at | ESR1 | ESR1 | 0.035 | 0.831 | -0.0253 | 0.7336 | discordant | No |
| 215552_s_at | ESR1 | ESR1 | 0.924 | 0.0419 | 0.6748 | 0.001934 | concordant | Yes |
| 217163_at | ESR1 | ESR1 | 0.001 | 0.966 | -0.0071 | 0.1695 | discordant | No |
| 217190_x_at | ESR1 | ESR1 | 0.170 | 0.18 | 0.0532 | 0.08792 | concordant | No |
| 204053_x_at | PTEN | PTEN | -2.030 | 9.78e-13 | 0.3799 | 0.1861 | discordant | No |
| 204054_at | PTEN | PTEN | -1.120 | 3.73e-08 | -0.3069 | 0.006888 | concordant | Yes |
| 211711_s_at | PTEN | PTEN | -2.370 | 4.41e-15 | -0.9506 | 0.001212 | concordant | Yes |
| 217492_s_at | PTEN | PTENP1///PTEN | -2.020 | 1.66e-13 | 0.7118 | 0.04205 | discordant | Yes |
| 225363_at | PTEN | PTEN | -1.210 | 9.55e-05 | -1.0096 | 3.98e-06 | concordant | Yes |
| 228006_at | PTEN | PTEN | -1.980 | 4.34e-12 | -1.5964 | 4.32e-09 | concordant | Yes |
| 233254_x_at | PTEN | PTEN | 0.241 | 0.375 | 0.0065 | 0.9757 | concordant | No |
| 242622_x_at | PTEN | PTEN | 0.348 | 0.0853 | -0.6158 | 0.004498 | discordant | No |
| 1565483_at | EGFR | EGFR | 0.116 | 0.826 | -0.0381 | 0.4642 | discordant | No |
| 1565484_x_at | EGFR | EGFR | 0.163 | 0.547 | -0.2185 | 0.1729 | discordant | No |
| 201983_s_at | EGFR | EGFR | -2.510 | 3.89e-05 | -1.8977 | 0.001164 | concordant | Yes |
| 201984_s_at | EGFR | EGFR | -2.090 | 9.35e-08 | -0.9937 | 0.01963 | concordant | Yes |
| 210984_x_at | EGFR | EGFR | -0.357 | 0.000592 | 0.0960 | 0.2101 | discordant | No |
| 211550_at | EGFR | EGFR | 0.002 | 0.928 | -0.1512 | 0.001484 | discordant | No |
| 211551_at | EGFR | EGFR | -0.111 | 0.179 | -0.0262 | 0.03735 | concordant | No |
| 211607_x_at | EGFR | EGFR | -0.226 | 0.0445 | 0.0778 | 0.2783 | discordant | No |
| 224999_at | EGFR | EGFR | -2.390 | 7.14e-08 | -1.9580 | 0.0005327 | concordant | Yes |

### Matched-panel BH comparison

Recalculation across exactly 38 probes yielded 18 significant probes in GSE42568 and 27 in GSE65194. Categories were 15 significant in both, 3 only in GSE42568, 12 only in GSE65194 and 8 in neither. The original significance counts were unchanged by this matched-panel recalculation in these outputs. The two datasets still differ in cohort and preprocessing.

| Probe | Gene | GSE42568 matched-38 BH p | GSE65194 matched-38 BH p | GSE42568 significant | GSE65194 significant |
| --- | --- | --- | --- | --- | --- |
| 204531_s_at | BRCA1 | 0.1205 | 0.05263 | No | No |
| 211851_x_at | BRCA1 | 0.956 | 0.004399 | No | Yes |
| 208368_s_at | BRCA2 | 0.0006688 | 0.001212 | Yes | Yes |
| 214727_at | BRCA2 | 0.57 | 0.0122 | No | Yes |
| 201746_at | TP53 | 0.009809 | 0.0005365 | Yes | Yes |
| 211300_s_at | TP53 | 0.6009 | 1.86e-07 | No | Yes |
| 210930_s_at | ERBB2 | 0.57 | 4.32e-09 | No | Yes |
| 216836_s_at | ERBB2 | 9.37e-09 | 2.04e-06 | Yes | Yes |
| 234354_x_at | ERBB2 | 0.12 | 2.76e-07 | No | Yes |
| 204369_at | PIK3CA | 5.75e-09 | 0.0009249 | Yes | Yes |
| 231854_at | PIK3CA | 2.01e-16 | 0.01963 | Yes | Yes |
| 235980_at | PIK3CA | 3.22e-09 | 0.0004011 | Yes | Yes |
| 205225_at | ESR1 | 0.07277 | 0.1861 | No | No |
| 211233_x_at | ESR1 | 0.57 | 0.001484 | No | Yes |
| 211234_x_at | ESR1 | 0.3448 | 0.0009249 | No | Yes |
| 211235_s_at | ESR1 | 0.3405 | 0.001212 | No | Yes |
| 211627_x_at | ESR1 | 0.6322 | 0.001034 | No | Yes |
| 215551_at | ESR1 | 0.7676 | 0.7336 | No | No |
| 215552_s_at | ESR1 | 0.02772 | 0.001934 | Yes | Yes |
| 217163_at | ESR1 | 0.956 | 0.1695 | No | No |
| 217190_x_at | ESR1 | 0.1205 | 0.08792 | No | No |
| 204053_x_at | PTEN | 2.37e-13 | 0.1861 | Yes | No |
| 204054_at | PTEN | 9.37e-09 | 0.006888 | Yes | Yes |
| 211711_s_at | PTEN | 1.46e-15 | 0.001212 | Yes | Yes |
| 217492_s_at | PTEN | 4.79e-14 | 0.04205 | Yes | Yes |
| 225363_at | PTEN | 3.48e-05 | 3.98e-06 | Yes | Yes |
| 228006_at | PTEN | 9.50e-13 | 4.32e-09 | Yes | Yes |
| 233254_x_at | PTEN | 0.3086 | 0.9757 | No | No |
| 242622_x_at | PTEN | 0.0594 | 0.004498 | No | Yes |
| 1565483_at | EGFR | 0.7676 | 0.4642 | No | No |
| 1565484_x_at | EGFR | 0.475 | 0.1729 | No | No |
| 201983_s_at | EGFR | 1.41e-05 | 0.001164 | Yes | Yes |
| 201984_s_at | EGFR | 2.10e-08 | 0.01963 | Yes | Yes |
| 210984_x_at | EGFR | 0.0002481 | 0.2101 | Yes | No |
| 211550_at | EGFR | 0.9131 | 0.001484 | No | Yes |
| 211551_at | EGFR | 0.1205 | 0.03735 | No | Yes |
| 211607_x_at | EGFR | 0.02808 | 0.2783 | Yes | No |
| 224999_at | EGFR | 1.73e-08 | 0.0005327 | Yes | Yes |

### Tumour versus healthy in GSE65194

The 38-probe Welch analysis identified **27/38** target probes at BH-adjusted p < 0.05. Each full-probe limma fit identified **33,763** (repA) and **33,552** (repB) probes of 54,673 at BH-adjusted p < 0.05; within the target panel, **19/38** and **19/38**, respectively, were significant. The limma–Welch sensitivity report records 38/38 directional agreement for each assignment. The different BH families (38 versus 54,673) matter when comparing significance.

The 38 target probes by gene are listed below. The per-gene counts summarize individual probe directions/significance; they are not gene-level expression estimates.

| Target gene | Probes | Tumour higher (Welch) | Tumour lower (Welch) | Welch significant | limma repA significant | limma repB significant |
| --- | --- | --- | --- | --- | --- | --- |
| 204531_s_at | 5 | 0 | 0 | 0 | 0 | 0 |
| 211851_x_at | 5 | 0 | 0 | 0 | 0 | 0 |
| 208368_s_at | 5 | 0 | 0 | 0 | 0 | 0 |
| 214727_at | 5 | 0 | 0 | 0 | 0 | 0 |
| 201746_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211300_s_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 210930_s_at | 5 | 0 | 0 | 0 | 0 | 0 |
| 216836_s_at | 5 | 0 | 0 | 0 | 0 | 0 |
| 234354_x_at | 5 | 0 | 0 | 0 | 0 | 0 |
| 204369_at | 6 | 0 | 0 | 0 | 0 | 0 |
| 231854_at | 6 | 0 | 0 | 0 | 0 | 0 |
| 235980_at | 6 | 0 | 0 | 0 | 0 | 0 |
| 205225_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211233_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211234_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211235_s_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211627_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 215551_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 215552_s_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 217163_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 217190_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 204053_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 204054_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211711_s_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 217492_s_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 225363_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 228006_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 233254_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 242622_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 1565483_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 1565484_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 201983_s_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 201984_s_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 210984_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211550_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211551_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 211607_x_at | 4 | 0 | 0 | 0 | 0 | 0 |
| 224999_at | 4 | 0 | 0 | 0 | 0 | 0 |

| Probe | Gene | Welch difference | Welch BH p (38) | limma repA logFC | repA BH p (54,673) | limma repB logFC | repB BH p (54,673) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 204531_s_at | BRCA1 | 0.4305 | 0.05263 | 0.4316 | 0.2072 | 0.4294 | 0.2097 |
| 211851_x_at | BRCA1 | 0.1564 | 0.004399 | 0.1565 | 0.3319 | 0.1564 | 0.3622 |
| 208368_s_at | BRCA2 | 0.7524 | 0.001212 | 0.7510 | 0.01268 | 0.7538 | 0.01364 |
| 214727_at | BRCA2 | 0.2055 | 0.0122 | 0.2026 | 0.3693 | 0.2085 | 0.3662 |
| 201746_at | TP53 | -0.8253 | 0.0005365 | -0.8251 | 0.0375 | -0.8255 | 0.03774 |
| 211300_s_at | TP53 | 1.2248 | 1.86e-07 | 1.2247 | 0.004015 | 1.2248 | 0.003946 |
| 210930_s_at | ERBB2 | 2.0906 | 4.32e-09 | 2.0857 | 0.02761 | 2.0956 | 0.03026 |
| 216836_s_at | ERBB2 | 1.8863 | 2.04e-06 | 1.8864 | 0.02097 | 1.8861 | 0.02153 |
| 234354_x_at | ERBB2 | 1.2108 | 2.76e-07 | 1.2084 | 0.0594 | 1.2131 | 0.06157 |
| 204369_at | PIK3CA | 2.9393 | 0.0009249 | 2.9377 | 1.75e-16 | 2.9408 | 9.68e-17 |
| 231854_at | PIK3CA | -0.3928 | 0.01963 | -0.4007 | 1.67e-06 | -0.3849 | 0.0003339 |
| 235980_at | PIK3CA | -1.3486 | 0.0004011 | -1.3524 | 5.24e-12 | -1.3448 | 4.85e-12 |
| 205225_at | ESR1 | -0.9293 | 0.1861 | -0.9277 | 0.5713 | -0.9310 | 0.5721 |
| 211233_x_at | ESR1 | 0.1449 | 0.001484 | 0.1449 | 0.3967 | 0.1449 | 0.3988 |
| 211234_x_at | ESR1 | 0.1524 | 0.0009249 | 0.1524 | 0.3466 | 0.1524 | 0.3483 |
| 211235_s_at | ESR1 | 0.3442 | 0.001212 | 0.3442 | 0.2781 | 0.3442 | 0.28 |
| 211627_x_at | ESR1 | -0.8792 | 0.001034 | -0.8789 | 2.03e-20 | -0.8794 | 4.75e-21 |
| 215551_at | ESR1 | -0.0253 | 0.7336 | -0.0253 | 0.8886 | -0.0253 | 0.8904 |
| 215552_s_at | ESR1 | 0.6748 | 0.001934 | 0.6748 | 0.1986 | 0.6748 | 0.1997 |
| 217163_at | ESR1 | -0.0071 | 0.1695 | -0.0071 | 0.1432 | -0.0072 | 0.136 |
| 217190_x_at | ESR1 | 0.0532 | 0.08792 | 0.0532 | 0.6151 | 0.0532 | 0.6179 |
| 204053_x_at | PTEN | 0.3799 | 0.1861 | 0.3802 | 0.2002 | 0.3797 | 0.1917 |
| 204054_at | PTEN | -0.3069 | 0.006888 | -0.3066 | 0.2326 | -0.3072 | 0.2517 |
| 211711_s_at | PTEN | -0.9506 | 0.001212 | -0.9500 | 0.002491 | -0.9512 | 0.003098 |
| 217492_s_at | PTEN | 0.7118 | 0.04205 | 0.7120 | 0.01157 | 0.7116 | 0.009763 |
| 225363_at | PTEN | -1.0096 | 3.98e-06 | -1.0094 | 7.81e-05 | -1.0097 | 5.77e-05 |
| 228006_at | PTEN | -1.5964 | 4.32e-09 | -1.5962 | 1.54e-05 | -1.5967 | 1.65e-05 |
| 233254_x_at | PTEN | 0.0065 | 0.9757 | 0.0076 | 0.9786 | 0.0054 | 0.984 |
| 242622_x_at | PTEN | -0.6158 | 0.004498 | -0.6106 | 1.71e-06 | -0.6210 | 5.22e-06 |
| 1565483_at | EGFR | -0.0381 | 0.4642 | -0.0422 | 0.5231 | -0.0340 | 0.665 |
| 1565484_x_at | EGFR | -0.2185 | 0.1729 | -0.2219 | 7.41e-07 | -0.2151 | 0.0002265 |
| 201983_s_at | EGFR | -1.8977 | 0.001164 | -1.8973 | 0.007189 | -1.8981 | 0.007076 |
| 201984_s_at | EGFR | -0.9937 | 0.01963 | -0.9931 | 0.06502 | -0.9942 | 0.06821 |
| 210984_x_at | EGFR | 0.0960 | 0.2101 | 0.0964 | 0.6798 | 0.0956 | 0.6882 |
| 211550_at | EGFR | -0.1512 | 0.001484 | -0.1515 | 8.24e-08 | -0.1508 | 5.63e-09 |
| 211551_at | EGFR | -0.0262 | 0.03735 | -0.0262 | 3.03e-05 | -0.0262 | 2.32e-05 |
| 211607_x_at | EGFR | 0.0778 | 0.2783 | 0.0782 | 0.7228 | 0.0774 | 0.7315 |
| 224999_at | EGFR | -1.9580 | 0.0005327 | -1.9576 | 0.0003464 | -1.9584 | 0.0002528 |

### limma repA/repB sensitivity

repA and repB are alternative selections from 23 paired patients. The probe-level comparison includes each logFC, difference, raw p-value and full-matrix BH p-value.

| Probe | Gene | repA logFC | repB logFC | repA−repB | repA BH p | repB BH p | repA significant | repB significant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 204531_s_at | BRCA1 | 0.4316 | 0.4294 | 0.0022 | 0.2072 | 0.2097 | No | No |
| 211851_x_at | BRCA1 | 0.1565 | 0.1564 | 0.0002 | 0.3319 | 0.3622 | No | No |
| 208368_s_at | BRCA2 | 0.7510 | 0.7538 | -0.0028 | 0.01268 | 0.01364 | Yes | Yes |
| 214727_at | BRCA2 | 0.2026 | 0.2085 | -0.0059 | 0.3693 | 0.3662 | No | No |
| 201746_at | TP53 | -0.8251 | -0.8255 | 0.0005 | 0.0375 | 0.03774 | Yes | Yes |
| 211300_s_at | TP53 | 1.2247 | 1.2248 | -0.0001 | 0.004015 | 0.003946 | Yes | Yes |
| 210930_s_at | ERBB2 | 2.0857 | 2.0956 | -0.0099 | 0.02761 | 0.03026 | Yes | Yes |
| 216836_s_at | ERBB2 | 1.8864 | 1.8861 | 0.0003 | 0.02097 | 0.02153 | Yes | Yes |
| 234354_x_at | ERBB2 | 1.2084 | 1.2131 | -0.0046 | 0.0594 | 0.06157 | No | No |
| 204369_at | PIK3CA | 2.9377 | 2.9408 | -0.0031 | 1.75e-16 | 9.68e-17 | Yes | Yes |
| 231854_at | PIK3CA | -0.4007 | -0.3849 | -0.0157 | 1.67e-06 | 0.0003339 | Yes | Yes |
| 235980_at | PIK3CA | -1.3524 | -1.3448 | -0.0075 | 5.24e-12 | 4.85e-12 | Yes | Yes |
| 205225_at | ESR1 | -0.9277 | -0.9310 | 0.0033 | 0.5713 | 0.5721 | No | No |
| 211233_x_at | ESR1 | 0.1449 | 0.1449 | 0.0000 | 0.3967 | 0.3988 | No | No |
| 211234_x_at | ESR1 | 0.1524 | 0.1524 | 0.0000 | 0.3466 | 0.3483 | No | No |
| 211235_s_at | ESR1 | 0.3442 | 0.3442 | -0.0000 | 0.2781 | 0.28 | No | No |
| 211627_x_at | ESR1 | -0.8789 | -0.8794 | 0.0006 | 2.03e-20 | 4.75e-21 | Yes | Yes |
| 215551_at | ESR1 | -0.0253 | -0.0253 | 0.0000 | 0.8886 | 0.8904 | No | No |
| 215552_s_at | ESR1 | 0.6748 | 0.6748 | 0.0000 | 0.1986 | 0.1997 | No | No |
| 217163_at | ESR1 | -0.0071 | -0.0072 | 0.0001 | 0.1432 | 0.136 | No | No |
| 217190_x_at | ESR1 | 0.0532 | 0.0532 | 0.0000 | 0.6151 | 0.6179 | No | No |
| 204053_x_at | PTEN | 0.3802 | 0.3797 | 0.0005 | 0.2002 | 0.1917 | No | No |
| 204054_at | PTEN | -0.3066 | -0.3072 | 0.0005 | 0.2326 | 0.2517 | No | No |
| 211711_s_at | PTEN | -0.9500 | -0.9512 | 0.0012 | 0.002491 | 0.003098 | Yes | Yes |
| 217492_s_at | PTEN | 0.7120 | 0.7116 | 0.0004 | 0.01157 | 0.009763 | Yes | Yes |
| 225363_at | PTEN | -1.0094 | -1.0097 | 0.0004 | 7.81e-05 | 5.77e-05 | Yes | Yes |
| 228006_at | PTEN | -1.5962 | -1.5967 | 0.0005 | 1.54e-05 | 1.65e-05 | Yes | Yes |
| 233254_x_at | PTEN | 0.0076 | 0.0054 | 0.0022 | 0.9786 | 0.984 | No | No |
| 242622_x_at | PTEN | -0.6106 | -0.6210 | 0.0104 | 1.71e-06 | 5.22e-06 | Yes | Yes |
| 1565483_at | EGFR | -0.0422 | -0.0340 | -0.0082 | 0.5231 | 0.665 | No | No |
| 1565484_x_at | EGFR | -0.2219 | -0.2151 | -0.0068 | 7.41e-07 | 0.0002265 | Yes | Yes |
| 201983_s_at | EGFR | -1.8973 | -1.8981 | 0.0008 | 0.007189 | 0.007076 | Yes | Yes |
| 201984_s_at | EGFR | -0.9931 | -0.9942 | 0.0010 | 0.06502 | 0.06821 | No | No |
| 210984_x_at | EGFR | 0.0964 | 0.0956 | 0.0007 | 0.6798 | 0.6882 | No | No |
| 211550_at | EGFR | -0.1515 | -0.1508 | -0.0007 | 8.24e-08 | 5.63e-09 | Yes | Yes |
| 211551_at | EGFR | -0.0262 | -0.0262 | 0.0000 | 3.03e-05 | 2.32e-05 | Yes | Yes |
| 211607_x_at | EGFR | 0.0782 | 0.0774 | 0.0008 | 0.7228 | 0.7315 | No | No |
| 224999_at | EGFR | -1.9576 | -1.9584 | 0.0008 | 0.0003464 | 0.0002528 | Yes | Yes |

### Patient-level subtype comparisons

The patient-level Welch and limma files each contain all 228 probe-by-contrast comparisons. Overall, limma identified **100/228** significant tests and Welch **95/228**. **87** were significant by both, **13** by limma only, **8** by Welch only and **120** by neither. Direction agreed for **228/228** tests. The table gives per-contrast counts.

| Contrast | limma significant | Welch significant | Both | limma only | Welch only | Neither | Direction agrees |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TNBC - Her2 | 18 | 12 | 12 | 6 | 0 | 20 | 38/38 |
| TNBC - Luminal A | 24 | 19 | 19 | 5 | 0 | 14 | 38/38 |
| TNBC - Luminal B | 21 | 22 | 20 | 1 | 2 | 15 | 38/38 |
| Her2 - Luminal A | 17 | 17 | 16 | 1 | 1 | 20 | 38/38 |
| Her2 - Luminal B | 13 | 15 | 13 | 0 | 2 | 23 | 38/38 |
| Luminal A - Luminal B | 7 | 10 | 7 | 0 | 3 | 28 | 38/38 |

### ESR1 and ERBB2 subtype patterns

Across the six contrasts, the number of ERBB2 probes higher in the first-named subtype is: TNBC - Her2: 0/3 higher in the first subtype; TNBC - Luminal A: 1/3 higher in the first subtype; TNBC - Luminal B: 0/3 higher in the first subtype; Her2 - Luminal A: 3/3 higher in the first subtype; Her2 - Luminal B: 3/3 higher in the first subtype; Luminal A - Luminal B: 0/3 higher in the first subtype. The exact-probe results and adjusted p-values below retain each probe separately. These group-level differences do not establish receptor status for individual samples.

ESR1 has nine probes with heterogeneous results rather than one interchangeable gene-level value. In the Luminal A − Luminal B contrast, 0 limma and 0 Welch ESR1 probes meet adjusted p < 0.05; for other contrasts, significance varies by probe. Exact signs and test counts are shown below without inferring clinical ER-positive status from Luminal labels.

| Gene | Contrast | Probe directions | limma/Welch significant probes | Direction agreement |
| --- | --- | --- | --- | --- |
| ESR1 | TNBC - Her2 | 1 higher in first / 8 higher in second | 0/0 | 9/9 |
| ESR1 | TNBC - Luminal A | 1 higher in first / 8 higher in second | 7/6 | 9/9 |
| ESR1 | TNBC - Luminal B | 1 higher in first / 8 higher in second | 5/6 | 9/9 |
| ESR1 | Her2 - Luminal A | 1 higher in first / 8 higher in second | 6/5 | 9/9 |
| ESR1 | Her2 - Luminal B | 1 higher in first / 8 higher in second | 5/5 | 9/9 |
| ESR1 | Luminal A - Luminal B | 7 higher in first / 2 higher in second | 0/0 | 9/9 |
| ERBB2 | TNBC - Her2 | 0 higher in first / 3 higher in second | 3/3 | 3/3 |
| ERBB2 | TNBC - Luminal A | 1 higher in first / 2 higher in second | 1/1 | 3/3 |
| ERBB2 | TNBC - Luminal B | 0 higher in first / 3 higher in second | 3/3 | 3/3 |
| ERBB2 | Her2 - Luminal A | 3 higher in first / 0 higher in second | 3/3 | 3/3 |
| ERBB2 | Her2 - Luminal B | 3 higher in first / 0 higher in second | 3/3 | 3/3 |
| ERBB2 | Luminal A - Luminal B | 0 higher in first / 3 higher in second | 3/3 | 3/3 |

#### ESR1 and ERBB2 probe-level subtype results

| Contrast | Probe | Gene | limma logFC | Welch difference | limma BH p (228) | Welch BH p (228) | Direction agrees | Significance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TNBC - Her2 | 210930_s_at | ERBB2 | -6.2394 | -6.2394 | 6.19e-38 | 5.33e-31 | Yes | Both |
| TNBC - Her2 | 216836_s_at | ERBB2 | -5.1905 | -5.1905 | 1.59e-38 | 5.33e-31 | Yes | Both |
| TNBC - Her2 | 234354_x_at | ERBB2 | -3.5398 | -3.5398 | 2.53e-21 | 1.66e-11 | Yes | Both |
| TNBC - Her2 | 205225_at | ESR1 | -0.3468 | -0.3468 | 0.571 | 0.5826 | Yes | Neither |
| TNBC - Her2 | 211233_x_at | ESR1 | -0.0010 | -0.0010 | 0.9979 | 0.6817 | Yes | Neither |
| TNBC - Her2 | 211234_x_at | ESR1 | -0.0025 | -0.0025 | 0.9979 | 0.4089 | Yes | Neither |
| TNBC - Her2 | 211235_s_at | ESR1 | -0.0038 | -0.0038 | 0.9979 | 0.7657 | Yes | Neither |
| TNBC - Her2 | 211627_x_at | ESR1 | -0.0014 | -0.0014 | 0.9979 | 0.989 | Yes | Neither |
| TNBC - Her2 | 215551_at | ESR1 | 0.0060 | 0.0060 | 0.9979 | 0.4605 | Yes | Neither |
| TNBC - Her2 | 215552_s_at | ESR1 | -0.0163 | -0.0163 | 0.9979 | 0.5098 | Yes | Neither |
| TNBC - Her2 | 217163_at | ESR1 | -0.0008 | -0.0008 | 0.8535 | 0.7639 | Yes | Neither |
| TNBC - Her2 | 217190_x_at | ESR1 | -0.0388 | -0.0388 | 0.6916 | 0.2398 | Yes | Neither |
| TNBC - Luminal A | 210930_s_at | ERBB2 | -0.3265 | -0.3265 | 0.4618 | 0.1288 | Yes | Neither |
| TNBC - Luminal A | 216836_s_at | ERBB2 | -1.0428 | -1.0428 | 0.0005336 | 0.0002557 | Yes | Both |
| TNBC - Luminal A | 234354_x_at | ERBB2 | 0.0008 | 0.0008 | 0.9979 | 0.9955 | Yes | Neither |
| TNBC - Luminal A | 205225_at | ESR1 | -8.2106 | -8.2106 | 9.28e-38 | 3.70e-28 | Yes | Both |
| TNBC - Luminal A | 211233_x_at | ESR1 | -0.4538 | -0.4538 | 0.0002077 | 0.01127 | Yes | Both |
| TNBC - Luminal A | 211234_x_at | ESR1 | -0.4696 | -0.4696 | 3.65e-05 | 0.005835 | Yes | Both |
| TNBC - Luminal A | 211235_s_at | ESR1 | -1.1887 | -1.1887 | 1.21e-08 | 0.0001591 | Yes | Both |
| TNBC - Luminal A | 211627_x_at | ESR1 | -0.0553 | -0.0553 | 0.3868 | 0.5228 | Yes | Neither |
| TNBC - Luminal A | 215551_at | ESR1 | -0.2382 | -0.2382 | 0.04365 | 0.1312 | Yes | Limma only |
| TNBC - Luminal A | 215552_s_at | ESR1 | -2.2449 | -2.2449 | 1.14e-12 | 1.61e-06 | Yes | Both |
| TNBC - Luminal A | 217163_at | ESR1 | 0.0049 | 0.0049 | 0.169 | 0.1222 | Yes | Neither |
| TNBC - Luminal A | 217190_x_at | ESR1 | -0.2454 | -0.2454 | 0.0008457 | 0.02992 | Yes | Both |
| TNBC - Luminal B | 210930_s_at | ERBB2 | -3.3937 | -3.3937 | 5.67e-18 | 3.62e-08 | Yes | Both |
| TNBC - Luminal B | 216836_s_at | ERBB2 | -4.2681 | -4.2681 | 1.19e-31 | 5.14e-17 | Yes | Both |
| TNBC - Luminal B | 234354_x_at | ERBB2 | -1.8880 | -1.8880 | 1.97e-08 | 6.40e-05 | Yes | Both |
| TNBC - Luminal B | 205225_at | ESR1 | -7.9556 | -7.9556 | 5.39e-37 | 2.21e-27 | Yes | Both |
| TNBC - Luminal B | 211233_x_at | ESR1 | -0.2157 | -0.2157 | 0.08825 | 0.0676 | Yes | Neither |
| TNBC - Luminal B | 211234_x_at | ESR1 | -0.2931 | -0.2931 | 0.01071 | 0.007954 | Yes | Both |
| TNBC - Luminal B | 211235_s_at | ESR1 | -0.7444 | -0.7444 | 0.0003273 | 0.0008979 | Yes | Both |
| TNBC - Luminal B | 211627_x_at | ESR1 | -0.0192 | -0.0192 | 0.8163 | 0.7044 | Yes | Neither |
| TNBC - Luminal B | 215551_at | ESR1 | -0.3215 | -0.3215 | 0.005284 | 0.009071 | Yes | Both |
| TNBC - Luminal B | 215552_s_at | ESR1 | -1.9061 | -1.9061 | 5.00e-10 | 1.49e-06 | Yes | Both |
| TNBC - Luminal B | 217163_at | ESR1 | 0.0024 | 0.0024 | 0.567 | 0.4525 | Yes | Neither |
| TNBC - Luminal B | 217190_x_at | ESR1 | -0.1063 | -0.1063 | 0.1751 | 0.04448 | Yes | Welch only |
| Her2 - Luminal A | 210930_s_at | ERBB2 | 5.9129 | 5.9129 | 4.81e-33 | 7.04e-27 | Yes | Both |
| Her2 - Luminal A | 216836_s_at | ERBB2 | 4.1477 | 4.1477 | 9.45e-28 | 1.32e-26 | Yes | Both |
| Her2 - Luminal A | 234354_x_at | ERBB2 | 3.5406 | 3.5406 | 4.05e-19 | 1.27e-11 | Yes | Both |
| Her2 - Luminal A | 205225_at | ESR1 | -7.8638 | -7.8638 | 2.00e-33 | 4.11e-23 | Yes | Both |
| Her2 - Luminal A | 211233_x_at | ESR1 | -0.4528 | -0.4528 | 0.0005512 | 0.01132 | Yes | Both |
| Her2 - Luminal A | 211234_x_at | ESR1 | -0.4671 | -0.4671 | 0.0001414 | 0.005929 | Yes | Both |
| Her2 - Luminal A | 211235_s_at | ESR1 | -1.1849 | -1.1849 | 1.10e-07 | 0.000162 | Yes | Both |
| Her2 - Luminal A | 211627_x_at | ESR1 | -0.0539 | -0.0539 | 0.4458 | 0.5312 | Yes | Neither |
| Her2 - Luminal A | 215551_at | ESR1 | -0.2442 | -0.2442 | 0.05608 | 0.1222 | Yes | Neither |
| Her2 - Luminal A | 215552_s_at | ESR1 | -2.2285 | -2.2285 | 2.92e-11 | 1.70e-06 | Yes | Both |
| Her2 - Luminal A | 217163_at | ESR1 | 0.0057 | 0.0057 | 0.127 | 0.07919 | Yes | Neither |
| Her2 - Luminal A | 217190_x_at | ESR1 | -0.2067 | -0.2067 | 0.01024 | 0.07426 | Yes | Limma only |
| Her2 - Luminal B | 210930_s_at | ERBB2 | 2.8457 | 2.8457 | 1.87e-12 | 1.61e-06 | Yes | Both |
| Her2 - Luminal B | 216836_s_at | ERBB2 | 0.9224 | 0.9224 | 0.004527 | 0.009194 | Yes | Both |
| Her2 - Luminal B | 234354_x_at | ERBB2 | 1.6519 | 1.6519 | 4.64e-06 | 0.00305 | Yes | Both |
| Her2 - Luminal B | 205225_at | ESR1 | -7.6088 | -7.6088 | 1.24e-32 | 1.38e-22 | Yes | Both |
| Her2 - Luminal B | 211233_x_at | ESR1 | -0.2147 | -0.2147 | 0.1179 | 0.0676 | Yes | Neither |
| Her2 - Luminal B | 211234_x_at | ESR1 | -0.2907 | -0.2907 | 0.01956 | 0.008176 | Yes | Both |
| Her2 - Luminal B | 211235_s_at | ESR1 | -0.7407 | -0.7407 | 0.0008592 | 0.0009363 | Yes | Both |
| Her2 - Luminal B | 211627_x_at | ESR1 | -0.0178 | -0.0178 | 0.848 | 0.7283 | Yes | Neither |
| Her2 - Luminal B | 215551_at | ESR1 | -0.3275 | -0.3275 | 0.008753 | 0.008002 | Yes | Both |
| Her2 - Luminal B | 215552_s_at | ESR1 | -1.8898 | -1.8898 | 7.63e-09 | 1.61e-06 | Yes | Both |
| Her2 - Luminal B | 217163_at | ESR1 | 0.0032 | 0.0032 | 0.4415 | 0.3272 | Yes | Neither |
| Her2 - Luminal B | 217190_x_at | ESR1 | -0.0675 | -0.0675 | 0.4712 | 0.286 | Yes | Neither |
| Luminal A - Luminal B | 210930_s_at | ERBB2 | -3.0672 | -3.0672 | 9.86e-14 | 3.12e-07 | Yes | Both |
| Luminal A - Luminal B | 216836_s_at | ERBB2 | -3.2254 | -3.2254 | 6.58e-20 | 2.38e-12 | Yes | Both |
| Luminal A - Luminal B | 234354_x_at | ERBB2 | -1.8888 | -1.8888 | 2.04e-07 | 6.40e-05 | Yes | Both |
| Luminal A - Luminal B | 205225_at | ESR1 | 0.2550 | 0.2550 | 0.7262 | 0.5974 | Yes | Neither |
| Luminal A - Luminal B | 211233_x_at | ESR1 | 0.2381 | 0.2381 | 0.08179 | 0.2753 | Yes | Neither |
| Luminal A - Luminal B | 211234_x_at | ESR1 | 0.1765 | 0.1765 | 0.1809 | 0.4077 | Yes | Neither |
| Luminal A - Luminal B | 211235_s_at | ESR1 | 0.4443 | 0.4443 | 0.05608 | 0.2242 | Yes | Neither |
| Luminal A - Luminal B | 211627_x_at | ESR1 | 0.0361 | 0.0361 | 0.63 | 0.7044 | Yes | Neither |
| Luminal A - Luminal B | 215551_at | ESR1 | -0.0833 | -0.0833 | 0.5972 | 0.7044 | Yes | Neither |
| Luminal A - Luminal B | 215552_s_at | ESR1 | 0.3387 | 0.3387 | 0.3761 | 0.5458 | Yes | Neither |
| Luminal A - Luminal B | 217163_at | ESR1 | -0.0025 | -0.0025 | 0.5841 | 0.5312 | Yes | Neither |
| Luminal A - Luminal B | 217190_x_at | ESR1 | 0.1391 | 0.1391 | 0.09007 | 0.2753 | Yes | Neither |

#### Complete 228-row subtype result

All target genes and probes, all six contrast directions, both methods' adjusted p-values, and significance categories are preserved in this table.

| Contrast | Probe | Gene | limma logFC | Welch difference | limma BH p (228) | Welch BH p (228) | Direction agrees | Significance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TNBC - Her2 | 204531_s_at | BRCA1 | 0.1227 | 0.1227 | 0.7274 | 0.7044 | Yes | Neither |
| TNBC - Her2 | 211851_x_at | BRCA1 | 0.2333 | 0.2333 | 0.0704 | 0.08404 | Yes | Neither |
| TNBC - Her2 | 208368_s_at | BRCA2 | 0.2307 | 0.2307 | 0.3481 | 0.4077 | Yes | Neither |
| TNBC - Her2 | 214727_at | BRCA2 | 0.0435 | 0.0435 | 0.8535 | 0.8737 | Yes | Neither |
| TNBC - Her2 | 201746_at | TP53 | 0.6825 | 0.6825 | 0.03664 | 0.09481 | Yes | Limma only |
| TNBC - Her2 | 211300_s_at | TP53 | 0.7483 | 0.7483 | 0.03002 | 0.0676 | Yes | Limma only |
| TNBC - Her2 | 210930_s_at | ERBB2 | -6.2394 | -6.2394 | 6.19e-38 | 5.33e-31 | Yes | Both |
| TNBC - Her2 | 216836_s_at | ERBB2 | -5.1905 | -5.1905 | 1.59e-38 | 5.33e-31 | Yes | Both |
| TNBC - Her2 | 234354_x_at | ERBB2 | -3.5398 | -3.5398 | 2.53e-21 | 1.66e-11 | Yes | Both |
| TNBC - Her2 | 204369_at | PIK3CA | 0.6168 | 0.6168 | 0.006815 | 0.008002 | Yes | Both |
| TNBC - Her2 | 231854_at | PIK3CA | 0.1697 | 0.1697 | 0.01185 | 0.04357 | Yes | Both |
| TNBC - Her2 | 235980_at | PIK3CA | 0.3152 | 0.3152 | 0.02989 | 0.05782 | Yes | Limma only |
| TNBC - Her2 | 205225_at | ESR1 | -0.3468 | -0.3468 | 0.571 | 0.5826 | Yes | Neither |
| TNBC - Her2 | 211233_x_at | ESR1 | -0.0010 | -0.0010 | 0.9979 | 0.6817 | Yes | Neither |
| TNBC - Her2 | 211234_x_at | ESR1 | -0.0025 | -0.0025 | 0.9979 | 0.4089 | Yes | Neither |
| TNBC - Her2 | 211235_s_at | ESR1 | -0.0038 | -0.0038 | 0.9979 | 0.7657 | Yes | Neither |
| TNBC - Her2 | 211627_x_at | ESR1 | -0.0014 | -0.0014 | 0.9979 | 0.989 | Yes | Neither |
| TNBC - Her2 | 215551_at | ESR1 | 0.0060 | 0.0060 | 0.9979 | 0.4605 | Yes | Neither |
| TNBC - Her2 | 215552_s_at | ESR1 | -0.0163 | -0.0163 | 0.9979 | 0.5098 | Yes | Neither |
| TNBC - Her2 | 217163_at | ESR1 | -0.0008 | -0.0008 | 0.8535 | 0.7639 | Yes | Neither |
| TNBC - Her2 | 217190_x_at | ESR1 | -0.0388 | -0.0388 | 0.6916 | 0.2398 | Yes | Neither |
| TNBC - Her2 | 204053_x_at | PTEN | -0.7402 | -0.7402 | 0.0004666 | 0.002937 | Yes | Both |
| TNBC - Her2 | 204054_at | PTEN | -0.3472 | -0.3472 | 0.09811 | 0.1492 | Yes | Neither |
| TNBC - Her2 | 211711_s_at | PTEN | -0.7213 | -0.7213 | 0.004162 | 0.009071 | Yes | Both |
| TNBC - Her2 | 217492_s_at | PTEN | -0.6856 | -0.6856 | 0.001011 | 0.004151 | Yes | Both |
| TNBC - Her2 | 225363_at | PTEN | -0.7018 | -0.7018 | 7.87e-06 | 0.0001634 | Yes | Both |
| TNBC - Her2 | 228006_at | PTEN | -0.8402 | -0.8402 | 0.0005336 | 0.00374 | Yes | Both |
| TNBC - Her2 | 233254_x_at | PTEN | -0.5562 | -0.5562 | 0.0003143 | 0.001557 | Yes | Both |
| TNBC - Her2 | 242622_x_at | PTEN | -0.2262 | -0.2262 | 0.023 | 0.05373 | Yes | Limma only |
| TNBC - Her2 | 1565483_at | EGFR | -0.0212 | -0.0212 | 0.7676 | 0.4503 | Yes | Neither |
| TNBC - Her2 | 1565484_x_at | EGFR | -0.0051 | -0.0051 | 0.8535 | 0.8768 | Yes | Neither |
| TNBC - Her2 | 201983_s_at | EGFR | 1.0618 | 1.0618 | 0.026 | 0.04763 | Yes | Both |
| TNBC - Her2 | 201984_s_at | EGFR | 0.7323 | 0.7323 | 0.03785 | 0.1372 | Yes | Limma only |
| TNBC - Her2 | 210984_x_at | EGFR | -0.1442 | -0.1442 | 0.4696 | 0.6211 | Yes | Neither |
| TNBC - Her2 | 211550_at | EGFR | 0.0061 | 0.0061 | 0.8471 | 0.678 | Yes | Neither |
| TNBC - Her2 | 211551_at | EGFR | -0.0005 | -0.0005 | 0.9559 | 0.8737 | Yes | Neither |
| TNBC - Her2 | 211607_x_at | EGFR | -0.1125 | -0.1125 | 0.5709 | 0.6817 | Yes | Neither |
| TNBC - Her2 | 224999_at | EGFR | 0.8134 | 0.8134 | 0.03664 | 0.08603 | Yes | Limma only |
| TNBC - Luminal A | 204531_s_at | BRCA1 | 0.8158 | 0.8158 | 0.001497 | 0.004371 | Yes | Both |
| TNBC - Luminal A | 211851_x_at | BRCA1 | 0.3192 | 0.3192 | 0.0116 | 0.01202 | Yes | Both |
| TNBC - Luminal A | 208368_s_at | BRCA2 | 1.2358 | 1.2358 | 1.04e-08 | 3.23e-09 | Yes | Both |
| TNBC - Luminal A | 214727_at | BRCA2 | 0.4484 | 0.4484 | 0.009255 | 0.004396 | Yes | Both |
| TNBC - Luminal A | 201746_at | TP53 | 0.2064 | 0.2064 | 0.6185 | 0.5312 | Yes | Neither |
| TNBC - Luminal A | 211300_s_at | TP53 | 0.6483 | 0.6483 | 0.06592 | 0.07129 | Yes | Neither |
| TNBC - Luminal A | 210930_s_at | ERBB2 | -0.3265 | -0.3265 | 0.4618 | 0.1288 | Yes | Neither |
| TNBC - Luminal A | 216836_s_at | ERBB2 | -1.0428 | -1.0428 | 0.0005336 | 0.0002557 | Yes | Both |
| TNBC - Luminal A | 234354_x_at | ERBB2 | 0.0008 | 0.0008 | 0.9979 | 0.9955 | Yes | Neither |
| TNBC - Luminal A | 204369_at | PIK3CA | 0.3581 | 0.3581 | 0.1365 | 0.06872 | Yes | Neither |
| TNBC - Luminal A | 231854_at | PIK3CA | 0.0927 | 0.0927 | 0.2017 | 0.2133 | Yes | Neither |
| TNBC - Luminal A | 235980_at | PIK3CA | -0.0397 | -0.0397 | 0.8535 | 0.8 | Yes | Neither |
| TNBC - Luminal A | 205225_at | ESR1 | -8.2106 | -8.2106 | 9.28e-38 | 3.70e-28 | Yes | Both |
| TNBC - Luminal A | 211233_x_at | ESR1 | -0.4538 | -0.4538 | 0.0002077 | 0.01127 | Yes | Both |
| TNBC - Luminal A | 211234_x_at | ESR1 | -0.4696 | -0.4696 | 3.65e-05 | 0.005835 | Yes | Both |
| TNBC - Luminal A | 211235_s_at | ESR1 | -1.1887 | -1.1887 | 1.21e-08 | 0.0001591 | Yes | Both |
| TNBC - Luminal A | 211627_x_at | ESR1 | -0.0553 | -0.0553 | 0.3868 | 0.5228 | Yes | Neither |
| TNBC - Luminal A | 215551_at | ESR1 | -0.2382 | -0.2382 | 0.04365 | 0.1312 | Yes | Limma only |
| TNBC - Luminal A | 215552_s_at | ESR1 | -2.2449 | -2.2449 | 1.14e-12 | 1.61e-06 | Yes | Both |
| TNBC - Luminal A | 217163_at | ESR1 | 0.0049 | 0.0049 | 0.169 | 0.1222 | Yes | Neither |
| TNBC - Luminal A | 217190_x_at | ESR1 | -0.2454 | -0.2454 | 0.0008457 | 0.02992 | Yes | Both |
| TNBC - Luminal A | 204053_x_at | PTEN | -0.4486 | -0.4486 | 0.04088 | 0.06498 | Yes | Limma only |
| TNBC - Luminal A | 204054_at | PTEN | -0.4209 | -0.4209 | 0.04365 | 0.05811 | Yes | Limma only |
| TNBC - Luminal A | 211711_s_at | PTEN | -0.7091 | -0.7091 | 0.005284 | 0.009801 | Yes | Both |
| TNBC - Luminal A | 217492_s_at | PTEN | -0.4598 | -0.4598 | 0.03495 | 0.05296 | Yes | Limma only |
| TNBC - Luminal A | 225363_at | PTEN | -1.2859 | -1.2859 | 1.20e-14 | 6.07e-11 | Yes | Both |
| TNBC - Luminal A | 228006_at | PTEN | -1.7514 | -1.7514 | 2.83e-12 | 3.38e-09 | Yes | Both |
| TNBC - Luminal A | 233254_x_at | PTEN | -0.7759 | -0.7759 | 5.11e-07 | 1.10e-05 | Yes | Both |
| TNBC - Luminal A | 242622_x_at | PTEN | -0.1847 | -0.1847 | 0.0704 | 0.07444 | Yes | Neither |
| TNBC - Luminal A | 1565483_at | EGFR | -0.0753 | -0.0753 | 0.1564 | 0.4032 | Yes | Neither |
| TNBC - Luminal A | 1565484_x_at | EGFR | 0.0004 | 0.0004 | 0.9979 | 0.9955 | Yes | Neither |
| TNBC - Luminal A | 201983_s_at | EGFR | 1.8761 | 1.8761 | 5.73e-05 | 1.29e-05 | Yes | Both |
| TNBC - Luminal A | 201984_s_at | EGFR | 2.1202 | 2.1202 | 1.12e-09 | 5.89e-11 | Yes | Both |
| TNBC - Luminal A | 210984_x_at | EGFR | 0.0668 | 0.0668 | 0.7891 | 0.4867 | Yes | Neither |
| TNBC - Luminal A | 211550_at | EGFR | 0.0449 | 0.0449 | 0.0293 | 0.07109 | Yes | Limma only |
| TNBC - Luminal A | 211551_at | EGFR | 0.0110 | 0.0110 | 0.008894 | 0.01686 | Yes | Both |
| TNBC - Luminal A | 211607_x_at | EGFR | 0.0715 | 0.0715 | 0.7507 | 0.4867 | Yes | Neither |
| TNBC - Luminal A | 224999_at | EGFR | 1.2600 | 1.2600 | 0.0009666 | 0.0002456 | Yes | Both |
| TNBC - Luminal B | 204531_s_at | BRCA1 | 0.1072 | 0.1072 | 0.7713 | 0.7328 | Yes | Neither |
| TNBC - Luminal B | 211851_x_at | BRCA1 | 0.2291 | 0.2291 | 0.07507 | 0.1528 | Yes | Neither |
| TNBC - Luminal B | 208368_s_at | BRCA2 | 0.6374 | 0.6374 | 0.003178 | 0.005929 | Yes | Both |
| TNBC - Luminal B | 214727_at | BRCA2 | 0.2863 | 0.2863 | 0.1018 | 0.1201 | Yes | Neither |
| TNBC - Luminal B | 201746_at | TP53 | 0.5886 | 0.5886 | 0.07507 | 0.06872 | Yes | Neither |
| TNBC - Luminal B | 211300_s_at | TP53 | 1.1735 | 1.1735 | 0.0004666 | 0.0007236 | Yes | Both |
| TNBC - Luminal B | 210930_s_at | ERBB2 | -3.3937 | -3.3937 | 5.67e-18 | 3.62e-08 | Yes | Both |
| TNBC - Luminal B | 216836_s_at | ERBB2 | -4.2681 | -4.2681 | 1.19e-31 | 5.14e-17 | Yes | Both |
| TNBC - Luminal B | 234354_x_at | ERBB2 | -1.8880 | -1.8880 | 1.97e-08 | 6.40e-05 | Yes | Both |
| TNBC - Luminal B | 204369_at | PIK3CA | 0.6338 | 0.6338 | 0.005284 | 0.0161 | Yes | Both |
| TNBC - Luminal B | 231854_at | PIK3CA | 0.0378 | 0.0378 | 0.6623 | 0.6625 | Yes | Neither |
| TNBC - Luminal B | 235980_at | PIK3CA | 0.0008 | 0.0008 | 0.9979 | 0.9955 | Yes | Neither |
| TNBC - Luminal B | 205225_at | ESR1 | -7.9556 | -7.9556 | 5.39e-37 | 2.21e-27 | Yes | Both |
| TNBC - Luminal B | 211233_x_at | ESR1 | -0.2157 | -0.2157 | 0.08825 | 0.0676 | Yes | Neither |
| TNBC - Luminal B | 211234_x_at | ESR1 | -0.2931 | -0.2931 | 0.01071 | 0.007954 | Yes | Both |
| TNBC - Luminal B | 211235_s_at | ESR1 | -0.7444 | -0.7444 | 0.0003273 | 0.0008979 | Yes | Both |
| TNBC - Luminal B | 211627_x_at | ESR1 | -0.0192 | -0.0192 | 0.8163 | 0.7044 | Yes | Neither |
| TNBC - Luminal B | 215551_at | ESR1 | -0.3215 | -0.3215 | 0.005284 | 0.009071 | Yes | Both |
| TNBC - Luminal B | 215552_s_at | ESR1 | -1.9061 | -1.9061 | 5.00e-10 | 1.49e-06 | Yes | Both |
| TNBC - Luminal B | 217163_at | ESR1 | 0.0024 | 0.0024 | 0.567 | 0.4525 | Yes | Neither |
| TNBC - Luminal B | 217190_x_at | ESR1 | -0.1063 | -0.1063 | 0.1751 | 0.04448 | Yes | Welch only |
| TNBC - Luminal B | 204053_x_at | PTEN | -0.6067 | -0.6067 | 0.004527 | 0.0143 | Yes | Both |
| TNBC - Luminal B | 204054_at | PTEN | -0.2870 | -0.2870 | 0.1843 | 0.2274 | Yes | Neither |
| TNBC - Luminal B | 211711_s_at | PTEN | -0.5378 | -0.5378 | 0.0358 | 0.0676 | Yes | Limma only |
| TNBC - Luminal B | 217492_s_at | PTEN | -0.5464 | -0.5464 | 0.01024 | 0.02763 | Yes | Both |
| TNBC - Luminal B | 225363_at | PTEN | -1.0192 | -1.0192 | 1.75e-10 | 4.97e-08 | Yes | Both |
| TNBC - Luminal B | 228006_at | PTEN | -1.3833 | -1.3833 | 1.13e-08 | 1.14e-06 | Yes | Both |
| TNBC - Luminal B | 233254_x_at | PTEN | -0.6431 | -0.6431 | 2.80e-05 | 9.69e-05 | Yes | Both |
| TNBC - Luminal B | 242622_x_at | PTEN | -0.2093 | -0.2093 | 0.03592 | 0.03363 | Yes | Both |
| TNBC - Luminal B | 1565483_at | EGFR | -0.0217 | -0.0217 | 0.7617 | 0.04046 | Yes | Welch only |
| TNBC - Luminal B | 1565484_x_at | EGFR | 0.0001 | 0.0001 | 0.9979 | 0.9955 | Yes | Neither |
| TNBC - Luminal B | 201983_s_at | EGFR | 3.4230 | 3.4230 | 7.61e-13 | 1.55e-10 | Yes | Both |
| TNBC - Luminal B | 201984_s_at | EGFR | 2.2629 | 2.2629 | 7.64e-11 | 8.43e-12 | Yes | Both |
| TNBC - Luminal B | 210984_x_at | EGFR | 0.0648 | 0.0648 | 0.79 | 0.506 | Yes | Neither |
| TNBC - Luminal B | 211550_at | EGFR | 0.0123 | 0.0123 | 0.627 | 0.6083 | Yes | Neither |
| TNBC - Luminal B | 211551_at | EGFR | 0.0031 | 0.0031 | 0.5425 | 0.4525 | Yes | Neither |
| TNBC - Luminal B | 211607_x_at | EGFR | 0.0534 | 0.0534 | 0.8203 | 0.6115 | Yes | Neither |
| TNBC - Luminal B | 224999_at | EGFR | 2.3262 | 2.3262 | 1.00e-09 | 3.41e-10 | Yes | Both |
| Her2 - Luminal A | 204531_s_at | BRCA1 | 0.6931 | 0.6931 | 0.01291 | 0.002773 | Yes | Both |
| Her2 - Luminal A | 211851_x_at | BRCA1 | 0.0859 | 0.0859 | 0.6185 | 0.1603 | Yes | Neither |
| Her2 - Luminal A | 208368_s_at | BRCA2 | 1.0050 | 1.0050 | 1.30e-05 | 1.22e-05 | Yes | Both |
| Her2 - Luminal A | 214727_at | BRCA2 | 0.4049 | 0.4049 | 0.03006 | 0.01122 | Yes | Both |
| Her2 - Luminal A | 201746_at | TP53 | -0.4762 | -0.4762 | 0.2007 | 0.1676 | Yes | Neither |
| Her2 - Luminal A | 211300_s_at | TP53 | -0.1000 | -0.1000 | 0.8535 | 0.7916 | Yes | Neither |
| Her2 - Luminal A | 210930_s_at | ERBB2 | 5.9129 | 5.9129 | 4.81e-33 | 7.04e-27 | Yes | Both |
| Her2 - Luminal A | 216836_s_at | ERBB2 | 4.1477 | 4.1477 | 9.45e-28 | 1.32e-26 | Yes | Both |
| Her2 - Luminal A | 234354_x_at | ERBB2 | 3.5406 | 3.5406 | 4.05e-19 | 1.27e-11 | Yes | Both |
| Her2 - Luminal A | 204369_at | PIK3CA | -0.2587 | -0.2587 | 0.3539 | 0.2877 | Yes | Neither |
| Her2 - Luminal A | 231854_at | PIK3CA | -0.0770 | -0.0770 | 0.3539 | 0.1147 | Yes | Neither |
| Her2 - Luminal A | 235980_at | PIK3CA | -0.3549 | -0.3549 | 0.02349 | 0.004559 | Yes | Both |
| Her2 - Luminal A | 205225_at | ESR1 | -7.8638 | -7.8638 | 2.00e-33 | 4.11e-23 | Yes | Both |
| Her2 - Luminal A | 211233_x_at | ESR1 | -0.4528 | -0.4528 | 0.0005512 | 0.01132 | Yes | Both |
| Her2 - Luminal A | 211234_x_at | ESR1 | -0.4671 | -0.4671 | 0.0001414 | 0.005929 | Yes | Both |
| Her2 - Luminal A | 211235_s_at | ESR1 | -1.1849 | -1.1849 | 1.10e-07 | 0.000162 | Yes | Both |
| Her2 - Luminal A | 211627_x_at | ESR1 | -0.0539 | -0.0539 | 0.4458 | 0.5312 | Yes | Neither |
| Her2 - Luminal A | 215551_at | ESR1 | -0.2442 | -0.2442 | 0.05608 | 0.1222 | Yes | Neither |
| Her2 - Luminal A | 215552_s_at | ESR1 | -2.2285 | -2.2285 | 2.92e-11 | 1.70e-06 | Yes | Both |
| Her2 - Luminal A | 217163_at | ESR1 | 0.0057 | 0.0057 | 0.127 | 0.07919 | Yes | Neither |
| Her2 - Luminal A | 217190_x_at | ESR1 | -0.2067 | -0.2067 | 0.01024 | 0.07426 | Yes | Limma only |
| Her2 - Luminal A | 204053_x_at | PTEN | 0.2916 | 0.2916 | 0.2469 | 0.02193 | Yes | Welch only |
| Her2 - Luminal A | 204054_at | PTEN | -0.0737 | -0.0737 | 0.8214 | 0.678 | Yes | Neither |
| Her2 - Luminal A | 211711_s_at | PTEN | 0.0122 | 0.0122 | 0.9979 | 0.9793 | Yes | Neither |
| Her2 - Luminal A | 217492_s_at | PTEN | 0.2259 | 0.2259 | 0.3878 | 0.09481 | Yes | Neither |
| Her2 - Luminal A | 225363_at | PTEN | -0.5841 | -0.5841 | 0.0006298 | 5.37e-06 | Yes | Both |
| Her2 - Luminal A | 228006_at | PTEN | -0.9112 | -0.9112 | 0.0005336 | 2.23e-06 | Yes | Both |
| Her2 - Luminal A | 233254_x_at | PTEN | -0.2197 | -0.2197 | 0.228 | 0.1354 | Yes | Neither |
| Her2 - Luminal A | 242622_x_at | PTEN | 0.0416 | 0.0416 | 0.7891 | 0.7146 | Yes | Neither |
| Her2 - Luminal A | 1565483_at | EGFR | -0.0541 | -0.0541 | 0.3835 | 0.5535 | Yes | Neither |
| Her2 - Luminal A | 1565484_x_at | EGFR | 0.0055 | 0.0055 | 0.8535 | 0.8737 | Yes | Neither |
| Her2 - Luminal A | 201983_s_at | EGFR | 0.8143 | 0.8143 | 0.1253 | 0.1861 | Yes | Neither |
| Her2 - Luminal A | 201984_s_at | EGFR | 1.3879 | 1.3879 | 0.0001944 | 0.002086 | Yes | Both |
| Her2 - Luminal A | 210984_x_at | EGFR | 0.2110 | 0.2110 | 0.3003 | 0.4471 | Yes | Neither |
| Her2 - Luminal A | 211550_at | EGFR | 0.0388 | 0.0388 | 0.08179 | 0.1056 | Yes | Neither |
| Her2 - Luminal A | 211551_at | EGFR | 0.0116 | 0.0116 | 0.01042 | 0.01074 | Yes | Both |
| Her2 - Luminal A | 211607_x_at | EGFR | 0.1840 | 0.1840 | 0.3442 | 0.4543 | Yes | Neither |
| Her2 - Luminal A | 224999_at | EGFR | 0.4466 | 0.4466 | 0.3442 | 0.4077 | Yes | Neither |
| Her2 - Luminal B | 204531_s_at | BRCA1 | -0.0155 | -0.0155 | 0.9979 | 0.9793 | Yes | Neither |
| Her2 - Luminal B | 211851_x_at | BRCA1 | -0.0041 | -0.0041 | 0.9979 | 0.9912 | Yes | Neither |
| Her2 - Luminal B | 208368_s_at | BRCA2 | 0.4067 | 0.4067 | 0.08912 | 0.1113 | Yes | Neither |
| Her2 - Luminal B | 214727_at | BRCA2 | 0.2428 | 0.2428 | 0.2134 | 0.1959 | Yes | Neither |
| Her2 - Luminal B | 201746_at | TP53 | -0.0939 | -0.0939 | 0.8535 | 0.8186 | Yes | Neither |
| Her2 - Luminal B | 211300_s_at | TP53 | 0.4252 | 0.4252 | 0.2984 | 0.1999 | Yes | Neither |
| Her2 - Luminal B | 210930_s_at | ERBB2 | 2.8457 | 2.8457 | 1.87e-12 | 1.61e-06 | Yes | Both |
| Her2 - Luminal B | 216836_s_at | ERBB2 | 0.9224 | 0.9224 | 0.004527 | 0.009194 | Yes | Both |
| Her2 - Luminal B | 234354_x_at | ERBB2 | 1.6519 | 1.6519 | 4.64e-06 | 0.00305 | Yes | Both |
| Her2 - Luminal B | 204369_at | PIK3CA | 0.0170 | 0.0170 | 0.9979 | 0.9829 | Yes | Neither |
| Her2 - Luminal B | 231854_at | PIK3CA | -0.1319 | -0.1319 | 0.07673 | 0.01202 | Yes | Welch only |
| Her2 - Luminal B | 235980_at | PIK3CA | -0.3145 | -0.3145 | 0.04337 | 0.03377 | Yes | Both |
| Her2 - Luminal B | 205225_at | ESR1 | -7.6088 | -7.6088 | 1.24e-32 | 1.38e-22 | Yes | Both |
| Her2 - Luminal B | 211233_x_at | ESR1 | -0.2147 | -0.2147 | 0.1179 | 0.0676 | Yes | Neither |
| Her2 - Luminal B | 211234_x_at | ESR1 | -0.2907 | -0.2907 | 0.01956 | 0.008176 | Yes | Both |
| Her2 - Luminal B | 211235_s_at | ESR1 | -0.7407 | -0.7407 | 0.0008592 | 0.0009363 | Yes | Both |
| Her2 - Luminal B | 211627_x_at | ESR1 | -0.0178 | -0.0178 | 0.848 | 0.7283 | Yes | Neither |
| Her2 - Luminal B | 215551_at | ESR1 | -0.3275 | -0.3275 | 0.008753 | 0.008002 | Yes | Both |
| Her2 - Luminal B | 215552_s_at | ESR1 | -1.8898 | -1.8898 | 7.63e-09 | 1.61e-06 | Yes | Both |
| Her2 - Luminal B | 217163_at | ESR1 | 0.0032 | 0.0032 | 0.4415 | 0.3272 | Yes | Neither |
| Her2 - Luminal B | 217190_x_at | ESR1 | -0.0675 | -0.0675 | 0.4712 | 0.286 | Yes | Neither |
| Her2 - Luminal B | 204053_x_at | PTEN | 0.1335 | 0.1335 | 0.645 | 0.4089 | Yes | Neither |
| Her2 - Luminal B | 204054_at | PTEN | 0.0602 | 0.0602 | 0.8535 | 0.7545 | Yes | Neither |
| Her2 - Luminal B | 211711_s_at | PTEN | 0.1834 | 0.1834 | 0.5896 | 0.4142 | Yes | Neither |
| Her2 - Luminal B | 217492_s_at | PTEN | 0.1392 | 0.1392 | 0.627 | 0.4221 | Yes | Neither |
| Her2 - Luminal B | 225363_at | PTEN | -0.3173 | -0.3173 | 0.07321 | 0.01193 | Yes | Welch only |
| Her2 - Luminal B | 228006_at | PTEN | -0.5431 | -0.5431 | 0.04327 | 0.003597 | Yes | Both |
| Her2 - Luminal B | 233254_x_at | PTEN | -0.0870 | -0.0870 | 0.6943 | 0.5442 | Yes | Neither |
| Her2 - Luminal B | 242622_x_at | PTEN | 0.0169 | 0.0169 | 0.9304 | 0.8796 | Yes | Neither |
| Her2 - Luminal B | 1565483_at | EGFR | -0.0005 | -0.0005 | 0.9979 | 0.9955 | Yes | Neither |
| Her2 - Luminal B | 1565484_x_at | EGFR | 0.0052 | 0.0052 | 0.8589 | 0.8737 | Yes | Neither |
| Her2 - Luminal B | 201983_s_at | EGFR | 2.3612 | 2.3612 | 1.82e-06 | 0.0001868 | Yes | Both |
| Her2 - Luminal B | 201984_s_at | EGFR | 1.5307 | 1.5307 | 3.04e-05 | 0.0007397 | Yes | Both |
| Her2 - Luminal B | 210984_x_at | EGFR | 0.2091 | 0.2091 | 0.3003 | 0.4503 | Yes | Neither |
| Her2 - Luminal B | 211550_at | EGFR | 0.0062 | 0.0062 | 0.8535 | 0.7885 | Yes | Neither |
| Her2 - Luminal B | 211551_at | EGFR | 0.0037 | 0.0037 | 0.4898 | 0.3579 | Yes | Neither |
| Her2 - Luminal B | 211607_x_at | EGFR | 0.1659 | 0.1659 | 0.3881 | 0.5068 | Yes | Neither |
| Her2 - Luminal B | 224999_at | EGFR | 1.5129 | 1.5129 | 0.0002021 | 0.002268 | Yes | Both |
| Luminal A - Luminal B | 204531_s_at | BRCA1 | -0.7086 | -0.7086 | 0.01104 | 0.00141 | Yes | Both |
| Luminal A - Luminal B | 211851_x_at | BRCA1 | -0.0900 | -0.0900 | 0.5972 | 0.4341 | Yes | Neither |
| Luminal A - Luminal B | 208368_s_at | BRCA2 | -0.5984 | -0.5984 | 0.01088 | 0.00305 | Yes | Both |
| Luminal A - Luminal B | 214727_at | BRCA2 | -0.1621 | -0.1621 | 0.4558 | 0.1299 | Yes | Neither |
| Luminal A - Luminal B | 201746_at | TP53 | 0.3823 | 0.3823 | 0.3344 | 0.09481 | Yes | Neither |
| Luminal A - Luminal B | 211300_s_at | TP53 | 0.5251 | 0.5251 | 0.1789 | 0.04357 | Yes | Welch only |
| Luminal A - Luminal B | 210930_s_at | ERBB2 | -3.0672 | -3.0672 | 9.86e-14 | 3.12e-07 | Yes | Both |
| Luminal A - Luminal B | 216836_s_at | ERBB2 | -3.2254 | -3.2254 | 6.58e-20 | 2.38e-12 | Yes | Both |
| Luminal A - Luminal B | 234354_x_at | ERBB2 | -1.8888 | -1.8888 | 2.04e-07 | 6.40e-05 | Yes | Both |
| Luminal A - Luminal B | 204369_at | PIK3CA | 0.2757 | 0.2757 | 0.3218 | 0.3272 | Yes | Neither |
| Luminal A - Luminal B | 231854_at | PIK3CA | -0.0549 | -0.0549 | 0.545 | 0.07405 | Yes | Neither |
| Luminal A - Luminal B | 235980_at | PIK3CA | 0.0404 | 0.0404 | 0.8535 | 0.7639 | Yes | Neither |
| Luminal A - Luminal B | 205225_at | ESR1 | 0.2550 | 0.2550 | 0.7262 | 0.5974 | Yes | Neither |
| Luminal A - Luminal B | 211233_x_at | ESR1 | 0.2381 | 0.2381 | 0.08179 | 0.2753 | Yes | Neither |
| Luminal A - Luminal B | 211234_x_at | ESR1 | 0.1765 | 0.1765 | 0.1809 | 0.4077 | Yes | Neither |
| Luminal A - Luminal B | 211235_s_at | ESR1 | 0.4443 | 0.4443 | 0.05608 | 0.2242 | Yes | Neither |
| Luminal A - Luminal B | 211627_x_at | ESR1 | 0.0361 | 0.0361 | 0.63 | 0.7044 | Yes | Neither |
| Luminal A - Luminal B | 215551_at | ESR1 | -0.0833 | -0.0833 | 0.5972 | 0.7044 | Yes | Neither |
| Luminal A - Luminal B | 215552_s_at | ESR1 | 0.3387 | 0.3387 | 0.3761 | 0.5458 | Yes | Neither |
| Luminal A - Luminal B | 217163_at | ESR1 | -0.0025 | -0.0025 | 0.5841 | 0.5312 | Yes | Neither |
| Luminal A - Luminal B | 217190_x_at | ESR1 | 0.1391 | 0.1391 | 0.09007 | 0.2753 | Yes | Neither |
| Luminal A - Luminal B | 204053_x_at | PTEN | -0.1581 | -0.1581 | 0.5896 | 0.2653 | Yes | Neither |
| Luminal A - Luminal B | 204054_at | PTEN | 0.1339 | 0.1339 | 0.6285 | 0.3862 | Yes | Neither |
| Luminal A - Luminal B | 211711_s_at | PTEN | 0.1713 | 0.1713 | 0.6188 | 0.4429 | Yes | Neither |
| Luminal A - Luminal B | 217492_s_at | PTEN | -0.0866 | -0.0866 | 0.79 | 0.6083 | Yes | Neither |
| Luminal A - Luminal B | 225363_at | PTEN | 0.2668 | 0.2668 | 0.1429 | 0.01558 | Yes | Welch only |
| Luminal A - Luminal B | 228006_at | PTEN | 0.3681 | 0.3681 | 0.1946 | 0.01846 | Yes | Welch only |
| Luminal A - Luminal B | 233254_x_at | PTEN | 0.1327 | 0.1327 | 0.5219 | 0.3517 | Yes | Neither |
| Luminal A - Luminal B | 242622_x_at | PTEN | -0.0247 | -0.0247 | 0.8636 | 0.7864 | Yes | Neither |
| Luminal A - Luminal B | 1565483_at | EGFR | 0.0536 | 0.0536 | 0.3868 | 0.5442 | Yes | Neither |
| Luminal A - Luminal B | 1565484_x_at | EGFR | -0.0003 | -0.0003 | 0.9979 | 0.9955 | Yes | Neither |
| Luminal A - Luminal B | 201983_s_at | EGFR | 1.5468 | 1.5468 | 0.002203 | 0.00413 | Yes | Both |
| Luminal A - Luminal B | 201984_s_at | EGFR | 0.1427 | 0.1427 | 0.79 | 0.4341 | Yes | Neither |
| Luminal A - Luminal B | 210984_x_at | EGFR | -0.0020 | -0.0020 | 0.9979 | 0.989 | Yes | Neither |
| Luminal A - Luminal B | 211550_at | EGFR | -0.0326 | -0.0326 | 0.1548 | 0.2897 | Yes | Neither |
| Luminal A - Luminal B | 211551_at | EGFR | -0.0079 | -0.0079 | 0.08838 | 0.1321 | Yes | Neither |
| Luminal A - Luminal B | 211607_x_at | EGFR | -0.0181 | -0.0181 | 0.9749 | 0.6817 | Yes | Neither |
| Luminal A - Luminal B | 224999_at | EGFR | 1.0662 | 1.0662 | 0.01053 | 0.003887 | Yes | Both |

## 8. Biological interpretation

The panel spans genes commonly studied in breast-cancer biology, including DNA-damage response (BRCA1/BRCA2), tumour-suppressor and cell-cycle biology (TP53/PTEN), growth-factor signalling (ERBB2/EGFR), PI3K-pathway signalling (PIK3CA), and estrogen-related signalling (ESR1). The observed probe-level signals can motivate follow-up questions about tumour/healthy differences and subtype-associated expression. Agreement across exact probes in an independent cohort is supportive of directional reproducibility for those probes, while discordance and probe-to-probe heterogeneity caution against broad gene-level summaries.

These are associations in retrospective public expression datasets. They do not establish causal mechanisms, protein abundance, pathway activity, treatment response, patient-level receptor status or diagnostic performance. No literature references or clinical annotations beyond the inspected project inputs are asserted here. Subtype labels are used as supplied; ER-positive status is not inferred solely from Luminal A or Luminal B labels.

## 9. Limitations

- **One independent cohort:** GSE42568 is the only independent validation cohort included.
- **Targeted panel:** validation covers 38 probes selected for eight genes, not genome-wide validation of those genes or the transcriptome.
- **Probe ambiguity:** `217492_s_at` is annotated `PTENP1 /// PTEN`; it qualifies under the target-panel annotation but cannot be attributed uniquely to PTEN.
- **Preprocessing differences:** both cohorts report GC-RMA-related processing, but GSE65194 additionally documents batch/hybridization fixed-effect correction. Exact CDF and explicit transformation history are unresolved, and the cohorts' scales must not be assumed identical.
- **Sample sizes and composition:** GSE42568 has 104 cancer and 17 normal samples; GSE65194 has 130 tumour patients and only 11 healthy controls. The smaller healthy group limits precision.
- **Method differences:** the GSE42568 source is GEO2R/limma, GSE65194 tumour-versus-healthy includes Welch and local limma analyses, and subtype comparisons use Welch and limma on patient-level data.
- **Different testing families:** GSE65194 tumour/healthy Welch uses BH over 38 probes; its limma fits use 54,673 probes; subtype comparisons use 228 tests per method; original cross-cohort adjusted p-values come from different source families. The matched-panel BH output adjusts each cohort's raw values across 38 exact probes.
- **Technical replicate provenance:** the GEO processing note says technical replicates were averaged, but the deposited array matrix contains separate columns for all 23 pairs. The project's probe-wise patient average is derived from these array-level values; the earlier processing history is unresolved. repA/repB sensitivity selections are not independent cohorts.
- **No clinical or laboratory validation:** the results are not diagnostic, prognostic, or treatment guidance and have not been validated experimentally.
- **Clinical receptor status:** no explicit ER-positive classification is inferred from subtype names; no unsupported clinical status is assigned.
- **Exploratory design:** the panel and analyses do not demonstrate causation, clinical utility, or independent prospective performance.

## 10. Conclusion

This analysis describes exact-probe tumour/healthy and subtype patterns in GSE65194 and compares tumour/normal directions with an independent GPL570 cohort. The outputs show 24/38 cross-cohort direction concordance, 228/228 subtype-direction agreement between Welch and limma, and clear probe-level heterogeneity that warrants preserving probe identities. These findings support exploratory follow-up only. Replication in additional cohorts, clarification of processing provenance, clinical annotation, and laboratory validation remain necessary before making biological or clinical claims.

## 11. Reproducibility

### Inputs and software

- GEO series: [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568) and [GSE65194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65194).
- GSE65194 input: `data/raw/GSE65194_series_matrix.txt.gz`; derived inputs: `data/processed/GSE65194_patient_level_matrix.tsv.gz` and `data/processed/GSE65194_patient_level_manifest.tsv`.
- Existing cohort-level GSE42568 GEO2R/limma output and report are external to this v2 workspace; the v2 replication and matched-panel TSVs retain the probe-level values used here.
- Software recorded in source reports: R 4.6.2, limma 3.68.5.; the final builder uses the existing `.venv` Python 3.15.2 and pandas 3.1.0. No packages were installed and no scientific analysis was rerun for this report/dashboard build.

### Analysis scripts and commands

- `scripts/01_load_and_qc.py` — initial matrix/metadata QC.
- `scripts/02_build_patient_level_matrix.py` — derived patient-level aggregation.
- `scripts/03_differential_expression.py` — targeted Welch analyses and report.
- `scripts/04_limma_sensitivity_analysis.R` — alternative repA/repB full-probe limma fits.
- `scripts/05_limma_subtype_analysis.R` — patient-level limma subtype contrasts.
- Builder command: `\.venv\Scripts\python.exe scripts\06_build_final_report_dashboard.py`
- Limma commands previously used: `Rscript scripts/04_limma_sensitivity_analysis.R` and `Rscript scripts/05_limma_subtype_analysis.R` with R 4.6.2 and limma 3.68.5.

### Existing result files used

[`GSE42568_vs_GSE65194_replication.tsv`](GSE42568_vs_GSE65194_replication.tsv), [`GSE42568_vs_GSE65194_matched_panel_BH.tsv`](GSE42568_vs_GSE65194_matched_panel_BH.tsv), [`GSE65194_tumor_vs_healthy_welch.tsv`](GSE65194_tumor_vs_healthy_welch.tsv), [`GSE65194_limma_repA_vs_repB.tsv`](GSE65194_limma_repA_vs_repB.tsv), [`GSE65194_subtype_pairwise_welch.tsv`](GSE65194_subtype_pairwise_welch.tsv), [`GSE65194_limma_subtype_pairwise.tsv`](GSE65194_limma_subtype_pairwise.tsv), [`GSE65194_limma_repA_full_table.tsv`](../data/processed/GSE65194_limma_repA_full_table.tsv), [`GSE65194_limma_repB_full_table.tsv`](../data/processed/GSE65194_limma_repB_full_table.tsv)

The standalone dashboard is [`../dashboard/index.html`](../dashboard/index.html). Its visualization data are embedded from the validated result tables; it uses no external scripts, fonts, APIs, CDNs or runtime data fetches.
