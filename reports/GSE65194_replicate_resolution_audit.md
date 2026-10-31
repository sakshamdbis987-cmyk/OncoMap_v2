# GSE65194 replicate-resolution audit

**Scope:** Read-only verification of the GEO series matrix, sample manifest,
compatibility audit, clinical metadata, and prior QC outputs. No expression
values or source files were changed. No replicate arrays were averaged or
removed, and no differential-expression analysis was performed.

## Verified facts

### Replicate pair validation

The matrix contains 178 GSM columns. The manifest identifies 153 tumour arrays
from 130 patients: 23 patients have two arrays and 107 patients have one.
For each of the 23 pairs below, the two GSMs map to the same `tumor_id`; both
matrix titles agree exactly with the corresponding manifest `sample_title`,
and their `_repA`/`_repB` title suffixes agree with the manifest replicate
labels. The listed Pearson correlations are from the existing QC replicate
metrics, calculated across all 54,673 probes.

| Patient ID | Array A GSM | Array A matrix title | Label | Array B GSM | Array B matrix title | Label | Sample group | Pearson correlation |
|---|---|---|---|---|---|---|---|---:|
| TUM001 | GSM1588970 | TNBC, TUM001_repA | repA | GSM1588971 | TNBC, TUM001_repB | repB | TNBC | 0.99386507 |
| TUM002 | GSM1588972 | Her2, TUM002_repA | repA | GSM1588973 | Her2, TUM002_repB | repB | Her2 | 0.99680579 |
| TUM006 | GSM1588977 | TNBC, TUM006_repA | repA | GSM1588978 | TNBC, TUM006_repB | repB | TNBC | 0.99635015 |
| TUM010 | GSM1588982 | TNBC, TUM010_repA | repA | GSM1588983 | TNBC, TUM010_repB | repB | TNBC | 0.99402480 |
| TUM013 | GSM1588986 | TNBC, TUM013_repA | repA | GSM1588987 | TNBC, TUM013_repB | repB | TNBC | 0.99396011 |
| TUM015 | GSM1588989 | Her2, TUM015_repA | repA | GSM1588990 | Her2, TUM015_repB | repB | Her2 | 0.99574838 |
| TUM016 | GSM1588991 | Her2, TUM016_repA | repA | GSM1588992 | Her2, TUM016_repB | repB | Her2 | 0.99680104 |
| TUM017 | GSM1588993 | TNBC, TUM017_repA | repA | GSM1588994 | TNBC, TUM017_repB | repB | TNBC | 0.99671069 |
| TUM018 | GSM1588995 | TNBC, TUM018_repA | repA | GSM1588996 | TNBC, TUM018_repB | repB | TNBC | 0.99542977 |
| TUM020 | GSM1588998 | Her2, TUM020_repA | repA | GSM1588999 | Her2, TUM020_repB | repB | Her2 | 0.99695276 |
| TUM021 | GSM1589000 | Her2, TUM021_repA | repA | GSM1589001 | Her2, TUM021_repB | repB | Her2 | 0.99654635 |
| TUM026 | GSM1589006 | Her2, TUM026_repA | repA | GSM1589007 | Her2, TUM026_repB | repB | Her2 | 0.99640705 |
| TUM028 | GSM1589009 | TNBC, TUM028_repA | repA | GSM1589010 | TNBC, TUM028_repB | repB | TNBC | 0.99595858 |
| TUM032 | GSM1589014 | TNBC, TUM032_repA | repA | GSM1589015 | TNBC, TUM032_repB | repB | TNBC | 0.99483256 |
| TUM035 | GSM1589018 | Her2, TUM035_repA | repA | GSM1589019 | Her2, TUM035_repB | repB | Her2 | 0.99672584 |
| TUM038 | GSM1589022 | Her2, TUM038_repA | repA | GSM1589023 | Her2, TUM038_repB | repB | Her2 | 0.99611082 |
| TUM043 | GSM1589028 | TNBC, TUM043_repA | repA | GSM1589029 | TNBC, TUM043_repB | repB | TNBC | 0.99553655 |
| TUM045 | GSM1589031 | TNBC, TUM045_repA | repA | GSM1589032 | TNBC, TUM045_repB | repB | TNBC | 0.99635094 |
| TUM048 | GSM1589035 | TNBC, TUM048_repA | repA | GSM1589036 | TNBC, TUM048_repB | repB | TNBC | 0.99667885 |
| TUM051 | GSM1589039 | TNBC, TUM051_repA | repA | GSM1589040 | TNBC, TUM051_repB | repB | TNBC | 0.99622609 |
| TUM058 | GSM1589047 | Her2, TUM058_repA | repA | GSM1589048 | Her2, TUM058_repB | repB | Her2 | 0.99657850 |
| TUM064 | GSM1589054 | TNBC, TUM064_repA | repA | GSM1589055 | TNBC, TUM064_repB | repB | TNBC | 0.99405899 |
| TUM067 | GSM1589058 | TNBC, TUM067_repA | repA | GSM1589059 | TNBC, TUM067_repB | repB | TNBC | 0.99572814 |

There are exactly **23** validated pairs. Their correlations range from
**0.99386507 to 0.99695276**. These correlations indicate strong similarity
between paired array profiles, but do not by themselves establish whether or
how those arrays were averaged during an earlier processing step.

### Pooled-median-IQR flagged tumour arrays

The existing outlier review identifies two flagged tumour arrays. The screen
uses the pooled distribution of all 178 array medians and the 1.5×IQR rule;
both flags are above the pooled upper fence.

| Flagged GSM | Patient | Title / group | Replicate label | Partner | Pair correlation | Interpretation |
|---|---|---|---|---|---:|---|
| GSM1588986 | TUM013 | TNBC, TUM013_repA / TNBC | repA | GSM1588987 (repB) | 0.99396011 | The partner is also in the pair metrics and is highly correlated. The median flag alone does not establish a technical problem. |
| GSM1589121 | TUM123 | Luminal B, TUM123 / Luminal B | None | No second array for TUM123 | Not available | There is no replicate correlation for this single-array patient; inspect with batch and other technical metadata if warranted. |

Neither is classified as failed. The pooled screen combines healthy, cell-line,
and tumour samples, so a flag can reflect group-level distribution differences
instead of a technical failure.

### Patient counts by subtype after duplicate arrays are accounted for

The manifest has 55 TNBC, 39 Her2, 29 Luminal A, and 30 Luminal B tumour arrays.
The 23 duplicate arrays comprise 14 TNBC pairs and 9 Her2 pairs; Luminal A and
Luminal B have no duplicated patients. Counting each `tumor_id` once gives:

| Preserved manifest sample group | Arrays | Duplicate pairs | Unique tumour patients |
|---|---:|---:|---:|
| TNBC (Basal-like terminology in clinical metadata) | 55 | 14 | 41 |
| Her2 | 39 | 9 | 30 |
| Luminal A | 29 | 0 | 29 |
| Luminal B | 30 | 0 | 30 |
| **Total tumour cohort** | **153** | **23** | **130** |

This agrees with the GEO series summary, which describes 130 breast cancer
samples as 41 TNBC, 30 Her2, 30 Luminal B, and 29 Luminal A. The clinical file
contains 130 GSM records mapping to 130 unique tumour IDs. Its subtype labels
are Basal-like (41), Her2 + (30), luminal A (29), and luminal B (30), consistent
with the preserved matrix/manifest categories after the vocabulary difference
for TNBC/Basal-like and Her2/HER2+ is noted.

## GEO replicate-averaging statement and remaining discrepancy

**Verified GEO processing text:** the `!Sample_data_processing` annotation says:

> Samples with technical replicates were subsequently averaged.

The same processing description mentions GC-RMA normalization and correction
for batch and hybridization effects. Despite that statement, the downloadable
series matrix contains 178 GSM columns, including both separate GSM accessions
for all 23 pairs. The GEO overall design also describes 178 arrays, of which
153 arrays represent 130 tumour samples plus 23 technical duplicates, alongside
11 normal samples and 14 cell lines.

**What is established:** the downloaded matrix is array-level in its columns;
it has not collapsed the 23 replicate pairs to 130 tumour columns. It therefore
retains 178 sample columns, not one column per biological tumour patient.

**Unresolved:** the available description does not establish exactly where or
whether replicate averaging occurred. It could describe a downstream analysis
table distinct from the deposited series matrix, or the processing annotation
could be stale/inconsistently applied. Neither possibility is verified here.
High pairwise correlations cannot resolve this provenance question. Until
clarified, do not treat the matrix columns as 178 independent tumour patients
or assume that replicate averaging has already occurred.

## Recommendation for one observation per patient

For downstream patient-level comparisons, the biological unit should be the
patient (`tumor_id`), not the array. First seek confirmation of whether the
deposited values were already replicate-averaged and verify that the paired
arrays are technical replicates on the same processing scale.

If these are confirmed to be separate measurements of the same patient and
the deposited expression values have not already been combined, a defensible
approach is to summarize each probe's replicate measurements into **one
patient-level value** (for example, the arithmetic mean of the two values on
the verified processed expression scale). Preserve the original array matrix
and retain a mapping/QC record of which GSMs contributed to each patient
summary. Single-array patients contribute their one observed array. Then use
one row/column per patient in comparisons and report the resulting patient
counts by subtype.

If provenance cannot establish that the replicates are equivalent technical
measurements, do not silently average them: resolve the source-processing
history first, or use a statistical model that explicitly accounts for arrays
nested within patients. No aggregation or downstream comparison was performed
in this audit.

## Limitations

- The compatibility audit reports that the exact CDF is unresolved: the GEO
  description says standard AffyCDF **or** Brainarray custom CDF.
- The replicate correlation metrics are descriptive and do not certify
  technical quality or verify prior averaging.
- The pooled median-IQR screen is not group-adjusted and is not a failure test.
- This audit is about sample identity and replicate resolution only. It does
  not biologically validate the dataset or perform differential-expression
  analysis.
