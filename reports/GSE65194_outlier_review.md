# GSE65194 outlier review

## Scope and interpretation

This is a read-only follow-up of the preliminary QC outputs. It uses the
existing sample summary, replicate metrics, manifest, and QC report; expression
values and source files were not changed. No samples were removed, replicate
arrays averaged, or differential-expression analyses performed.

The original QC script screened each array's **median expression across all
54,673 probes** against a single pooled distribution of the 178 array medians.
It calculated the 25th percentile (Q1), 75th percentile (Q3), and
IQR = Q3 − Q1, then flagged values strictly below Q1 − 1.5 × IQR or strictly
above Q3 + 1.5 × IQR.

| Threshold component | Value |
|---|---:|
| Q1 of array medians | 3.1959969000 |
| Q3 of array medians | 3.2403202250 |
| IQR | 0.0443233250 |
| Lower fence (Q1 − 1.5 × IQR) | 3.1295119125 |
| Upper fence (Q3 + 1.5 × IQR) | 3.3068052125 |

All 15 flagged arrays are above the upper fence. None are below the lower
fence. These are **statistical outliers under this pooled screening rule**, not
confirmed technical failures.

Subtype labels below retain the manifest's `sample_group` values. TNBC is
described as TNBC/Basal-like in the dataset annotations; HER2 is shown as
HER2 while the manifest label is `Her2`. Healthy and cell-line samples do not
have tumour subtype or tumour patient IDs.

## Flagged arrays

| GSM accession | Patient ID | Sample title | Sample group | Subtype | Replicate label | Median expression | IQR thresholds used (lower, upper) | Technical replicate partner | Pair Pearson correlation |
|---|---|---|---|---|---|---:|---|---|---:|
| GSM1588986 | TUM013 | TNBC, TUM013_repA | TNBC | TNBC / Basal-like | repA | 3.3118500 | 3.1295119125, 3.3068052125 | GSM1588987 (repB) | 0.99396011 |
| GSM1589121 | TUM123 | Luminal B, TUM123 | Luminal B | Luminal B | — | 3.3235764 | 3.1295119125, 3.3068052125 | None; single array | Not available |
| GSM1589129 | — | CellLine, 184B5 | CellLine | Not applicable (cell line) | — | 3.3082159 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589130 | — | Healthy, H3 | Healthy | Not applicable (healthy) | — | 3.6045629 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589132 | — | Healthy, H13 | Healthy | Not applicable (healthy) | — | 3.7058247 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589136 | — | Healthy, H10 | Healthy | Not applicable (healthy) | — | 3.3431339 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589139 | — | Healthy, H18 | Healthy | Not applicable (healthy) | — | 3.6158316 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589142 | — | Healthy, H12 | Healthy | Not applicable (healthy) | — | 3.4635968 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589144 | — | Healthy, H11 | Healthy | Not applicable (healthy) | — | 3.4031765 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589145 | — | Healthy, H15 | Healthy | Not applicable (healthy) | — | 3.4131943 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589147 | — | CellLine, BT-549 | CellLine | Not applicable (cell line) | — | 3.3301442 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589148 | — | Healthy, H4 | Healthy | Not applicable (healthy) | — | 3.5424792 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589150 | — | Healthy, H14 | Healthy | Not applicable (healthy) | — | 3.3542882 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589151 | — | Healthy, H5 | Healthy | Not applicable (healthy) | — | 3.3533328 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |
| GSM1589152 | — | CellLine, MDA-MB-157 | CellLine | Not applicable (cell line) | — | 3.3131541 | 3.1295119125, 3.3068052125 | None; not a tumour patient | Not available |

## Distribution across sample groups

| Sample group | Flagged | Total arrays | Flagged within group |
|---|---:|---:|---:|
| Healthy | 10 | 11 | 90.9% |
| CellLine | 3 | 14 | 21.4% |
| TNBC | 1 | 55 | 1.8% |
| Luminal B | 1 | 30 | 3.3% |
| Her2 | 0 | 39 | 0.0% |
| Luminal A | 0 | 29 | 0.0% |
| **All groups** | **15** | **178** | **8.4%** |

The flags are strongly concentrated in the non-tumour groups: 10 of 11 healthy
arrays and 3 of 14 cell-line arrays are flagged, compared with 2 of 153 tumour
arrays. This is an unusual distribution, but the rule pools array medians from
different sample types and biological groups into one reference distribution.
The concentration may therefore reflect differences between sample types
rather than poor array quality. This screen is not group-adjusted and does not
identify the cause of the higher medians.

## Closer-inspection priorities

- **GSM1588986 (TUM013, TNBC/Basal-like, repA):** inspect alongside its partner
  GSM1588987. The pair correlation is 0.99396011, which is high, so the median
  flag alone does not suggest a clear disagreement between the two arrays.
  Reviewing both arrays' technical metadata and their per-array distribution
  is reasonable; high correlation does not certify either array.
- **GSM1589121 (TUM123, Luminal B):** this is the other flagged tumour array
  and has no replicate partner in the manifest, so replicate concordance cannot
  be checked. It merits contextual inspection, but the median screen alone is
  not evidence of failure.
- **Healthy arrays:** the 10 flagged healthy samples merit contextual review
  together, especially against the one healthy sample not flagged. Since
  healthy samples make up almost all flags and the threshold was pooled across
  all groups, first consider sample-type distribution and relevant batch or
  processing metadata rather than treating these as 10 independent failures.
- **Cell-line arrays:** the three flagged cell lines may be reviewed in the
  context of the other cell-line arrays and their metadata. Their flags are
  likewise not evidence of technical failure by themselves.

No sample is designated as failed. Confirming a technical failure would require
additional evidence (for example, raw-array diagnostics, processing and batch
records, or independent technical QC), which this median-only investigation
does not provide.
