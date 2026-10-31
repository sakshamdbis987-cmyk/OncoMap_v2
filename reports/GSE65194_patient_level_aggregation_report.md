# GSE65194 patient-level aggregation report

## Summary

This script created an explicit patient-level derived expression matrix from
the deposited array-level GEO values. It did not modify the original GEO
matrix, source sample manifest, or prior QC outputs.

## Input and output dimensions

| Dataset | Probes | Observations |
|---|---:|---:|
| Deposited expression matrix | 54673 | 178 GSM arrays |
| Derived patient-level matrix | 54673 | 141 observations |

- Probe identifier column: `ID_REF`
- Output matrix columns plus probe ID: 142
- Tumour patients: 130
- Healthy controls: 11
- Cell-line arrays excluded: 14
- Missing/non-finite values in derived matrix: 0
- Non-numeric values: none; inputs and output were validated as finite numeric values.

## Validation results

- **PASS:** 54673 probes and 178 unique GSM arrays loaded by the marker-aware parser from Script 01.
- **PASS:** all matrix GSM IDs match exactly one manifest row; metadata was aligned by GSM.
- **PASS:** 153 tumour arrays represent 130 unique tumour patients.
- **PASS:** 23 patients have exactly two arrays and 107 have one array; no patient has more than two.
- **PASS:** 11 healthy controls and 14 cell-line arrays identified.
- **PASS:** all 23 replicate pairs have matching patient IDs and sample groups/subtypes, distinct repA/repB labels, matching matrix titles, and complete numeric values.
- **PASS:** each recalculated pair correlation agrees with the prior QC metric within 1e-06.
- **PASS:** derived matrix has 54673 probes × 141 observations; every observation has exactly one manifest record in the same order.
- **PASS:** each of the 23 pairs was aggregated once; all singleton tumour arrays and healthy controls were retained unchanged.

### Tumour patients by preserved subtype label

| Manifest sample group / subtype | Unique patients |
|---|---:|
| TNBC | 41 |
| Her2 | 30 |
| Luminal A | 29 |
| Luminal B | 30 |
| **Total** | **130** |

### Replicate-correlation cross-check

| Patient ID | GSM A | GSM B | Sample group | Recalculated correlation | Existing QC correlation | Absolute difference |
|---|---|---|---|---:|---:|---:|
| TUM001 | GSM1588970 | GSM1588971 | TNBC | 0.99386507 | 0.99386507 | 1.05e-09 |
| TUM002 | GSM1588972 | GSM1588973 | Her2 | 0.99680579 | 0.99680579 | 4.22e-09 |
| TUM006 | GSM1588977 | GSM1588978 | TNBC | 0.99635015 | 0.99635015 | 4.51e-09 |
| TUM010 | GSM1588982 | GSM1588983 | TNBC | 0.99402480 | 0.99402480 | 1.45e-09 |
| TUM013 | GSM1588986 | GSM1588987 | TNBC | 0.99396011 | 0.99396011 | 4.54e-10 |
| TUM015 | GSM1588989 | GSM1588990 | Her2 | 0.99574838 | 0.99574838 | 2.94e-09 |
| TUM016 | GSM1588991 | GSM1588992 | Her2 | 0.99680104 | 0.99680104 | 1.03e-09 |
| TUM017 | GSM1588993 | GSM1588994 | TNBC | 0.99671069 | 0.99671069 | 1.68e-09 |
| TUM018 | GSM1588995 | GSM1588996 | TNBC | 0.99542977 | 0.99542977 | 3e-09 |
| TUM020 | GSM1588998 | GSM1588999 | Her2 | 0.99695276 | 0.99695276 | 4.66e-09 |
| TUM021 | GSM1589000 | GSM1589001 | Her2 | 0.99654635 | 0.99654635 | 3.66e-09 |
| TUM026 | GSM1589006 | GSM1589007 | Her2 | 0.99640705 | 0.99640705 | 2.97e-09 |
| TUM028 | GSM1589009 | GSM1589010 | TNBC | 0.99595858 | 0.99595858 | 3.15e-09 |
| TUM032 | GSM1589014 | GSM1589015 | TNBC | 0.99483256 | 0.99483256 | 7.59e-10 |
| TUM035 | GSM1589018 | GSM1589019 | Her2 | 0.99672584 | 0.99672584 | 1.38e-09 |
| TUM038 | GSM1589022 | GSM1589023 | Her2 | 0.99611082 | 0.99611082 | 2.76e-10 |
| TUM043 | GSM1589028 | GSM1589029 | TNBC | 0.99553655 | 0.99553655 | 5.12e-10 |
| TUM045 | GSM1589031 | GSM1589032 | TNBC | 0.99635094 | 0.99635094 | 2.79e-09 |
| TUM048 | GSM1589035 | GSM1589036 | TNBC | 0.99667885 | 0.99667885 | 2.69e-09 |
| TUM051 | GSM1589039 | GSM1589040 | TNBC | 0.99622609 | 0.99622609 | 3.51e-09 |
| TUM058 | GSM1589047 | GSM1589048 | Her2 | 0.99657850 | 0.99657850 | 2.9e-09 |
| TUM064 | GSM1589054 | GSM1589055 | TNBC | 0.99405899 | 0.99405899 | 4.65e-09 |
| TUM067 | GSM1589058 | GSM1589059 | TNBC | 0.99572814 | 0.99572814 | 1.41e-09 |

## Aggregation method

For each paired tumour patient, the two deposited processed expression values
for every probe were combined using the arithmetic mean:

`patient_value_for_probe = (array_A_value + array_B_value) / 2`

The inputs are treated as the processed log2 expression values described in
the project audit; no additional log transformation or normalization was
applied. Each singleton tumour patient contributes its original array values
unchanged. Each healthy breast control remains an individual GSM-named
observation with values unchanged. The 14 cell-line arrays are not included
in this patient-level tumour/healthy dataset.

The aggregation methods are recorded per observation in the output manifest as
`arithmetic_mean_of_processed_log2_values`,
`retained_unchanged_single_array`, or
`retained_unchanged_healthy_control`.

## Flagged source arrays retained

The two tumour arrays flagged by the existing pooled median-IQR screen were
not removed:

- `GSM1588986`: TUM013 (paired with GSM1588987); retained and included.
- `GSM1589121`: TUM123 (single-array patient); retained and included.

GSM1588986 is included in the TUM013 pair average with GSM1588987. GSM1589121
is retained unchanged as the single-array observation for TUM123.

## GEO processing provenance

The GEO processing description says that samples with technical replicates
were subsequently averaged, but the downloadable series matrix still contains
separate GSM columns for the replicate arrays. This script explicitly creates
a patient-level derived matrix from those deposited array-level values. It
does **not** establish or claim that the original GEO processing already
averaged the deposited columns; that provenance discrepancy remains
unresolved.

## Source GSMs for every output observation

`sample_group` and tumour `subtype` preserve the manifest labels. Healthy
controls have `sample_type=healthy_control`, `sample_group=Healthy`, a blank
tumour patient ID and blank subtype, and use their GSM as `sample_id`.

| Sample ID | Sample type | Patient ID | Sample group | Subtype | Source GSM IDs | Source array count | Aggregation method |
|---|---|---|---|---|---|---:|---|
| TUM001 | tumour | TUM001 | TNBC | TNBC | GSM1588970;GSM1588971 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM002 | tumour | TUM002 | Her2 | Her2 | GSM1588972;GSM1588973 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM003 | tumour | TUM003 | TNBC | TNBC | GSM1588974 | 1 | retained_unchanged_single_array |
| TUM004 | tumour | TUM004 | TNBC | TNBC | GSM1588975 | 1 | retained_unchanged_single_array |
| TUM005 | tumour | TUM005 | TNBC | TNBC | GSM1588976 | 1 | retained_unchanged_single_array |
| TUM006 | tumour | TUM006 | TNBC | TNBC | GSM1588977;GSM1588978 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM007 | tumour | TUM007 | TNBC | TNBC | GSM1588979 | 1 | retained_unchanged_single_array |
| TUM008 | tumour | TUM008 | TNBC | TNBC | GSM1588980 | 1 | retained_unchanged_single_array |
| TUM009 | tumour | TUM009 | Her2 | Her2 | GSM1588981 | 1 | retained_unchanged_single_array |
| TUM010 | tumour | TUM010 | TNBC | TNBC | GSM1588982;GSM1588983 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM011 | tumour | TUM011 | TNBC | TNBC | GSM1588984 | 1 | retained_unchanged_single_array |
| TUM012 | tumour | TUM012 | TNBC | TNBC | GSM1588985 | 1 | retained_unchanged_single_array |
| TUM013 | tumour | TUM013 | TNBC | TNBC | GSM1588986;GSM1588987 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM014 | tumour | TUM014 | TNBC | TNBC | GSM1588988 | 1 | retained_unchanged_single_array |
| TUM015 | tumour | TUM015 | Her2 | Her2 | GSM1588989;GSM1588990 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM016 | tumour | TUM016 | Her2 | Her2 | GSM1588991;GSM1588992 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM017 | tumour | TUM017 | TNBC | TNBC | GSM1588993;GSM1588994 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM018 | tumour | TUM018 | TNBC | TNBC | GSM1588995;GSM1588996 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM019 | tumour | TUM019 | Her2 | Her2 | GSM1588997 | 1 | retained_unchanged_single_array |
| TUM020 | tumour | TUM020 | Her2 | Her2 | GSM1588998;GSM1588999 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM021 | tumour | TUM021 | Her2 | Her2 | GSM1589000;GSM1589001 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM022 | tumour | TUM022 | Her2 | Her2 | GSM1589002 | 1 | retained_unchanged_single_array |
| TUM023 | tumour | TUM023 | TNBC | TNBC | GSM1589003 | 1 | retained_unchanged_single_array |
| TUM024 | tumour | TUM024 | TNBC | TNBC | GSM1589004 | 1 | retained_unchanged_single_array |
| TUM025 | tumour | TUM025 | TNBC | TNBC | GSM1589005 | 1 | retained_unchanged_single_array |
| TUM026 | tumour | TUM026 | Her2 | Her2 | GSM1589006;GSM1589007 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM027 | tumour | TUM027 | Her2 | Her2 | GSM1589008 | 1 | retained_unchanged_single_array |
| TUM028 | tumour | TUM028 | TNBC | TNBC | GSM1589009;GSM1589010 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM029 | tumour | TUM029 | TNBC | TNBC | GSM1589011 | 1 | retained_unchanged_single_array |
| TUM030 | tumour | TUM030 | TNBC | TNBC | GSM1589012 | 1 | retained_unchanged_single_array |
| TUM031 | tumour | TUM031 | Her2 | Her2 | GSM1589013 | 1 | retained_unchanged_single_array |
| TUM032 | tumour | TUM032 | TNBC | TNBC | GSM1589014;GSM1589015 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM033 | tumour | TUM033 | TNBC | TNBC | GSM1589016 | 1 | retained_unchanged_single_array |
| TUM034 | tumour | TUM034 | TNBC | TNBC | GSM1589017 | 1 | retained_unchanged_single_array |
| TUM035 | tumour | TUM035 | Her2 | Her2 | GSM1589018;GSM1589019 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM036 | tumour | TUM036 | Her2 | Her2 | GSM1589020 | 1 | retained_unchanged_single_array |
| TUM037 | tumour | TUM037 | Her2 | Her2 | GSM1589021 | 1 | retained_unchanged_single_array |
| TUM038 | tumour | TUM038 | Her2 | Her2 | GSM1589022;GSM1589023 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM039 | tumour | TUM039 | TNBC | TNBC | GSM1589024 | 1 | retained_unchanged_single_array |
| TUM040 | tumour | TUM040 | Her2 | Her2 | GSM1589025 | 1 | retained_unchanged_single_array |
| TUM041 | tumour | TUM041 | Her2 | Her2 | GSM1589026 | 1 | retained_unchanged_single_array |
| TUM042 | tumour | TUM042 | Her2 | Her2 | GSM1589027 | 1 | retained_unchanged_single_array |
| TUM043 | tumour | TUM043 | TNBC | TNBC | GSM1589028;GSM1589029 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM044 | tumour | TUM044 | Her2 | Her2 | GSM1589030 | 1 | retained_unchanged_single_array |
| TUM045 | tumour | TUM045 | TNBC | TNBC | GSM1589031;GSM1589032 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM046 | tumour | TUM046 | Her2 | Her2 | GSM1589033 | 1 | retained_unchanged_single_array |
| TUM047 | tumour | TUM047 | Her2 | Her2 | GSM1589034 | 1 | retained_unchanged_single_array |
| TUM048 | tumour | TUM048 | TNBC | TNBC | GSM1589035;GSM1589036 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM049 | tumour | TUM049 | TNBC | TNBC | GSM1589037 | 1 | retained_unchanged_single_array |
| TUM050 | tumour | TUM050 | TNBC | TNBC | GSM1589038 | 1 | retained_unchanged_single_array |
| TUM051 | tumour | TUM051 | TNBC | TNBC | GSM1589039;GSM1589040 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM052 | tumour | TUM052 | TNBC | TNBC | GSM1589041 | 1 | retained_unchanged_single_array |
| TUM053 | tumour | TUM053 | TNBC | TNBC | GSM1589042 | 1 | retained_unchanged_single_array |
| TUM054 | tumour | TUM054 | TNBC | TNBC | GSM1589043 | 1 | retained_unchanged_single_array |
| TUM055 | tumour | TUM055 | Her2 | Her2 | GSM1589044 | 1 | retained_unchanged_single_array |
| TUM056 | tumour | TUM056 | TNBC | TNBC | GSM1589045 | 1 | retained_unchanged_single_array |
| TUM057 | tumour | TUM057 | Her2 | Her2 | GSM1589046 | 1 | retained_unchanged_single_array |
| TUM058 | tumour | TUM058 | Her2 | Her2 | GSM1589047;GSM1589048 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM059 | tumour | TUM059 | Her2 | Her2 | GSM1589049 | 1 | retained_unchanged_single_array |
| TUM060 | tumour | TUM060 | Her2 | Her2 | GSM1589050 | 1 | retained_unchanged_single_array |
| TUM061 | tumour | TUM061 | Her2 | Her2 | GSM1589051 | 1 | retained_unchanged_single_array |
| TUM062 | tumour | TUM062 | Her2 | Her2 | GSM1589052 | 1 | retained_unchanged_single_array |
| TUM063 | tumour | TUM063 | Her2 | Her2 | GSM1589053 | 1 | retained_unchanged_single_array |
| TUM064 | tumour | TUM064 | TNBC | TNBC | GSM1589054;GSM1589055 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM065 | tumour | TUM065 | TNBC | TNBC | GSM1589056 | 1 | retained_unchanged_single_array |
| TUM066 | tumour | TUM066 | Her2 | Her2 | GSM1589057 | 1 | retained_unchanged_single_array |
| TUM067 | tumour | TUM067 | TNBC | TNBC | GSM1589058;GSM1589059 | 2 | arithmetic_mean_of_processed_log2_values |
| TUM068 | tumour | TUM068 | TNBC | TNBC | GSM1589060 | 1 | retained_unchanged_single_array |
| TUM069 | tumour | TUM069 | TNBC | TNBC | GSM1589061 | 1 | retained_unchanged_single_array |
| TUM070 | tumour | TUM070 | TNBC | TNBC | GSM1589062 | 1 | retained_unchanged_single_array |
| TUM071 | tumour | TUM071 | TNBC | TNBC | GSM1589063 | 1 | retained_unchanged_single_array |
| TUM072 | tumour | TUM072 | Luminal A | Luminal A | GSM1589064 | 1 | retained_unchanged_single_array |
| TUM073 | tumour | TUM073 | Luminal A | Luminal A | GSM1589065 | 1 | retained_unchanged_single_array |
| TUM074 | tumour | TUM074 | Luminal A | Luminal A | GSM1589066 | 1 | retained_unchanged_single_array |
| TUM075 | tumour | TUM075 | Luminal A | Luminal A | GSM1589067 | 1 | retained_unchanged_single_array |
| TUM076 | tumour | TUM076 | Luminal A | Luminal A | GSM1589068 | 1 | retained_unchanged_single_array |
| TUM077 | tumour | TUM077 | Luminal A | Luminal A | GSM1589069 | 1 | retained_unchanged_single_array |
| TUM078 | tumour | TUM078 | Luminal A | Luminal A | GSM1589070 | 1 | retained_unchanged_single_array |
| TUM079 | tumour | TUM079 | Luminal A | Luminal A | GSM1589071 | 1 | retained_unchanged_single_array |
| TUM080 | tumour | TUM080 | Luminal A | Luminal A | GSM1589072 | 1 | retained_unchanged_single_array |
| TUM081 | tumour | TUM081 | Luminal A | Luminal A | GSM1589073 | 1 | retained_unchanged_single_array |
| TUM082 | tumour | TUM082 | Luminal A | Luminal A | GSM1589074 | 1 | retained_unchanged_single_array |
| TUM083 | tumour | TUM083 | Luminal A | Luminal A | GSM1589075 | 1 | retained_unchanged_single_array |
| TUM084 | tumour | TUM084 | Luminal A | Luminal A | GSM1589076 | 1 | retained_unchanged_single_array |
| TUM085 | tumour | TUM085 | Luminal A | Luminal A | GSM1589077 | 1 | retained_unchanged_single_array |
| TUM086 | tumour | TUM086 | Luminal A | Luminal A | GSM1589078 | 1 | retained_unchanged_single_array |
| TUM087 | tumour | TUM087 | Luminal A | Luminal A | GSM1589079 | 1 | retained_unchanged_single_array |
| TUM088 | tumour | TUM088 | Luminal A | Luminal A | GSM1589080 | 1 | retained_unchanged_single_array |
| TUM089 | tumour | TUM089 | Luminal A | Luminal A | GSM1589081 | 1 | retained_unchanged_single_array |
| TUM090 | tumour | TUM090 | Luminal A | Luminal A | GSM1589082 | 1 | retained_unchanged_single_array |
| TUM091 | tumour | TUM091 | Luminal A | Luminal A | GSM1589083 | 1 | retained_unchanged_single_array |
| TUM092 | tumour | TUM092 | Luminal A | Luminal A | GSM1589084 | 1 | retained_unchanged_single_array |
| TUM093 | tumour | TUM093 | Luminal A | Luminal A | GSM1589085 | 1 | retained_unchanged_single_array |
| TUM094 | tumour | TUM094 | Luminal A | Luminal A | GSM1589086 | 1 | retained_unchanged_single_array |
| TUM095 | tumour | TUM095 | Luminal A | Luminal A | GSM1589087 | 1 | retained_unchanged_single_array |
| TUM096 | tumour | TUM096 | Luminal A | Luminal A | GSM1589088 | 1 | retained_unchanged_single_array |
| TUM097 | tumour | TUM097 | Luminal A | Luminal A | GSM1589089 | 1 | retained_unchanged_single_array |
| TUM098 | tumour | TUM098 | Luminal A | Luminal A | GSM1589090 | 1 | retained_unchanged_single_array |
| TUM099 | tumour | TUM099 | Luminal A | Luminal A | GSM1589091 | 1 | retained_unchanged_single_array |
| TUM100 | tumour | TUM100 | Luminal A | Luminal A | GSM1589092 | 1 | retained_unchanged_single_array |
| TUM101 | tumour | TUM101 | Luminal B | Luminal B | GSM1589093 | 1 | retained_unchanged_single_array |
| TUM102 | tumour | TUM102 | Luminal B | Luminal B | GSM1589094 | 1 | retained_unchanged_single_array |
| TUM103 | tumour | TUM103 | Luminal B | Luminal B | GSM1589095 | 1 | retained_unchanged_single_array |
| TUM104 | tumour | TUM104 | Luminal B | Luminal B | GSM1589097 | 1 | retained_unchanged_single_array |
| TUM105 | tumour | TUM105 | Luminal B | Luminal B | GSM1589098 | 1 | retained_unchanged_single_array |
| TUM106 | tumour | TUM106 | Luminal B | Luminal B | GSM1589099 | 1 | retained_unchanged_single_array |
| TUM107 | tumour | TUM107 | Luminal B | Luminal B | GSM1589101 | 1 | retained_unchanged_single_array |
| TUM108 | tumour | TUM108 | Luminal B | Luminal B | GSM1589103 | 1 | retained_unchanged_single_array |
| TUM109 | tumour | TUM109 | Luminal B | Luminal B | GSM1589105 | 1 | retained_unchanged_single_array |
| TUM110 | tumour | TUM110 | Luminal B | Luminal B | GSM1589107 | 1 | retained_unchanged_single_array |
| TUM111 | tumour | TUM111 | Luminal B | Luminal B | GSM1589109 | 1 | retained_unchanged_single_array |
| TUM112 | tumour | TUM112 | Luminal B | Luminal B | GSM1589110 | 1 | retained_unchanged_single_array |
| TUM113 | tumour | TUM113 | Luminal B | Luminal B | GSM1589111 | 1 | retained_unchanged_single_array |
| TUM114 | tumour | TUM114 | Luminal B | Luminal B | GSM1589112 | 1 | retained_unchanged_single_array |
| TUM115 | tumour | TUM115 | Luminal B | Luminal B | GSM1589113 | 1 | retained_unchanged_single_array |
| TUM116 | tumour | TUM116 | Luminal B | Luminal B | GSM1589114 | 1 | retained_unchanged_single_array |
| TUM117 | tumour | TUM117 | Luminal B | Luminal B | GSM1589115 | 1 | retained_unchanged_single_array |
| TUM118 | tumour | TUM118 | Luminal B | Luminal B | GSM1589116 | 1 | retained_unchanged_single_array |
| TUM119 | tumour | TUM119 | Luminal B | Luminal B | GSM1589117 | 1 | retained_unchanged_single_array |
| TUM120 | tumour | TUM120 | Luminal B | Luminal B | GSM1589118 | 1 | retained_unchanged_single_array |
| TUM121 | tumour | TUM121 | Luminal B | Luminal B | GSM1589119 | 1 | retained_unchanged_single_array |
| TUM122 | tumour | TUM122 | Luminal B | Luminal B | GSM1589120 | 1 | retained_unchanged_single_array |
| TUM123 | tumour | TUM123 | Luminal B | Luminal B | GSM1589121 | 1 | retained_unchanged_single_array |
| TUM124 | tumour | TUM124 | Luminal B | Luminal B | GSM1589122 | 1 | retained_unchanged_single_array |
| TUM125 | tumour | TUM125 | Luminal B | Luminal B | GSM1589123 | 1 | retained_unchanged_single_array |
| TUM126 | tumour | TUM126 | Luminal B | Luminal B | GSM1589124 | 1 | retained_unchanged_single_array |
| TUM127 | tumour | TUM127 | Luminal B | Luminal B | GSM1589125 | 1 | retained_unchanged_single_array |
| TUM128 | tumour | TUM128 | Luminal B | Luminal B | GSM1589126 | 1 | retained_unchanged_single_array |
| TUM129 | tumour | TUM129 | Luminal B | Luminal B | GSM1589127 | 1 | retained_unchanged_single_array |
| TUM130 | tumour | TUM130 | Luminal B | Luminal B | GSM1589128 | 1 | retained_unchanged_single_array |
| GSM1589130 | healthy_control | — | Healthy | — | GSM1589130 | 1 | retained_unchanged_healthy_control |
| GSM1589132 | healthy_control | — | Healthy | — | GSM1589132 | 1 | retained_unchanged_healthy_control |
| GSM1589135 | healthy_control | — | Healthy | — | GSM1589135 | 1 | retained_unchanged_healthy_control |
| GSM1589136 | healthy_control | — | Healthy | — | GSM1589136 | 1 | retained_unchanged_healthy_control |
| GSM1589139 | healthy_control | — | Healthy | — | GSM1589139 | 1 | retained_unchanged_healthy_control |
| GSM1589142 | healthy_control | — | Healthy | — | GSM1589142 | 1 | retained_unchanged_healthy_control |
| GSM1589144 | healthy_control | — | Healthy | — | GSM1589144 | 1 | retained_unchanged_healthy_control |
| GSM1589145 | healthy_control | — | Healthy | — | GSM1589145 | 1 | retained_unchanged_healthy_control |
| GSM1589148 | healthy_control | — | Healthy | — | GSM1589148 | 1 | retained_unchanged_healthy_control |
| GSM1589150 | healthy_control | — | Healthy | — | GSM1589150 | 1 | retained_unchanged_healthy_control |
| GSM1589151 | healthy_control | — | Healthy | — | GSM1589151 | 1 | retained_unchanged_healthy_control |

## Analysis limitations

This is a documented preparation step only. No differential-expression tests,
subtype comparisons, sample exclusions, or biological validation were
performed. High technical-replicate correlations support profile similarity
but do not resolve the historical GEO averaging statement.
