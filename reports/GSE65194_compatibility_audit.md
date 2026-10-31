# GSE65194 dataset download and compatibility audit

**Audit date:** 2026-10-31<br>
**Scope:** Download, integrity, dataset metadata, sample composition, probe compatibility, and processing compatibility only. No statistical analysis, normalization, replicate removal, or subtype comparisons were performed.

## Corrected clinical-record mapping in the manifest

The original manifest included a misleading `clinical_record_id` column that copied each array's own GSM accession, which implied a one-to-one clinical record for every array. That was incorrect because the GEO clinical file contains one patient-level record per tumour patient, while the matrix contains 178 arrays including technical replicates. The manifest has been corrected to separate patient-level clinical linkage from array-level direct GSM matching:

- `clinical_record_gsm_id`: patient-level clinical GSM for the tumour patient; blank for healthy/cell-line samples.
- `patient_clinical_record_matched`: `yes` for tumour arrays belonging to a clinically documented patient; `no` otherwise.
- `direct_gsm_clinical_match`: `yes` only when that specific array GSM is present in the clinical file; `no` for its replicate partner or for healthy/cell-line samples.

This preserves all 178 arrays and all metadata while clarifying that the clinical file is a 130-patient tumour registry, not a one-record-per-array table. The 23 extra tumour arrays are replicate partners and therefore match the same patient-level clinical record, but do not carry a direct clinical GSM match of their own.

## Sources and downloaded files

- Matrix: [NCBI GEO series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65194/matrix/GSE65194_series_matrix.txt.gz)
- Clinical annotation: [NCBI GEO clinical update](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65194/suppl/GSE65194_clinical_data_update.txt.gz)

Both URLs returned HTTP 200 before download. Only the two requested gzip files were downloaded; no raw CEL files or archive were downloaded. The compressed source files are retained unchanged.

| File | Compressed bytes | Decompressed bytes | Gzip test | SHA-256 |
|---|---:|---:|---|---|
| `data/raw/GSE65194_series_matrix.txt.gz` | 76464351 | 165378449 | Passed | `F92D2109FE0D64ED056D9673E8733B2AF25173230F3632E370201B6A4054CABA` |
| `data/metadata/GSE65194_clinical_data_update.txt.gz` | 6831 | 158250 | Passed | `16CC3D74954B1311476887DAF54FC2A98F2BF7774E2998499EC218C1A25EA58B` |

Both files are non-empty, decompressed fully without error, and downloaded byte counts matched HTTP content lengths.

## Matrix structure and processing

- **Platform:** GPL570 for all 178 samples.
- **Dimensions:** 54,673 unique probe rows × 178 sample columns (179 columns including `ID_REF`).
- No irregular row widths, duplicate probe IDs, missing measurements, or non-numeric measurements were found. All measurements are numeric.
- **Observed value range:** -1.2226913296569 to 19.9991463988896.
- Matrix rows contain probe IDs but no gene-symbol annotations.

The GEO processing description states:

> Data were analyzed using standard AffyCDF or Brainarray HGU133Plus2_Hs_ENTREZG version 13 custom chipset definition file. Data were normalized using GC-RMA. Batch and hybridization effects were corrected using a linear model in which they were treated as fixed effects. Samples with technical replicates were subsequently averaged.

GC-RMA ordinarily produces processed log2-scale expression; the numeric range is compatible with that scale. The GEO text does not separately spell out a `log2` transformation step, so the scale is supported by the named method, not an explicit log2 label. The exact CDF is unresolved because the description says “standard AffyCDF or Brainarray”.

**Replicate caveat:** the downloaded matrix nevertheless retains 178 separate GSM columns, including titles marked `_repA`/`_repB`. This does not look like a matrix with those arrays collapsed to one column per patient, despite the processing sentence. Keep duplicates distinct in this audit, but do not count them as independent patients or remove/average them yet.

## Sample metadata

The four `!Sample_characteristics_ch1` fields are `sample_group`, `array_batch`, `hybridation batch` (GEO spelling), and `cell line` (blank except for cell-line samples). Sample-group counts: TNBC 55 arrays; Her2 39; Luminal A 29; Luminal B 30; Healthy 11; CellLine 14.

All GSM IDs, exact titles, source names, and all sample characteristics are reported below. `[blank]` indicates an empty value.

| GSM | Sample title | Source name | Characteristics |
|---|---|---|---|
| GSM1588970 | TNBC, TUM001_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B3; hybridation batch: H1; [blank] |
| GSM1588971 | TNBC, TUM001_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588972 | Her2, TUM002_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B1; hybridation batch: H2; [blank] |
| GSM1588973 | Her2, TUM002_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588974 | TNBC, TUM003 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588975 | TNBC, TUM004 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588976 | TNBC, TUM005 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1588977 | TNBC, TUM006_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B2; hybridation batch: H2; [blank] |
| GSM1588978 | TNBC, TUM006_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1588979 | TNBC, TUM007 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588980 | TNBC, TUM008 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588981 | Her2, TUM009 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1588982 | TNBC, TUM010_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B1; hybridation batch: H1; [blank] |
| GSM1588983 | TNBC, TUM010_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588984 | TNBC, TUM011 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1588985 | TNBC, TUM012 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1588986 | TNBC, TUM013_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B4; hybridation batch: H3; [blank] |
| GSM1588987 | TNBC, TUM013_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1588988 | TNBC, TUM014 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1588989 | Her2, TUM015_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B4; hybridation batch: H2; [blank] |
| GSM1588990 | Her2, TUM015_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1588991 | Her2, TUM016_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B2; hybridation batch: H2; [blank] |
| GSM1588992 | Her2, TUM016_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588993 | TNBC, TUM017_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B4; hybridation batch: H2; [blank] |
| GSM1588994 | TNBC, TUM017_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588995 | TNBC, TUM018_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B4; hybridation batch: H1; [blank] |
| GSM1588996 | TNBC, TUM018_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588997 | Her2, TUM019 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1588998 | Her2, TUM020_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B1; hybridation batch: H3; [blank] |
| GSM1588999 | Her2, TUM020_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589000 | Her2, TUM021_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B4; hybridation batch: H1; [blank] |
| GSM1589001 | Her2, TUM021_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589002 | Her2, TUM022 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589003 | TNBC, TUM023 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589004 | TNBC, TUM024 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589005 | TNBC, TUM025 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1589006 | Her2, TUM026_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B5; hybridation batch: H3; [blank] |
| GSM1589007 | Her2, TUM026_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589008 | Her2, TUM027 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589009 | TNBC, TUM028_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B1; hybridation batch: H3; [blank] |
| GSM1589010 | TNBC, TUM028_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589011 | TNBC, TUM029 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1589012 | TNBC, TUM030 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589013 | Her2, TUM031 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1589014 | TNBC, TUM032_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B3; hybridation batch: H2; [blank] |
| GSM1589015 | TNBC, TUM032_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589016 | TNBC, TUM033 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589017 | TNBC, TUM034 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589018 | Her2, TUM035_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B5; hybridation batch: H1; [blank] |
| GSM1589019 | Her2, TUM035_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589020 | Her2, TUM036 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1589021 | Her2, TUM037 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589022 | Her2, TUM038_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B3; hybridation batch: H1; [blank] |
| GSM1589023 | Her2, TUM038_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589024 | TNBC, TUM039 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589025 | Her2, TUM040 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589026 | Her2, TUM041 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589027 | Her2, TUM042 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589028 | TNBC, TUM043_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B5; hybridation batch: H2; [blank] |
| GSM1589029 | TNBC, TUM043_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589030 | Her2, TUM044 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589031 | TNBC, TUM045_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B1; hybridation batch: H2; [blank] |
| GSM1589032 | TNBC, TUM045_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589033 | Her2, TUM046 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589034 | Her2, TUM047 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1589035 | TNBC, TUM048_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B2; hybridation batch: H1; [blank] |
| GSM1589036 | TNBC, TUM048_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589037 | TNBC, TUM049 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589038 | TNBC, TUM050 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589039 | TNBC, TUM051_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B2; hybridation batch: H3; [blank] |
| GSM1589040 | TNBC, TUM051_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589041 | TNBC, TUM052 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589042 | TNBC, TUM053 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589043 | TNBC, TUM054 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589044 | Her2, TUM055 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589045 | TNBC, TUM056 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589046 | Her2, TUM057 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589047 | Her2, TUM058_repA | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B2; hybridation batch: H3; [blank] |
| GSM1589048 | Her2, TUM058_repB | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589049 | Her2, TUM059 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589050 | Her2, TUM060 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B8; hybridation batch: H4; [blank] |
| GSM1589051 | Her2, TUM061 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589052 | Her2, TUM062 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589053 | Her2, TUM063 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589054 | TNBC, TUM064_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B5; hybridation batch: H3; [blank] |
| GSM1589055 | TNBC, TUM064_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589056 | TNBC, TUM065 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589057 | Her2, TUM066 | Breast tumor sample from IHC defined Er- Her2+ subgroup | sample_group: Her2; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589058 | TNBC, TUM067_repA | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B5; hybridation batch: H1; [blank] |
| GSM1589059 | TNBC, TUM067_repB | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589060 | TNBC, TUM068 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589061 | TNBC, TUM069 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589062 | TNBC, TUM070 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B6; hybridation batch: H4; [blank] |
| GSM1589063 | TNBC, TUM071 | Breast tumor sample from IHC defined Triple Negative subgroup | sample_group: TNBC; array_batch: B7; hybridation batch: H4; [blank] |
| GSM1589064 | Luminal A, TUM072 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B2; hybridation batch: H2; [blank] |
| GSM1589065 | Luminal A, TUM073 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B4; hybridation batch: H2; [blank] |
| GSM1589066 | Luminal A, TUM074 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B5; hybridation batch: H2; [blank] |
| GSM1589067 | Luminal A, TUM075 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B3; hybridation batch: H1; [blank] |
| GSM1589068 | Luminal A, TUM076 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B1; hybridation batch: H1; [blank] |
| GSM1589069 | Luminal A, TUM077 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B4; hybridation batch: H2; [blank] |
| GSM1589070 | Luminal A, TUM078 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B3; hybridation batch: H3; [blank] |
| GSM1589071 | Luminal A, TUM079 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B5; hybridation batch: H2; [blank] |
| GSM1589072 | Luminal A, TUM080 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B4; hybridation batch: H3; [blank] |
| GSM1589073 | Luminal A, TUM081 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B1; hybridation batch: H1; [blank] |
| GSM1589074 | Luminal A, TUM082 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B3; hybridation batch: H1; [blank] |
| GSM1589075 | Luminal A, TUM083 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B5; hybridation batch: H3; [blank] |
| GSM1589076 | Luminal A, TUM084 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B5; hybridation batch: H3; [blank] |
| GSM1589077 | Luminal A, TUM085 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B2; hybridation batch: H1; [blank] |
| GSM1589078 | Luminal A, TUM086 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B5; hybridation batch: H1; [blank] |
| GSM1589079 | Luminal A, TUM087 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B3; hybridation batch: H2; [blank] |
| GSM1589080 | Luminal A, TUM088 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B2; hybridation batch: H2; [blank] |
| GSM1589081 | Luminal A, TUM089 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B2; hybridation batch: H3; [blank] |
| GSM1589082 | Luminal A, TUM090 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B2; hybridation batch: H3; [blank] |
| GSM1589083 | Luminal A, TUM091 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B1; hybridation batch: H3; [blank] |
| GSM1589084 | Luminal A, TUM092 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B4; hybridation batch: H3; [blank] |
| GSM1589085 | Luminal A, TUM093 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B4; hybridation batch: H1; [blank] |
| GSM1589086 | Luminal A, TUM094 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B2; hybridation batch: H1; [blank] |
| GSM1589087 | Luminal A, TUM095 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B4; hybridation batch: H2; [blank] |
| GSM1589088 | Luminal A, TUM096 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B4; hybridation batch: H3; [blank] |
| GSM1589089 | Luminal A, TUM097 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B1; hybridation batch: H2; [blank] |
| GSM1589090 | Luminal A, TUM098 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B4; hybridation batch: H1; [blank] |
| GSM1589091 | Luminal A, TUM099 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B3; hybridation batch: H1; [blank] |
| GSM1589092 | Luminal A, TUM100 | Breast tumor sample from IHC defined Luminal A subgroup | sample_group: Luminal A; array_batch: B3; hybridation batch: H2; [blank] |
| GSM1589093 | Luminal B, TUM101 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B5; hybridation batch: H2; [blank] |
| GSM1589094 | Luminal B, TUM102 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B1; hybridation batch: H2; [blank] |
| GSM1589095 | Luminal B, TUM103 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B4; hybridation batch: H3; [blank] |
| GSM1589097 | Luminal B, TUM104 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B1; hybridation batch: H3; [blank] |
| GSM1589098 | Luminal B, TUM105 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B2; hybridation batch: H3; [blank] |
| GSM1589099 | Luminal B, TUM106 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B2; hybridation batch: H2; [blank] |
| GSM1589101 | Luminal B, TUM107 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B5; hybridation batch: H1; [blank] |
| GSM1589103 | Luminal B, TUM108 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B3; hybridation batch: H2; [blank] |
| GSM1589105 | Luminal B, TUM109 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B1; hybridation batch: H1; [blank] |
| GSM1589107 | Luminal B, TUM110 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B2; hybridation batch: H2; [blank] |
| GSM1589109 | Luminal B, TUM111 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B2; hybridation batch: H3; [blank] |
| GSM1589110 | Luminal B, TUM112 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B5; hybridation batch: H3; [blank] |
| GSM1589111 | Luminal B, TUM113 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B1; hybridation batch: H1; [blank] |
| GSM1589112 | Luminal B, TUM114 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B3; hybridation batch: H2; [blank] |
| GSM1589113 | Luminal B, TUM115 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B3; hybridation batch: H2; [blank] |
| GSM1589114 | Luminal B, TUM116 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B3; hybridation batch: H3; [blank] |
| GSM1589115 | Luminal B, TUM117 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B2; hybridation batch: H1; [blank] |
| GSM1589116 | Luminal B, TUM118 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B3; hybridation batch: H3; [blank] |
| GSM1589117 | Luminal B, TUM119 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B5; hybridation batch: H3; [blank] |
| GSM1589118 | Luminal B, TUM120 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B3; hybridation batch: H1; [blank] |
| GSM1589119 | Luminal B, TUM121 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B5; hybridation batch: H1; [blank] |
| GSM1589120 | Luminal B, TUM122 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B1; hybridation batch: H3; [blank] |
| GSM1589121 | Luminal B, TUM123 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B2; hybridation batch: H2; [blank] |
| GSM1589122 | Luminal B, TUM124 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B1; hybridation batch: H2; [blank] |
| GSM1589123 | Luminal B, TUM125 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B4; hybridation batch: H1; [blank] |
| GSM1589124 | Luminal B, TUM126 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B4; hybridation batch: H3; [blank] |
| GSM1589125 | Luminal B, TUM127 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B3; hybridation batch: H1; [blank] |
| GSM1589126 | Luminal B, TUM128 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B2; hybridation batch: H1; [blank] |
| GSM1589127 | Luminal B, TUM129 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B5; hybridation batch: H3; [blank] |
| GSM1589128 | Luminal B, TUM130 | Breast tumor sample from IHC defined Luminal B subgroup | sample_group: Luminal B; array_batch: B4; hybridation batch: H1; [blank] |
| GSM1589129 | CellLine, 184B5 | Sample from breast cancer derived cell line 184B5 | sample_group: CellLine; array_batch: B1; hybridation batch: H1; cell line: 184B5 |
| GSM1589130 | Healthy, H3 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B1; hybridation batch: H2; [blank] |
| GSM1589131 | CellLine, MDA-MB-436 | Sample from breast cancer derived cell line MDA-MB-436 | sample_group: CellLine; array_batch: B1; hybridation batch: H2; cell line: MDA-MB-436 |
| GSM1589132 | Healthy, H13 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B1; hybridation batch: H3; [blank] |
| GSM1589133 | CellLine, HCC1143 | Sample from breast cancer derived cell line HCC1143 | sample_group: CellLine; array_batch: B1; hybridation batch: H3; cell line: HCC1143 |
| GSM1589134 | CellLine, HCC1187 | Sample from breast cancer derived cell line HCC1187 | sample_group: CellLine; array_batch: B2; hybridation batch: H1; cell line: HCC1187 |
| GSM1589135 | Healthy, H9 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B2; hybridation batch: H1; [blank] |
| GSM1589136 | Healthy, H10 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B2; hybridation batch: H2; [blank] |
| GSM1589137 | CellLine, BT20 | Sample from breast cancer derived cell line BT20 | sample_group: CellLine; array_batch: B2; hybridation batch: H2; cell line: BT20 |
| GSM1589138 | CellLine, HCC1937 | Sample from breast cancer derived cell line HCC1937 | sample_group: CellLine; array_batch: B2; hybridation batch: H3; cell line: HCC1937 |
| GSM1589139 | Healthy, H18 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B2; hybridation batch: H3; [blank] |
| GSM1589140 | CellLine, MCF-12A | Sample from breast cancer derived cell line MCF-12A | sample_group: CellLine; array_batch: B3; hybridation batch: H1; cell line: MCF-12A |
| GSM1589141 | CellLine, HCC38 | Sample from breast cancer derived cell line HCC38 | sample_group: CellLine; array_batch: B3; hybridation batch: H2; cell line: HCC38 |
| GSM1589142 | Healthy, H12 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B3; hybridation batch: H3; [blank] |
| GSM1589143 | CellLine, Hs 578T | Sample from breast cancer derived cell line Hs 578T | sample_group: CellLine; array_batch: B3; hybridation batch: H3; cell line: Hs 578T |
| GSM1589144 | Healthy, H11 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B4; hybridation batch: H1; [blank] |
| GSM1589145 | Healthy, H15 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B4; hybridation batch: H1; [blank] |
| GSM1589146 | CellLine, MDA-MB-468 | Sample from breast cancer derived cell line MDA-MB-468 | sample_group: CellLine; array_batch: B4; hybridation batch: H2; cell line: MDA-MB-468 |
| GSM1589147 | CellLine, BT-549 | Sample from breast cancer derived cell line BT-549 | sample_group: CellLine; array_batch: B4; hybridation batch: H3; cell line: BT-549 |
| GSM1589148 | Healthy, H4 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B5; hybridation batch: H1; [blank] |
| GSM1589149 | CellLine, HCC70 | Sample from breast cancer derived cell line HCC70 | sample_group: CellLine; array_batch: B5; hybridation batch: H1; cell line: HCC70 |
| GSM1589150 | Healthy, H14 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B5; hybridation batch: H2; [blank] |
| GSM1589151 | Healthy, H5 | Healthy breast sample from mammoplasty surgery | sample_group: Healthy; array_batch: B5; hybridation batch: H2; [blank] |
| GSM1589152 | CellLine, MDA-MB-157 | Sample from breast cancer derived cell line MDA-MB-157 | sample_group: CellLine; array_batch: B5; hybridation batch: H2; cell line: MDA-MB-157 |
| GSM1589153 | CellLine, MDA-MB-231 | Sample from breast cancer derived cell line MDA-MB-231 | sample_group: CellLine; array_batch: B5; hybridation batch: H3; cell line: MDA-MB-231 |

## Clinical annotation

The clinical file has 130 data rows and 1,024 tab-delimited fields: 40 named fields and 984 blank-name fields. `Histological type` is a duplicate header. Exact named columns with zero-based positions:

- 0 `GSE65194_ID`; 1 `Age at diagnosis(yrs)`; 2 `Menopausal status`; 3 `Height(cm)`; 4 `weight(kg)`; 5 `lateral`; 6 `Tumor type`; 7 `Tumor size(mm)`; 8 `Percentage of tumor cells (%)`; 9 `Histological type` (histopathology); 10 `Margins of resection`; 11 `Necrosis insitu`; 12 `Necrosis in the invasive`; 13 `tumor embolism`; 14 `Mitotic count / 10 fields (X40)`; 15 `Number of lymph nodes analyzed`; 16 `Number of metastatic lymph nodes`; 17 `Total number of involved nodes`; 18 `RO_POSITIVE`; 19 `R0 labeled cells(%)`; 20 `RP_POSITIVE`; 21 `RP labeled cells(%)`; 22 `erbB2_POSITIVE`; 23 `erbB2 labeled cells(%)`; 24 `Histological type` (subtype values); 25 `Treatment with neoadjuvant chemotherapy`; 26 `Surgery`; 27 `Chemotherapy`; 28 `CHEMO_DRUG_GROUP`; 29 `Radiotherapy`; 30 `Hormone therapy`; 31 `HORMO_DRUG_GROUP`; 32 `HERCEPTIN`; 33 `Local recurrence`; 34 `Local recurrence time(month)`; 35 `Metastasis`; 36 `Metastasis time(month)`; 37 `contralateral`; 38 `Death event`; 39 `Death event time(months)`.

- Patient ID: no dedicated patient-ID column; matrix titles contain identifiers such as `TUM001`.
- GSM ID: `GSE65194_ID` (column 0), 130 unique IDs. All 130 matched matrix columns by exact GSM accession; none unmatched. No row-order matching was used.
- Subtype: second `Histological type` (column 24), values Basal-like 41, Her2 + 30, luminal A 29, luminal B 30. Direct GSM mapping agrees with matrix groups TNBC 41, Her2 30, Luminal A 29, Luminal B 30. First `Histological type` (column 9) is histopathology.
- Tissue/sample class: `Tumor type` (column 6) includes primary, recurrent, and contralateral categories. Matrix `sample_group`, title, and source distinguish tumour subtypes, healthy tissue, and cell lines.
- Normal/tumour status: the clinical file has no normal records or explicit normal/tumour flag. Matrix labels healthy samples `sample_group: Healthy`; cell lines are separately labelled `CellLine`.
- Technical replicate field: none in clinical data. Matrix titles carry `_repA`/`_repB` and shared TUM identifiers.

Clinical `Tumor type` counts: primary tumor 107; contralateral 10; `controlateral` 1 (GEO spelling); local recurrence 3; primary tumor bilateral 4; primary tumor second cancer 5.

## Sample composition

| Category | Arrays | Unique tumour patients |
|---|---:|---:|
| TNBC | 55 | 41 |
| HER2 / Her2 | 39 | 30 |
| Luminal B | 30 | 30 |
| Luminal A | 29 | 29 |
| **Tumour total** | **153** | **130** |
| Healthy normal breast | 11 | - |
| Cell lines | 14 | - |
| **All arrays** | **178** | - |

The 153 tumour arrays represent 130 unique patients. There are 23 duplicate patient IDs, each with two arrays (23 extra arrays); none has more than two. The official sample composition is supported. Clinical `Basal-like` labels correspond by GSM to matrix TNBC labels.

## Compatibility with original OncoMap GSE42568

- Both series-matrix records identify **GPL570**.
- GSE42568 is processed; its matrix says: “All microarray data were called using the GC robust multichip average method (GC-RMA). R and Bioconductor - GCRMA package.” The original matrix has 54,675 probe rows × 121 samples; numeric range 2.312811642 to 16.04690273.
- GSE65194 also reports GC-RMA, but additionally documents fixed-effect correction for array and hybridization batch. No batch-correction method is mentioned in GSE42568's processing text, so preprocessing is related but not documented as identical.
- GC-RMA supports log2-scale measurements, but log2 is not separately stated for GSE65194. The different observed ranges do not establish identical scales.
- The datasets can reasonably be analyzed independently and logFC directions compared only after defining compatible groups and contrast direction. Do **not** combine the matrices, assume identical scales, or repeat normalization automatically. Account for duplicates and keep cell lines separate from primary tissue in any future analysis.

### Exact target-probe compatibility

The original `GSE42568.top.table.tsv` was filtered by exact gene-symbol tokens, splitting multi-symbol annotations on `///`. All 38 probes are present in GSE65194 (**38 present, 0 missing**), with multiple probes for all eight genes. All 38 symbol strings agree with the original project's bundled GPL570 annotation. GSE65194's matrix itself has no gene-symbol annotations.

| Target gene | Exact probe IDs | Count |
|---|---|---:|
| BRCA1 | 204531_s_at, 211851_x_at | 2 |
| BRCA2 | 208368_s_at, 214727_at | 2 |
| TP53 | 201746_at, 211300_s_at | 2 |
| ERBB2 | 210930_s_at, 216836_s_at, 234354_x_at | 3 |
| PIK3CA | 204369_at, 231854_at, 235980_at | 3 |
| ESR1 | 205225_at, 211233_x_at, 211234_x_at, 211235_s_at, 211627_x_at, 215551_at, 215552_s_at, 217163_at, 217190_x_at | 9 |
| PTEN | 204053_x_at, 204054_at, 211711_s_at, 217492_s_at, 225363_at, 228006_at, 233254_x_at, 242622_x_at | 8 |
| EGFR | 1565483_at, 1565484_x_at, 201983_s_at, 201984_s_at, 210984_x_at, 211550_at, 211551_at, 211607_x_at, 224999_at | 9 |

`217492_s_at` is annotated `PTENP1///PTEN`; it qualifies by exact PTEN token but remains ambiguous and should not be attributed uniquely to PTEN.

## Outstanding uncertainties and recommendation

1. Exact CDF choice (“standard AffyCDF or Brainarray”) and explicit log2 transformation are unresolved.
2. Processing text says technical replicates were averaged, but the distributed matrix retains separate GSM columns for 23 replicate pairs. Clarify before patient-level statistics.
3. Batch and hybridization correction is documented only for GSE65194, not in the GSE42568 matrix description.
4. GSE65194 has no gene-symbol annotation in its matrix; coverage was checked using probe IDs and the original bundled GPL570 annotation.
5. The clinical file has 984 unnamed fields and a duplicate header; use documented field positions.

**Recommendation: suitable for OncoMap v2 as an independent validation dataset**, particularly for checking shared-probe direction consistency, if later work uses actual GSM/TUM identifiers, treats 130 patients—not 153 arrays—as the tumour biological sample count, keeps 14 cell lines separate, and defines the clinical cohort carefully. It is not appropriate to directly pool matrices or assume identical preprocessing. This audit stops before statistical analysis and awaits approval for Python analysis.

## Audit boundaries

No t-tests, fold changes, multiple-testing corrections, subtype comparisons, replicate removal/averaging, normalization, statistical-analysis code, or changes to the original OncoMap project were performed. The original project was read-only throughout.