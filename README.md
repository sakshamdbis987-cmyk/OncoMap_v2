# OncoMap v2 — Breast Cancer Expression and Molecular Subtype Explorer

OncoMap v2 is an exploratory, probe-level analysis of public breast cancer
expression datasets. It provides a standalone interactive dashboard and a
reproducible set of Python and R analysis scripts. The results are not
clinically validated and are not intended for diagnosis or treatment decisions.

## Research objectives

This project asks:

1. How do selected expression probes differ between tumour and healthy breast
   tissue?
2. How do those probes vary across four breast cancer subtype groups?
3. Are tumour-versus-normal directions consistent in an independent cohort?

The prespecified panel comprises eight genes—**BRCA1, BRCA2, TP53, ERBB2,
PIK3CA, ESR1, PTEN, and EGFR**—represented by 38 GPL570 probes. Results remain
probe-level: multiple probes are not averaged into a synthetic gene-expression
measurement.

## Datasets and preparation

- **[GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568):**
  104 cancer and 17 normal samples on the GPL570 platform. It is the independent
  cohort used for exact-probe direction and matched-panel comparisons.
- **[GSE65194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65194):**
  54,673 probes across 178 source arrays on GPL570. The cohort includes 130
  tumour patients, 11 healthy controls, and 14 cell-line arrays. The tumour
  subtypes represented among patients are TNBC (41), Her2 (30), Luminal A (29),
  and Luminal B (30).

GSE65194 has 153 tumour arrays for 130 patients: 23 technical-replicate pairs
and 107 singleton arrays. The GEO processing annotation says technical
replicates were averaged, but the deposited matrix retains separate columns for
the arrays in all 23 pairs. The patient-level matrix is therefore a
**project-derived dataset**: the two available values in each validated pair
were arithmetically averaged probe-by-probe. This earlier processing
discrepancy and uncertainty about replicate provenance remain limitations.
Cell lines are excluded from patient-level analyses; the 11 healthy controls
are used only for tumour-versus-healthy comparisons.

The supplied matrices are already processed expression measurements. The
analyses do not repeat log transformation, normalization, or batch correction.
Raw GEO expression data and the bulky derived patient-level matrix are not
included in this repository; see [Reproducibility](#reproducibility).

## Methods

- **Tumour versus healthy:** Welch's independent two-sample t-test on the
  patient-level matrix uses Benjamini–Hochberg (BH) adjustment across 38 target
  probes. Two local limma analyses use alternative repA and repB selections for
  paired patients. Each limma fit adjusts across all 54,673 probes. repA and
  repB are sensitivity analyses, not independent cohorts.
- **Subtype comparisons:** Welch and empirical-Bayes limma analyses compare all
  six pairs among TNBC, Her2, Luminal A, and Luminal B using 130 tumour
  patients. Each method applies BH adjustment to the complete 228-test family
  (38 probes × 6 contrasts). Positive effects indicate higher expression in
  the first-named group.
- **Independent cohort:** Results are joined by exact GPL570 probe ID. The
  matched-panel comparison recalculates BH adjustment separately in each
  cohort across the same 38 probes. No cohort matrices are pooled.

## Summary of findings

- Exact-probe tumour/normal directions were concordant for **24/38** probes
  between cohorts.
- In the matched-38 BH comparison, **15** probes were significant in both
  cohorts, **3** only in GSE42568, **12** only in GSE65194, and **8** in neither.
- For GSE65194 tumour versus healthy, Welch identified **27/38** target probes
  at adjusted *p* < 0.05; limma identified **19/38** for each of repA and repB.
  Direction agreed across methods for all 38 targets in both assignments.
- Across subtype comparisons, Welch identified **95/228** and limma
  **100/228** significant tests. **87** were significant under both methods,
  **13** under limma only, and **8** under Welch only. Effect directions agreed
  for all 228 probe-contrast tests.

These are exploratory associations in retrospective datasets. There is only
one independent comparison cohort; preprocessing, testing families, and
sample sizes differ. The targeted 38-probe panel is not genome-wide validation
of these genes. Probe `217492_s_at` is annotated `PTENP1 /// PTEN` and cannot be
attributed uniquely to PTEN. No laboratory, clinical, or diagnostic validation
is provided, and subtype labels do not establish clinical receptor status.
See the [final research report](reports/OncoMap_v2_Final_Research_Report.md)
for detailed results and limitations.

## Interactive dashboard

Open [dashboard/index.html](dashboard/index.html) locally, or host the project
with GitHub Pages. The standalone dashboard embeds its visualization data and
does not require a CDN, external API, or network access to render. Its GEO
links require internet access.

## Reproducibility

The project was tested with Python 3.15.2, NumPy 2.5.4, pandas 3.1.0, SciPy
1.18.1, Matplotlib 3.11.2, R 4.6.2, and limma 3.68.5. The Python scripts use
NumPy, pandas, SciPy, and Matplotlib; the R scripts require limma. No
dependency lockfile is currently provided.

Download the GSE65194 series matrix from GEO and place it at
`data/raw/GSE65194_series_matrix.txt.gz`. The validated
`data/metadata/GSE65194_sample_manifest.tsv` is included. To rebuild that
manifest, also download the GEO clinical update file to
`data/metadata/GSE65194_clinical_data_update.txt.gz` and run the manifest
builder. From the project root, run the analysis scripts in order:

```text
python data/metadata/build_sample_manifest.py  # optional: rebuild manifest
python scripts/01_load_and_qc.py
python scripts/02_build_patient_level_matrix.py
python scripts/03_differential_expression.py
Rscript scripts/04_limma_sensitivity_analysis.R
Rscript scripts/05_limma_subtype_analysis.R
python scripts/06_build_final_report_dashboard.py
```

The scripts expect their documented inputs and stop when required inputs are
missing or outputs already exist. In particular, the final builder refuses to
overwrite the report or dashboard. Preserve existing outputs before rerunning
the pipeline. The existing GSE42568 comparison outputs are retained in
`reports/`; the source full probe table used for that cohort is not distributed
here.

The six analysis scripts are:

1. `scripts/01_load_and_qc.py` — source matrix and sample quality checks.
2. `scripts/02_build_patient_level_matrix.py` — probe-wise patient-level
   replicate aggregation.
3. `scripts/03_differential_expression.py` — targeted Welch analyses and
   assignment/result outputs.
4. `scripts/04_limma_sensitivity_analysis.R` — repA/repB limma sensitivity
   analyses.
5. `scripts/05_limma_subtype_analysis.R` — patient-level limma subtype
   comparisons.
6. `scripts/06_build_final_report_dashboard.py` — validates existing results
   and builds the final report and dashboard.

Compact manifests, reports, and result tables are retained under `data/` and
`reports/`. The full limma tables are retained because they are linked from the
dashboard and document the genome-wide adjusted results.

## Main project outputs

- [Final research report](reports/OncoMap_v2_Final_Research_Report.md)
- [Interactive dashboard](dashboard/index.html)
- Cross-cohort comparisons:
  [replication](reports/GSE42568_vs_GSE65194_replication.tsv) and
  [matched-panel BH](reports/GSE42568_vs_GSE65194_matched_panel_BH.tsv)
- Tumour/healthy and subtype probe-level result tables are available in
  `reports/`; full repA/repB limma tables are in `data/processed/`.

## Scope and non-diagnostic use

OncoMap v2 is an exploratory computational research project. It does not
establish causation, protein abundance, pathway activity, treatment response,
diagnostic performance, or clinical utility. Its results should not be used for
individual medical decisions.
