"""Build a patient-level expression matrix from the deposited GEO array data.

Tumour technical replicates are averaged probe by probe in this explicit
derived-data step. Single-array tumour patients and healthy controls are kept
unchanged. Cell lines are excluded from this patient-level dataset.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = PROJECT_DIR / "scripts"
MATRIX_PATH = PROJECT_DIR / "data" / "raw" / "GSE65194_series_matrix.txt.gz"
MANIFEST_PATH = PROJECT_DIR / "data" / "metadata" / "GSE65194_sample_manifest.tsv"
QC_METRICS_PATH = PROJECT_DIR / "reports" / "qc_replicate_metrics.tsv"
REPLICATE_AUDIT_PATH = PROJECT_DIR / "reports" / "GSE65194_replicate_resolution_audit.md"
OUTPUT_DIR = PROJECT_DIR / "data" / "processed"
OUTPUT_MATRIX_PATH = OUTPUT_DIR / "GSE65194_patient_level_matrix.tsv.gz"
OUTPUT_MANIFEST_PATH = OUTPUT_DIR / "GSE65194_patient_level_manifest.tsv"
OUTPUT_REPORT_PATH = PROJECT_DIR / "reports" / "GSE65194_patient_level_aggregation_report.md"

EXPECTED_PROBES = 54_673
EXPECTED_ARRAYS = 178
EXPECTED_TUMOUR_ARRAYS = 153
EXPECTED_TUMOUR_PATIENTS = 130
EXPECTED_REPLICATE_PAIRS = 23
EXPECTED_SINGLETON_PATIENTS = 107
EXPECTED_HEALTHY = 11
EXPECTED_CELLLINE = 14
EXPECTED_SUBTYPE_PATIENTS = {
    "TNBC": 41,
    "Her2": 30,
    "Luminal A": 29,
    "Luminal B": 30,
}
CORRELATION_TOLERANCE = 1e-6


class AggregationValidationError(ValueError):
    """Raised when source data do not meet the required aggregation checks."""


def load_expression_with_existing_parser() -> tuple[pd.DataFrame, dict[str, object]]:
    """Use Script 01's GEO marker-aware matrix loader without duplicating it."""
    parser_path = SCRIPTS_DIR / "01_load_and_qc.py"
    if not parser_path.is_file():
        raise FileNotFoundError(f"Existing matrix parser not found: {parser_path}")

    spec = importlib.util.spec_from_file_location("oncomap_qc_loader", parser_path)
    if spec is None or spec.loader is None:
        raise AggregationValidationError(f"Could not load matrix parser: {parser_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.load_expression_matrix(MATRIX_PATH)


def read_required_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """Read inputs and validate the expected file and manifest schemas."""
    matrix, matrix_details = load_expression_with_existing_parser()
    if matrix_details["probe_count"] != EXPECTED_PROBES:
        raise AggregationValidationError(
            f"Expected {EXPECTED_PROBES} probes, found {matrix_details['probe_count']}."
        )
    if matrix_details["sample_count"] != EXPECTED_ARRAYS:
        raise AggregationValidationError(
            f"Expected {EXPECTED_ARRAYS} arrays, found {matrix_details['sample_count']}."
        )

    manifest = pd.read_csv(MANIFEST_PATH, sep="\t", dtype=str, keep_default_na=False)
    required_manifest_columns = {
        "gsm_id",
        "sample_title",
        "tumor_id",
        "replicate_label",
        "sample_group",
    }
    missing_columns = required_manifest_columns - set(manifest.columns)
    if missing_columns:
        raise AggregationValidationError(
            f"Manifest is missing required columns: {sorted(missing_columns)}."
        )

    if len(manifest) != EXPECTED_ARRAYS:
        raise AggregationValidationError(
            f"Expected {EXPECTED_ARRAYS} manifest records, found {len(manifest)}."
        )
    if manifest["gsm_id"].eq("").any() or manifest["gsm_id"].duplicated().any():
        raise AggregationValidationError("Manifest GSM IDs are blank or duplicated.")
    matrix_gsms = set(matrix.columns)
    manifest_gsms = set(manifest["gsm_id"])
    if matrix_gsms != manifest_gsms:
        raise AggregationValidationError(
            "Matrix and manifest GSM sets do not match. "
            f"Missing from manifest: {sorted(matrix_gsms - manifest_gsms)[:10]}; "
            f"not found in matrix: {sorted(manifest_gsms - matrix_gsms)[:10]}."
        )

    # Align metadata to the exact expression-matrix column order by GSM.
    aligned_manifest = manifest.set_index("gsm_id").loc[matrix.columns].reset_index()
    if aligned_manifest["gsm_id"].tolist() != matrix.columns.tolist():
        raise AggregationValidationError("Manifest alignment to matrix GSM order failed.")

    if not QC_METRICS_PATH.is_file():
        raise FileNotFoundError(f"Existing replicate QC metrics not found: {QC_METRICS_PATH}")
    qc_metrics = pd.read_csv(QC_METRICS_PATH, sep="\t")
    required_qc_columns = {
        "tumor_id",
        "gsm_id_array_1",
        "gsm_id_array_2",
        "pearson_correlation",
    }
    missing_qc_columns = required_qc_columns - set(qc_metrics.columns)
    if missing_qc_columns:
        raise AggregationValidationError(
            f"Replicate QC metrics are missing columns: {sorted(missing_qc_columns)}."
        )
    if len(qc_metrics) != EXPECTED_REPLICATE_PAIRS:
        raise AggregationValidationError(
            f"Expected {EXPECTED_REPLICATE_PAIRS} existing QC pair metrics, found {len(qc_metrics)}."
        )
    if qc_metrics["tumor_id"].duplicated().any():
        raise AggregationValidationError("Existing replicate QC metrics contain duplicate patient IDs.")

    if not REPLICATE_AUDIT_PATH.is_file():
        raise FileNotFoundError(f"Replicate resolution audit not found: {REPLICATE_AUDIT_PATH}")
    audit_text = REPLICATE_AUDIT_PATH.read_text(encoding="utf-8")
    averaging_statement = "Samples with technical replicates were subsequently averaged."
    if averaging_statement not in audit_text:
        raise AggregationValidationError(
            "The replicate audit does not contain the expected GEO processing statement."
        )

    return matrix, aligned_manifest, qc_metrics, audit_text


def verify_manifest_against_titles(manifest: pd.DataFrame) -> None:
    """Check that tumour IDs and replicate labels are consistent with titles."""
    for row in manifest.itertuples(index=False):
        group = row.sample_group
        is_tumour = group not in {"Healthy", "CellLine"}
        title_tumour_match = re.search(r"TUM\d+", row.sample_title)
        if is_tumour:
            if not row.tumor_id:
                raise AggregationValidationError(
                    f"Tumour GSM {row.gsm_id} has no patient ID in the manifest."
                )
            if not title_tumour_match or title_tumour_match.group(0) != row.tumor_id:
                raise AggregationValidationError(
                    f"Tumour ID/title mismatch for {row.gsm_id}: "
                    f"{row.tumor_id!r} vs {row.sample_title!r}."
                )
            title_rep_a = "_repA" in row.sample_title
            title_rep_b = "_repB" in row.sample_title
            if title_rep_a and title_rep_b:
                raise AggregationValidationError(
                    f"Sample title has both replicate suffixes: {row.sample_title!r}."
                )
            expected_label = "repA" if title_rep_a else "repB" if title_rep_b else ""
            if row.replicate_label != expected_label:
                raise AggregationValidationError(
                    f"Replicate label/title mismatch for {row.gsm_id}: "
                    f"{row.replicate_label!r} vs {row.sample_title!r}."
                )
        elif row.tumor_id:
            raise AggregationValidationError(
                f"Non-tumour sample {row.gsm_id} incorrectly has tumour ID {row.tumor_id!r}."
            )


def validate_replicate_pairs(
    matrix: pd.DataFrame, manifest: pd.DataFrame, qc_metrics: pd.DataFrame
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Validate all pairs and compare recalculated correlations with prior QC."""
    verify_manifest_against_titles(manifest)

    group_counts = manifest["sample_group"].value_counts().to_dict()
    expected_array_counts = {
        "TNBC": 55,
        "Her2": 39,
        "Luminal A": 29,
        "Luminal B": 30,
        "Healthy": EXPECTED_HEALTHY,
        "CellLine": EXPECTED_CELLLINE,
    }
    if group_counts != expected_array_counts:
        raise AggregationValidationError(
            f"Unexpected sample-group array counts: {group_counts}; "
            f"expected {expected_array_counts}."
        )

    tumour_mask = ~manifest["sample_group"].isin(["Healthy", "CellLine"])
    tumour_manifest = manifest.loc[tumour_mask]
    patient_groups = {
        patient_id: rows.copy()
        for patient_id, rows in tumour_manifest.groupby("tumor_id", sort=True)
    }
    if len(tumour_manifest) != EXPECTED_TUMOUR_ARRAYS:
        raise AggregationValidationError(
            f"Expected {EXPECTED_TUMOUR_ARRAYS} tumour arrays, found {len(tumour_manifest)}."
        )
    if len(patient_groups) != EXPECTED_TUMOUR_PATIENTS:
        raise AggregationValidationError(
            f"Expected {EXPECTED_TUMOUR_PATIENTS} tumour patients, found {len(patient_groups)}."
        )

    array_count_distribution = patient_groups and pd.Series(
        [len(group) for group in patient_groups.values()]
    ).value_counts().to_dict()
    singletons = sum(len(group) == 1 for group in patient_groups.values())
    pairs = sum(len(group) == 2 for group in patient_groups.values())
    invalid_sizes = {
        patient_id: len(group)
        for patient_id, group in patient_groups.items()
        if len(group) not in {1, 2}
    }
    if invalid_sizes:
        raise AggregationValidationError(
            f"Tumour patients must have one or two arrays: {invalid_sizes}."
        )
    if singletons != EXPECTED_SINGLETON_PATIENTS or pairs != EXPECTED_REPLICATE_PAIRS:
        raise AggregationValidationError(
            f"Expected {EXPECTED_SINGLETON_PATIENTS} singleton patients and "
            f"{EXPECTED_REPLICATE_PAIRS} pairs; found {singletons} and {pairs}."
        )

    subtype_patient_counts = (
        tumour_manifest.drop_duplicates("tumor_id")["sample_group"].value_counts().to_dict()
    )
    if subtype_patient_counts != EXPECTED_SUBTYPE_PATIENTS:
        raise AggregationValidationError(
            f"Unexpected tumour patient counts by subtype: {subtype_patient_counts}; "
            f"expected {EXPECTED_SUBTYPE_PATIENTS}."
        )

    metrics_by_patient = qc_metrics.set_index("tumor_id")
    if set(metrics_by_patient.index) != {
        patient_id for patient_id, rows in patient_groups.items() if len(rows) == 2
    }:
        raise AggregationValidationError(
            "Existing QC pair patient IDs do not exactly match the manifest replicate pairs."
        )

    pair_results: list[dict[str, object]] = []
    for patient_id, rows in sorted(patient_groups.items()):
        if len(rows) == 1:
            continue
        if rows["gsm_id"].duplicated().any():
            raise AggregationValidationError(f"Patient {patient_id} has duplicated GSM arrays.")
        if rows["sample_group"].nunique() != 1:
            raise AggregationValidationError(
                f"Replicate pair for {patient_id} has different sample groups/subtypes: "
                f"{rows['sample_group'].tolist()}."
            )
        labels = set(rows["replicate_label"])
        if labels != {"repA", "repB"}:
            raise AggregationValidationError(
                f"Replicate pair for {patient_id} must have distinct repA/repB labels; "
                f"found {sorted(labels)}."
            )

        rows = rows.sort_values("replicate_label")
        gsm_a, gsm_b = rows["gsm_id"].tolist()
        if gsm_a not in matrix.columns or gsm_b not in matrix.columns:
            raise AggregationValidationError(
                f"Replicate GSMs for {patient_id} are not both present in the matrix."
            )

        values_a = matrix[gsm_a].to_numpy(dtype=np.float64)
        values_b = matrix[gsm_b].to_numpy(dtype=np.float64)
        if not np.isfinite(values_a).all() or not np.isfinite(values_b).all():
            raise AggregationValidationError(
                f"Replicate arrays for {patient_id} contain missing or non-finite values."
            )
        calculated_correlation = float(np.corrcoef(values_a, values_b)[0, 1])
        existing = metrics_by_patient.loc[patient_id]
        existing_gsms = {existing["gsm_id_array_1"], existing["gsm_id_array_2"]}
        if existing_gsms != {gsm_a, gsm_b}:
            raise AggregationValidationError(
                f"QC metric GSM IDs for {patient_id} do not match manifest pair "
                f"{gsm_a}, {gsm_b}."
            )
        existing_correlation = float(existing["pearson_correlation"])
        correlation_difference = abs(calculated_correlation - existing_correlation)
        if correlation_difference > CORRELATION_TOLERANCE:
            raise AggregationValidationError(
                f"Calculated and existing QC correlations disagree for {patient_id}: "
                f"{calculated_correlation:.10f} vs {existing_correlation:.10f} "
                f"(difference {correlation_difference:.3g})."
            )

        pair_results.append(
            {
                "patient_id": patient_id,
                "gsm_a": gsm_a,
                "gsm_b": gsm_b,
                "sample_group": rows.iloc[0]["sample_group"],
                "calculated_correlation": calculated_correlation,
                "existing_correlation": existing_correlation,
                "correlation_difference": correlation_difference,
            }
        )

    if len(pair_results) != EXPECTED_REPLICATE_PAIRS:
        raise AggregationValidationError(
            f"Validated {len(pair_results)} replicate pairs; expected {EXPECTED_REPLICATE_PAIRS}."
        )

    return pair_results, {
        "group_array_counts": expected_array_counts,
        "tumour_array_count": len(tumour_manifest),
        "tumour_patient_count": len(patient_groups),
        "singleton_patient_count": singletons,
        "replicate_pair_count": pairs,
        "healthy_count": group_counts["Healthy"],
        "cell_line_count": group_counts["CellLine"],
        "subtype_patient_counts": subtype_patient_counts,
        "array_count_distribution": array_count_distribution,
        "patient_groups": patient_groups,
    }


def create_patient_level_data(
    matrix: pd.DataFrame, manifest: pd.DataFrame, sample_info: dict[str, object]
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Aggregate tumour pairs and retain singleton tumours and healthy controls."""
    manifest_by_gsm = manifest.set_index("gsm_id", drop=False)
    patient_groups = sample_info["patient_groups"]
    tumour_rows = manifest.loc[~manifest["sample_group"].isin(["Healthy", "CellLine"])]
    tumour_order = tumour_rows.drop_duplicates("tumor_id")["tumor_id"].tolist()

    output_columns: dict[str, np.ndarray] = {}
    metadata_rows: list[dict[str, object]] = []
    aggregated_pair_count = 0

    for patient_id in tumour_order:
        rows = patient_groups[patient_id]
        subtype_values = rows["sample_group"].unique().tolist()
        if len(subtype_values) != 1:
            raise AggregationValidationError(
                f"Cannot aggregate {patient_id}; it has inconsistent subtypes {subtype_values}."
            )
        subtype = subtype_values[0]
        source_rows = rows.sort_values("replicate_label")
        source_gsms = source_rows["gsm_id"].tolist()
        sample_id = patient_id

        if len(source_gsms) == 2:
            # The pair was fully validated before this probe-wise arithmetic mean.
            patient_values = matrix[source_gsms].mean(axis=1).to_numpy(dtype=np.float64)
            method = "arithmetic_mean_of_processed_log2_values"
            aggregated_pair_count += 1
        elif len(source_gsms) == 1:
            patient_values = matrix[source_gsms[0]].to_numpy(dtype=np.float64, copy=True)
            method = "retained_unchanged_single_array"
        else:
            raise AggregationValidationError(
                f"Unexpected array count for tumour patient {patient_id}: {len(source_gsms)}."
            )

        output_columns[sample_id] = patient_values
        metadata_rows.append(
            {
                "sample_id": sample_id,
                "sample_type": "tumour",
                "patient_id": patient_id,
                "sample_group": subtype,
                "subtype": subtype,
                "source_gsm_ids": ";".join(source_gsms),
                "source_array_count": len(source_gsms),
                "aggregation_method": method,
            }
        )

    # Healthy controls remain separate GSM-named observations; no cell lines enter.
    healthy_rows = manifest.loc[manifest["sample_group"] == "Healthy"]
    for row in healthy_rows.itertuples(index=False):
        sample_id = row.gsm_id
        output_columns[sample_id] = matrix[sample_id].to_numpy(dtype=np.float64, copy=True)
        metadata_rows.append(
            {
                "sample_id": sample_id,
                "sample_type": "healthy_control",
                "patient_id": "",
                "sample_group": row.sample_group,
                "subtype": "",
                "source_gsm_ids": row.gsm_id,
                "source_array_count": 1,
                "aggregation_method": "retained_unchanged_healthy_control",
            }
        )

    patient_matrix = pd.DataFrame(output_columns, index=matrix.index)
    patient_matrix.index.name = "ID_REF"
    output_manifest = pd.DataFrame(metadata_rows)

    excluded_cell_line_count = int((manifest["sample_group"] == "CellLine").sum())
    if aggregated_pair_count != EXPECTED_REPLICATE_PAIRS:
        raise AggregationValidationError(
            f"Aggregated {aggregated_pair_count} pairs; expected {EXPECTED_REPLICATE_PAIRS}."
        )
    if patient_matrix.shape != (EXPECTED_PROBES, 141):
        raise AggregationValidationError(
            f"Expected patient matrix shape {(EXPECTED_PROBES, 141)}, "
            f"found {patient_matrix.shape}."
        )
    if len(output_manifest) != patient_matrix.shape[1]:
        raise AggregationValidationError(
            "Patient-level matrix and manifest observation counts differ."
        )
    if output_manifest["sample_id"].duplicated().any():
        raise AggregationValidationError("Patient-level manifest sample IDs are duplicated.")
    if output_manifest["sample_id"].tolist() != patient_matrix.columns.tolist():
        raise AggregationValidationError(
            "Patient-level matrix columns and manifest sample order do not match."
        )
    if output_manifest["source_array_count"].sum() != EXPECTED_TUMOUR_ARRAYS + EXPECTED_HEALTHY:
        raise AggregationValidationError(
            "The output observations do not account for all tumour arrays and healthy samples."
        )
    if not np.isfinite(patient_matrix.to_numpy(dtype=np.float64)).all():
        raise AggregationValidationError(
            "The patient-level matrix contains missing or non-finite expression values."
        )

    return patient_matrix, output_manifest, excluded_cell_line_count


def make_report(
    matrix: pd.DataFrame,
    patient_matrix: pd.DataFrame,
    output_manifest: pd.DataFrame,
    pair_results: list[dict[str, object]],
    sample_info: dict[str, object],
    excluded_cell_line_count: int,
    matrix_info: dict[str, object],
) -> str:
    """Create the aggregation report, including each output's source GSM IDs."""
    subtype_rows = "\n".join(
        f"| {subtype} | {count} |"
        for subtype, count in EXPECTED_SUBTYPE_PATIENTS.items()
    )
    pair_table = "\n".join(
        "| {patient_id} | {gsm_a} | {gsm_b} | {sample_group} | "
        "{calculated_correlation:.8f} | {existing_correlation:.8f} | "
        "{correlation_difference:.3g} |".format(**pair)
        for pair in pair_results
    )
    source_table = "\n".join(
        "| {sample_id} | {sample_type} | {patient_id} | {sample_group} | {subtype} | "
        "{source_gsm_ids} | {source_array_count} | {aggregation_method} |".format(
            sample_id=row.sample_id,
            sample_type=row.sample_type,
            patient_id=row.patient_id or "—",
            sample_group=row.sample_group,
            subtype=row.subtype or "—",
            source_gsm_ids=row.source_gsm_ids,
            source_array_count=row.source_array_count,
            aggregation_method=row.aggregation_method,
        )
        for row in output_manifest.itertuples(index=False)
    )

    flagged_sources = {
        "GSM1588986": "TUM013 (paired with GSM1588987)",
        "GSM1589121": "TUM123 (single-array patient)",
    }
    retained_flagged = [
        f"- `{gsm}`: {description}; retained and included."
        for gsm, description in flagged_sources.items()
        if gsm in set(matrix.columns)
    ]

    return f"""# GSE65194 patient-level aggregation report

## Summary

This script created an explicit patient-level derived expression matrix from
the deposited array-level GEO values. It did not modify the original GEO
matrix, source sample manifest, or prior QC outputs.

## Input and output dimensions

| Dataset | Probes | Observations |
|---|---:|---:|
| Deposited expression matrix | {matrix.shape[0]} | {matrix.shape[1]} GSM arrays |
| Derived patient-level matrix | {patient_matrix.shape[0]} | {patient_matrix.shape[1]} observations |

- Probe identifier column: `ID_REF`
- Output matrix columns plus probe ID: {patient_matrix.shape[1] + 1}
- Tumour patients: {sample_info['tumour_patient_count']}
- Healthy controls: {sample_info['healthy_count']}
- Cell-line arrays excluded: {excluded_cell_line_count}
- Missing/non-finite values in derived matrix: 0
- Non-numeric values: none; inputs and output were validated as finite numeric values.

## Validation results

- **PASS:** {matrix_info['probe_count']} probes and {matrix_info['sample_count']} unique GSM arrays loaded by the marker-aware parser from Script 01.
- **PASS:** all matrix GSM IDs match exactly one manifest row; metadata was aligned by GSM.
- **PASS:** {sample_info['tumour_array_count']} tumour arrays represent {sample_info['tumour_patient_count']} unique tumour patients.
- **PASS:** {sample_info['replicate_pair_count']} patients have exactly two arrays and {sample_info['singleton_patient_count']} have one array; no patient has more than two.
- **PASS:** {sample_info['healthy_count']} healthy controls and {sample_info['cell_line_count']} cell-line arrays identified.
- **PASS:** all 23 replicate pairs have matching patient IDs and sample groups/subtypes, distinct repA/repB labels, matching matrix titles, and complete numeric values.
- **PASS:** each recalculated pair correlation agrees with the prior QC metric within {CORRELATION_TOLERANCE:g}.
- **PASS:** derived matrix has {patient_matrix.shape[0]} probes × {patient_matrix.shape[1]} observations; every observation has exactly one manifest record in the same order.
- **PASS:** each of the 23 pairs was aggregated once; all singleton tumour arrays and healthy controls were retained unchanged.

### Tumour patients by preserved subtype label

| Manifest sample group / subtype | Unique patients |
|---|---:|
{subtype_rows}
| **Total** | **{sum(EXPECTED_SUBTYPE_PATIENTS.values())}** |

### Replicate-correlation cross-check

| Patient ID | GSM A | GSM B | Sample group | Recalculated correlation | Existing QC correlation | Absolute difference |
|---|---|---|---|---:|---:|---:|
{pair_table}

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

{chr(10).join(retained_flagged)}

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
{source_table}

## Analysis limitations

This is a documented preparation step only. No differential-expression tests,
subtype comparisons, sample exclusions, or biological validation were
performed. High technical-replicate correlations support profile similarity
but do not resolve the historical GEO averaging statement.
"""


def write_new_outputs(
    patient_matrix: pd.DataFrame, output_manifest: pd.DataFrame, report_text: str
) -> None:
    """Write all outputs without replacing existing derived files."""
    final_paths = [OUTPUT_MATRIX_PATH, OUTPUT_MANIFEST_PATH, OUTPUT_REPORT_PATH]
    existing = [path for path in final_paths if path.exists()]
    if existing:
        raise AggregationValidationError(
            "Refusing to overwrite existing derived output(s): "
            + ", ".join(str(path) for path in existing)
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_PATH = OUTPUT_REPORT_PATH.parent
    REPORT_PATH.mkdir(parents=True, exist_ok=True)

    temporary_paths = [
        OUTPUT_MATRIX_PATH.with_name(OUTPUT_MATRIX_PATH.name + ".tmp"),
        OUTPUT_MANIFEST_PATH.with_name(OUTPUT_MANIFEST_PATH.name + ".tmp"),
        OUTPUT_REPORT_PATH.with_name(OUTPUT_REPORT_PATH.name + ".tmp"),
    ]
    existing_temporary = [path for path in temporary_paths if path.exists()]
    if existing_temporary:
        raise AggregationValidationError(
            "Refusing to overwrite temporary output(s): "
            + ", ".join(str(path) for path in existing_temporary)
        )

    try:
        patient_matrix.to_csv(
            temporary_paths[0],
            sep="\t",
            index=True,
            index_label="ID_REF",
            compression="gzip",
            float_format="%.10g",
        )
        output_manifest.to_csv(
            temporary_paths[1],
            sep="\t",
            index=False,
            lineterminator="\n",
        )
        temporary_paths[2].write_text(report_text, encoding="utf-8")

        # All files are complete before any final output path is created.
        for temporary_path, final_path in zip(temporary_paths, final_paths):
            temporary_path.replace(final_path)
    except Exception:
        for temporary_path in temporary_paths:
            if temporary_path.exists():
                temporary_path.unlink()
        raise


def main() -> None:
    print(f"Project: {PROJECT_DIR}")
    for output_path in [OUTPUT_MATRIX_PATH, OUTPUT_MANIFEST_PATH, OUTPUT_REPORT_PATH]:
        if output_path.exists():
            raise AggregationValidationError(
                f"Output already exists; refusing to overwrite: {output_path}"
            )

    matrix, manifest, qc_metrics, audit_text = read_required_inputs()
    print(f"Loaded expression matrix: {matrix.shape[0]} probes × {matrix.shape[1]} arrays.")
    print(f"Manifest columns: {', '.join(manifest.columns)}")
    print(f"Replicate audit loaded ({len(audit_text)} characters).")

    pair_results, sample_info = validate_replicate_pairs(matrix, manifest, qc_metrics)
    print(
        "Validated tumour arrays/patients: "
        f"{sample_info['tumour_array_count']} arrays, "
        f"{sample_info['tumour_patient_count']} patients, "
        f"{sample_info['replicate_pair_count']} pairs, "
        f"{sample_info['singleton_patient_count']} singletons."
    )

    patient_matrix, output_manifest, excluded_cell_line_count = create_patient_level_data(
        matrix, manifest, sample_info
    )
    matrix_info = {
        "probe_count": matrix.shape[0],
        "sample_count": matrix.shape[1],
    }
    report_text = make_report(
        matrix,
        patient_matrix,
        output_manifest,
        pair_results,
        sample_info,
        excluded_cell_line_count,
        matrix_info,
    )
    write_new_outputs(patient_matrix, output_manifest, report_text)

    print(f"Derived matrix dimensions: {patient_matrix.shape[0]} probes × {patient_matrix.shape[1]} observations.")
    print(
        "Subtype patient counts: "
        + ", ".join(
            f"{subtype}={count}"
            for subtype, count in sample_info["subtype_patient_counts"].items()
        )
    )
    print(
        f"Healthy controls: {sample_info['healthy_count']}; "
        f"cell-line arrays excluded: {excluded_cell_line_count}."
    )
    print("All 23 pairs aggregated exactly once; flagged tumour arrays retained.")
    print("Generated files:")
    print(f"  {OUTPUT_MATRIX_PATH}")
    print(f"  {OUTPUT_MANIFEST_PATH}")
    print(f"  {OUTPUT_REPORT_PATH}")
    print("No differential-expression analysis was performed.")


if __name__ == "__main__":
    try:
        main()
    except (AggregationValidationError, FileNotFoundError) as error:
        print(f"Validation error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
