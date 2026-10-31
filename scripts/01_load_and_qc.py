"""Load the GSE65194 expression matrix and create preliminary QC outputs.

This script describes the existing processed expression values. It does not
normalize data, average or remove samples, or perform differential-expression
analysis.
"""

from __future__ import annotations

import gzip
import hashlib
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


PROJECT_DIR = Path(__file__).resolve().parents[1]
MATRIX_PATH = PROJECT_DIR / "data" / "raw" / "GSE65194_series_matrix.txt.gz"
MANIFEST_PATH = PROJECT_DIR / "data" / "metadata" / "GSE65194_sample_manifest.tsv"
REPORT_DIR = PROJECT_DIR / "reports"

TABLE_BEGIN = "!series_matrix_table_begin"
TABLE_END = "!series_matrix_table_end"
EXPECTED_PROBES = 54_673
EXPECTED_SAMPLES = 178
EXPECTED_GROUP_COUNTS = {
    "TNBC": 55,
    "Her2": 39,
    "Luminal A": 29,
    "Luminal B": 30,
    "Healthy": 11,
    "CellLine": 14,
}
EXPECTED_PATIENTS = 130
EXPECTED_REPLICATE_PAIRS = 23

GROUP_COLORS = {
    "TNBC": "#d95f02",
    "Her2": "#7570b3",
    "Luminal A": "#1b9e77",
    "Luminal B": "#66a61e",
    "Healthy": "#1f78b4",
    "CellLine": "#e7298a",
}


class QCValidationError(ValueError):
    """Raised when an input does not meet the expected QC structure."""


def sha256_file(path: Path) -> str:
    """Return a file's SHA-256 digest without changing it."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_expression_table_bounds(path: Path) -> tuple[int, int]:
    """Find the marker line numbers and validate that the table is present."""
    begin_lines: list[int] = []
    end_lines: list[int] = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            marker = line.strip()
            if marker == TABLE_BEGIN:
                begin_lines.append(line_number)
            elif marker == TABLE_END:
                end_lines.append(line_number)

    if len(begin_lines) != 1 or len(end_lines) != 1:
        raise QCValidationError(
            f"Expected exactly one {TABLE_BEGIN!r} and one {TABLE_END!r} marker; "
            f"found {len(begin_lines)} begin marker(s) and {len(end_lines)} end marker(s)."
        )
    if end_lines[0] <= begin_lines[0] + 1:
        raise QCValidationError("The expression table markers are missing a table header or data.")

    return begin_lines[0], end_lines[0]


def load_expression_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    """Read only the tabular expression section of the compressed GEO file."""
    if not path.is_file():
        raise FileNotFoundError(f"Expression matrix not found: {path}")

    digest_before = sha256_file(path)
    begin_line, end_line = find_expression_table_bounds(path)
    table_line_count = end_line - begin_line - 1

    # Stream past the GEO metadata and let pandas read only the table rows.
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            if line_number == begin_line and line.strip() == TABLE_BEGIN:
                break
        else:
            raise QCValidationError(f"Could not reach {TABLE_BEGIN!r} in the matrix.")

        matrix = pd.read_csv(
            handle,
            sep="\t",
            nrows=table_line_count - 1,
            index_col=0,
            dtype=str,
        )

    if matrix.index.name != "ID_REF":
        raise QCValidationError(
            f"Expected the probe identifier column to be 'ID_REF'; got {matrix.index.name!r}."
        )

    probe_count, sample_count = matrix.shape
    if probe_count != EXPECTED_PROBES or sample_count != EXPECTED_SAMPLES:
        raise QCValidationError(
            "Expression matrix dimensions differ from the expected structure: "
            f"found {probe_count} probes × {sample_count} samples; expected "
            f"{EXPECTED_PROBES} × {EXPECTED_SAMPLES}. No rows or columns were discarded."
        )

    if matrix.index.isna().any():
        raise QCValidationError("The matrix contains missing probe identifiers.")
    probe_ids = pd.Index(matrix.index.astype(str))
    if (probe_ids.str.strip() == "").any():
        raise QCValidationError("The matrix contains missing or blank probe identifiers.")
    if not probe_ids.is_unique:
        duplicate_probes = probe_ids[probe_ids.duplicated()].unique().tolist()[:10]
        raise QCValidationError(f"Duplicate probe identifiers found, for example: {duplicate_probes}")

    sample_ids = [str(column) for column in matrix.columns]
    invalid_sample_ids = [sample_id for sample_id in sample_ids if not re.fullmatch(r"GSM\d+", sample_id)]
    if invalid_sample_ids:
        raise QCValidationError(
            "Expression sample columns are not all GSM accessions; "
            f"examples: {invalid_sample_ids[:10]}"
        )
    if len(set(sample_ids)) != EXPECTED_SAMPLES:
        raise QCValidationError("Expression matrix sample GSM identifiers are not all unique.")

    original_values = matrix
    numeric_matrix = matrix.apply(pd.to_numeric, errors="coerce")
    non_numeric = original_values.notna() & original_values.ne("") & numeric_matrix.isna()
    non_numeric_count = int(non_numeric.to_numpy().sum())
    missing_count = int(numeric_matrix.isna().to_numpy().sum())
    if non_numeric_count:
        examples = original_values.where(non_numeric).stack().head(5).tolist()
        raise QCValidationError(
            f"Found {non_numeric_count} non-numeric expression values; examples: {examples}"
        )
    if missing_count:
        raise QCValidationError(f"Found {missing_count} missing expression values in the matrix.")

    values = numeric_matrix.to_numpy(dtype=np.float64)
    if not np.isfinite(values).all():
        non_finite_count = int((~np.isfinite(values)).sum())
        raise QCValidationError(f"Found {non_finite_count} non-finite expression values.")

    digest_after = sha256_file(path)
    if digest_before != digest_after:
        raise QCValidationError("The compressed expression matrix changed while QC was running.")

    numeric_matrix.index = probe_ids
    numeric_matrix.columns = sample_ids
    details = {
        "probe_count": probe_count,
        "sample_count": sample_count,
        "unique_probe_count": int(probe_ids.nunique()),
        "unique_sample_count": len(set(sample_ids)),
        "missing_count": missing_count,
        "non_numeric_count": non_numeric_count,
        "minimum": float(values.min()),
        "maximum": float(values.max()),
        "sha256": digest_after,
    }
    return numeric_matrix, details


def load_and_align_manifest(path: Path, matrix_sample_ids: list[str]) -> tuple[pd.DataFrame, dict[str, object]]:
    """Validate manifest sample IDs and align records to matrix column order."""
    if not path.is_file():
        raise FileNotFoundError(f"Sample manifest not found: {path}")

    manifest = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    print(f"Manifest columns found: {', '.join(manifest.columns)}")

    required_columns = {
        "gsm_id",
        "tumor_id",
        "sample_group",
        "sample_title",
        "replicate_label",
    }
    missing_columns = sorted(required_columns - set(manifest.columns))
    if missing_columns:
        raise QCValidationError(
            f"Manifest is missing required columns: {missing_columns}. "
            "Expected gsm_id, tumor_id, sample_group, sample_title, and replicate_label."
        )

    if len(manifest) != EXPECTED_SAMPLES:
        raise QCValidationError(
            f"Manifest has {len(manifest)} records; expected {EXPECTED_SAMPLES}."
        )
    if manifest["gsm_id"].eq("").any():
        raise QCValidationError("Manifest contains blank GSM identifiers.")
    if manifest["gsm_id"].duplicated().any():
        duplicates = manifest.loc[manifest["gsm_id"].duplicated(), "gsm_id"].tolist()[:10]
        raise QCValidationError(f"Manifest contains duplicate GSM identifiers: {duplicates}")

    matrix_set = set(matrix_sample_ids)
    manifest_set = set(manifest["gsm_id"])
    missing_manifest = sorted(matrix_set - manifest_set)
    extra_manifest = sorted(manifest_set - matrix_set)
    if missing_manifest or extra_manifest:
        raise QCValidationError(
            "Matrix/manifest sample mapping is incorrect. "
            f"Missing manifest GSMs: {missing_manifest[:10]}; "
            f"manifest GSMs absent from matrix: {extra_manifest[:10]}."
        )

    aligned = manifest.set_index("gsm_id").loc[matrix_sample_ids].reset_index()
    if aligned["gsm_id"].tolist() != matrix_sample_ids:
        raise QCValidationError("Manifest rows could not be aligned exactly to matrix sample order.")

    return aligned, {
        "row_count": len(manifest),
        "unique_gsm_count": int(manifest["gsm_id"].nunique()),
        "columns": list(manifest.columns),
        "matrix_manifest_match": True,
    }


def validate_sample_composition(manifest: pd.DataFrame) -> dict[str, object]:
    """Check expected sample groups and tumour-patient replicate counts."""
    group_counts = manifest["sample_group"].value_counts().to_dict()
    expected_groups = set(EXPECTED_GROUP_COUNTS)
    actual_groups = set(group_counts)
    if actual_groups != expected_groups:
        raise QCValidationError(
            "Unexpected sample groups: "
            f"found {sorted(actual_groups)}, expected {sorted(expected_groups)}."
        )

    group_discrepancies = {
        group: {"expected": expected_count, "actual": int(group_counts.get(group, 0))}
        for group, expected_count in EXPECTED_GROUP_COUNTS.items()
        if int(group_counts.get(group, 0)) != expected_count
    }

    tumour_mask = ~manifest["sample_group"].isin(["Healthy", "CellLine"])
    tumour_rows = manifest.loc[tumour_mask].copy()
    patient_ids = tumour_rows["tumor_id"].str.strip()
    if patient_ids.eq("").any():
        missing_patient_gsms = tumour_rows.loc[patient_ids.eq(""), "gsm_id"].tolist()
        raise QCValidationError(
            f"Tumour arrays have blank patient IDs: {missing_patient_gsms[:10]}"
        )

    arrays_per_patient = tumour_rows.groupby("tumor_id", sort=True)["gsm_id"].count()
    patient_count = int(arrays_per_patient.size)
    single_array_patients = int((arrays_per_patient == 1).sum())
    replicate_pairs = int((arrays_per_patient == 2).sum())
    patients_over_two = arrays_per_patient[arrays_per_patient > 2]

    replicate_discrepancies: list[str] = []
    if patient_count != EXPECTED_PATIENTS:
        replicate_discrepancies.append(
            f"{patient_count} unique tumour patients (expected {EXPECTED_PATIENTS})"
        )
    if replicate_pairs != EXPECTED_REPLICATE_PAIRS:
        replicate_discrepancies.append(
            f"{replicate_pairs} replicate pairs (expected {EXPECTED_REPLICATE_PAIRS})"
        )
    if single_array_patients != 107:
        replicate_discrepancies.append(
            f"{single_array_patients} single-array patients (expected 107)"
        )
    if not patients_over_two.empty:
        replicate_discrepancies.append(
            f"patients with more than two arrays: {patients_over_two.to_dict()}"
        )

    return {
        "group_counts": {group: int(group_counts[group]) for group in EXPECTED_GROUP_COUNTS},
        "group_discrepancies": group_discrepancies,
        "tumour_array_count": int(tumour_mask.sum()),
        "patient_count": patient_count,
        "single_array_patients": single_array_patients,
        "replicate_pair_count": replicate_pairs,
        "replicate_discrepancies": replicate_discrepancies,
        "tumour_mask": tumour_mask,
    }


def build_sample_summary(matrix: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """Calculate descriptive expression summaries for each array."""
    summary = pd.DataFrame(
        {
            "gsm_id": matrix.columns,
            "sample_title": manifest["sample_title"].tolist(),
            "tumor_id": manifest["tumor_id"].tolist(),
            "sample_group": manifest["sample_group"].tolist(),
            "mean_expression": matrix.mean(axis=0).to_numpy(),
            "median_expression": matrix.median(axis=0).to_numpy(),
            "standard_deviation": matrix.std(axis=0, ddof=1).to_numpy(),
            "minimum": matrix.min(axis=0).to_numpy(),
            "maximum": matrix.max(axis=0).to_numpy(),
            "percentile_25": matrix.quantile(0.25, axis=0).to_numpy(),
            "percentile_75": matrix.quantile(0.75, axis=0).to_numpy(),
            "interquartile_range": (
                matrix.quantile(0.75, axis=0) - matrix.quantile(0.25, axis=0)
            ).to_numpy(),
            "missing_values": matrix.isna().sum(axis=0).to_numpy(),
        }
    )
    return summary


def calculate_replicate_metrics(
    matrix: pd.DataFrame, manifest: pd.DataFrame, sample_info: dict[str, object]
) -> pd.DataFrame:
    """Calculate pairwise descriptive metrics without averaging arrays."""
    if sample_info["replicate_pair_count"] != EXPECTED_REPLICATE_PAIRS:
        raise QCValidationError(
            "Replicate metrics were not calculated: expected exactly "
            f"{EXPECTED_REPLICATE_PAIRS} tumour replicate pairs, found "
            f"{sample_info['replicate_pair_count']}."
        )
    if sample_info["replicate_discrepancies"]:
        raise QCValidationError(
            "Tumour patient counts do not match expectations; replicate metrics were not "
            "calculated. " + "; ".join(sample_info["replicate_discrepancies"])
        )

    tumour_rows = manifest.loc[sample_info["tumour_mask"]]
    rows: list[dict[str, object]] = []

    for patient_id, patient_rows in tumour_rows.groupby("tumor_id", sort=True):
        if len(patient_rows) == 1:
            continue
        if len(patient_rows) != 2:
            raise QCValidationError(
                f"Patient {patient_id} has {len(patient_rows)} arrays; exactly two are required."
            )

        pair_rows = patient_rows.sort_values(
            "replicate_label", key=lambda labels: labels.map({"repA": 0, "repB": 1}).fillna(2)
        )
        gsm_a, gsm_b = pair_rows["gsm_id"].tolist()
        if gsm_a == gsm_b:
            raise QCValidationError(f"Patient {patient_id} has a duplicate GSM in its pair.")

        values_a = matrix[gsm_a].to_numpy(dtype=np.float64)
        values_b = matrix[gsm_b].to_numpy(dtype=np.float64)
        difference = values_a - values_b
        correlation = float(np.corrcoef(values_a, values_b)[0, 1])

        rows.append(
            {
                "tumor_id": patient_id,
                "gsm_id_array_1": gsm_a,
                "replicate_label_array_1": pair_rows.iloc[0]["replicate_label"],
                "gsm_id_array_2": gsm_b,
                "replicate_label_array_2": pair_rows.iloc[1]["replicate_label"],
                "probe_count": len(difference),
                "pearson_correlation": correlation,
                "median_absolute_difference": float(np.median(np.abs(difference))),
                "mean_signed_difference_array_1_minus_array_2": float(np.mean(difference)),
                "rmse": float(np.sqrt(np.mean(np.square(difference)))),
            }
        )

    metrics = pd.DataFrame(rows)
    if len(metrics) != EXPECTED_REPLICATE_PAIRS:
        raise QCValidationError(
            f"Calculated {len(metrics)} replicate pairs; expected {EXPECTED_REPLICATE_PAIRS}."
        )
    return metrics


def save_expression_plots(
    matrix: pd.DataFrame, manifest: pd.DataFrame, summary: pd.DataFrame
) -> None:
    """Write a distribution boxplot and a per-array median plot."""
    group_labels = manifest["sample_group"].tolist()
    colors = [GROUP_COLORS[group] for group in group_labels]
    sample_labels = manifest["sample_title"].tolist()

    figure, axis = plt.subplots(figsize=(34, 10))
    boxplot = axis.boxplot(
        matrix.to_numpy().T,
        patch_artist=True,
        showfliers=False,
        widths=0.65,
        medianprops={"color": "black", "linewidth": 0.8},
    )
    for box, color in zip(boxplot["boxes"], colors):
        box.set_facecolor(color)
        box.set_alpha(0.72)
    axis.set_title("Expression distributions across all 178 arrays")
    axis.set_xlabel("Array (GEO sample title)")
    axis.set_ylabel("Processed expression value")
    axis.set_xticks(np.arange(1, len(sample_labels) + 1))
    axis.set_xticklabels(sample_labels, rotation=90, fontsize=6)
    axis.grid(axis="y", alpha=0.25)
    legend_handles = [
        Patch(facecolor=GROUP_COLORS[group], label=group, alpha=0.72)
        for group in EXPECTED_GROUP_COUNTS
    ]
    axis.legend(handles=legend_handles, title="Manifest sample group", loc="upper right")
    figure.tight_layout()
    figure.savefig(REPORT_DIR / "qc_expression_distributions.png", dpi=160)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(24, 7))
    x_values = np.arange(len(summary))
    axis.scatter(
        x_values,
        summary["median_expression"],
        c=colors,
        s=24,
        alpha=0.85,
        edgecolors="none",
    )
    axis.set_title("Median expression for each array")
    axis.set_xlabel("Array in matrix order (GEO sample title)")
    axis.set_ylabel("Median processed expression value")
    axis.set_xticks(x_values)
    axis.set_xticklabels(sample_labels, rotation=90, fontsize=6)
    axis.grid(axis="y", alpha=0.25)
    axis.legend(handles=legend_handles, title="Manifest sample group", loc="best")
    figure.tight_layout()
    figure.savefig(REPORT_DIR / "qc_sample_medians.png", dpi=160)
    plt.close(figure)


def save_replicate_plot(metrics: pd.DataFrame) -> None:
    """Plot each tumour technical-replicate pair's Pearson correlation."""
    figure, axis = plt.subplots(figsize=(14, 7))
    x_values = np.arange(len(metrics))
    axis.scatter(x_values, metrics["pearson_correlation"], color="#386cb0", s=42)
    axis.set_title("Pearson correlation across probes for technical replicate pairs")
    axis.set_xlabel("Tumour patient ID")
    axis.set_ylabel("Pearson correlation")
    axis.set_xticks(x_values)
    axis.set_xticklabels(metrics["tumor_id"], rotation=60, ha="right")
    axis.set_ylim(-1.05, 1.05)
    axis.grid(axis="y", alpha=0.25)
    figure.tight_layout()
    figure.savefig(REPORT_DIR / "qc_replicate_correlations.png", dpi=160)
    plt.close(figure)


def format_table(rows: list[tuple[str, object]]) -> str:
    """Format simple report values as a Markdown table."""
    lines = ["| Item | Result |", "|---|---:|"]
    lines.extend(f"| {name} | {value} |" for name, value in rows)
    return "\n".join(lines)


def write_report(
    matrix: pd.DataFrame,
    matrix_info: dict[str, object],
    manifest_info: dict[str, object],
    sample_info: dict[str, object],
    summary: pd.DataFrame,
    replicate_metrics: pd.DataFrame,
    source_sha256: str,
) -> None:
    """Write a plain-language report of the completed structural QC."""
    group_lines = "\n".join(
        f"- {group}: {count}"
        for group, count in sample_info["group_counts"].items()
    )
    if sample_info["group_discrepancies"]:
        group_warning = "\n".join(
            f"- WARNING: {group}: expected {values['expected']}, found {values['actual']}."
            for group, values in sample_info["group_discrepancies"].items()
        )
    else:
        group_warning = "- All sample-group counts match the expected counts."

    medians = summary["median_expression"]
    median_q1 = float(medians.quantile(0.25))
    median_q3 = float(medians.quantile(0.75))
    median_iqr = median_q3 - median_q1
    lower_screen = median_q1 - 1.5 * median_iqr
    upper_screen = median_q3 + 1.5 * median_iqr
    review_candidates = summary.loc[
        (summary["median_expression"] < lower_screen)
        | (summary["median_expression"] > upper_screen),
        ["gsm_id", "sample_title", "median_expression"],
    ]

    if review_candidates.empty:
        outlier_text = (
            "No arrays fell outside the 1.5×IQR screening range for sample medians. "
            "This is only a simple descriptive screen; it does not prove that no unusual "
            "arrays exist."
        )
    else:
        candidate_items = "; ".join(
            f"{row.gsm_id} ({row.sample_title}, median {row.median_expression:.4f})"
            for row in review_candidates.itertuples(index=False)
        )
        outlier_text = (
            "The following arrays are candidates for visual/follow-up inspection because "
            "their sample median is outside the 1.5×IQR screening range across all arrays: "
            f"{candidate_items}. This is a review flag only, not an automatic outlier "
            "classification; no samples were removed."
        )

    distribution_rows = [
        ("Probe rows", matrix_info["probe_count"]),
        ("Unique probe identifiers", matrix_info["unique_probe_count"]),
        ("Array columns", matrix_info["sample_count"]),
        ("Unique GSM identifiers", matrix_info["unique_sample_count"]),
        ("Missing expression values", matrix_info["missing_count"]),
        ("Non-numeric expression values", matrix_info["non_numeric_count"]),
        ("Minimum expression value", f"{matrix_info['minimum']:.6f}"),
        ("Maximum expression value", f"{matrix_info['maximum']:.6f}"),
        ("Median of array medians", f"{float(medians.median()):.6f}"),
        ("Lowest array median", f"{float(medians.min()):.6f}"),
        ("Highest array median", f"{float(medians.max()):.6f}"),
        ("Median-array IQR", f"{median_iqr:.6f}"),
    ]

    replicate_headers = list(replicate_metrics.columns)
    replicate_table = "\n".join(
        [
            "| " + " | ".join(replicate_headers) + " |",
            "| " + " | ".join("---" for _ in replicate_headers) + " |",
            *[
                "| "
                + " | ".join(
                    f"{value:.6f}" if isinstance(value, (float, np.floating)) else str(value)
                    for value in row
                )
                + " |"
                for row in replicate_metrics.itertuples(index=False, name=None)
            ],
        ]
    )
    if sample_info["group_discrepancies"] or sample_info["replicate_discrepancies"]:
        overall_status = "DID NOT PASS all expected-count checks; see warnings below."
    else:
        overall_status = "PASSED the preliminary structural and descriptive QC checks."

    report = f"""# GSE65194 preliminary data-loading and QC report

## Summary

**Result:** {overall_status}

These are preliminary file-structure, sample-mapping, expression-summary, and
technical-replicate checks. Passing them does **not** biologically validate the
dataset, establish clinical validity, or confirm that the GEO preprocessing
description is correct.

No normalization, log transformation, sample removal, replicate averaging,
differential-expression analysis, or subtype expression comparison was performed.
Expression values were analyzed as stored in the GEO series matrix.

## Software

- Python: {sys.version.split()[0]}
- pandas: {pd.__version__}
- NumPy: {np.__version__}
- matplotlib: {matplotlib.__version__}
- SciPy: not used

## Matrix and manifest validation

{format_table([
    ("Matrix dimensions (probes × arrays)", f"{matrix_info['probe_count']} × {matrix_info['sample_count']}"),
    ("Unique probes", matrix_info["unique_probe_count"]),
    ("Unique matrix GSMs", matrix_info["unique_sample_count"]),
    ("Manifest records", manifest_info["row_count"]),
    ("Unique manifest GSMs", manifest_info["unique_gsm_count"]),
    ("Matrix/manifest exact GSM match", "PASS" if manifest_info["matrix_manifest_match"] else "FAIL"),
    ("Probe index", "ID_REF"),
    ("Matrix minimum", f"{matrix_info['minimum']:.6f}"),
    ("Matrix maximum", f"{matrix_info['maximum']:.6f}"),
    ("Source SHA-256", source_sha256),
])}

The manifest columns found were: `{", ".join(manifest_info["columns"])}`.
The `sample_group` values were preserved as stored. `tumor_id` was used as the
patient identifier; healthy and cell-line labels were kept distinct.

## Sample composition

{group_lines}

Expected-count check:

{group_warning}

- Tumour arrays: {sample_info['tumour_array_count']}
- Unique tumour patients: {sample_info['patient_count']} (expected {EXPECTED_PATIENTS})
- Patients with one array: {sample_info['single_array_patients']} (expected 107)
- Patients with two arrays: {sample_info['replicate_pair_count']} (expected {EXPECTED_REPLICATE_PAIRS})
- Patients with more than two arrays: 0

The matrix/manifest sample-order match was validated by GSM accession. The
manifest was aligned to the expression matrix column order without changing
either input file.

## Expression distribution findings

{format_table(distribution_rows)}

Missing-value check: **{matrix_info['missing_count']}** missing values were found
in the expression matrix. Non-numeric value check: **{matrix_info['non_numeric_count']}**
non-numeric values were found.

{outlier_text}

The distribution and sample-median figures use distinct colors and a legend for
the preserved manifest sample groups. They are descriptive plots and do not
perform subtype comparisons.

## Technical replicate findings

Exactly {len(replicate_metrics)} tumour replicate pairs were evaluated across all
{matrix_info['probe_count']} available probes per pair. Correlation is Pearson's
correlation; difference metrics use array 1 minus array 2 for the signed
difference. No replicate values were averaged.

{replicate_table}

High technical-replicate correlation does not by itself prove that the GEO
processing description's statement about replicate averaging is correct. The
matrix still contains separate GSM columns for these arrays.

## Warnings and unresolved questions

- The GEO record describes GC-RMA processing and mentions averaging technical
  replicates, while this matrix retains separate replicate GSM columns. These
  QC calculations do not resolve that discrepancy.
- The exact Affymetrix CDF used for processing remains unresolved in the GEO
  description, which mentions either a standard AffyCDF or a Brainarray custom
  definition.
- Sample-median screening is a heuristic only. Investigate any candidates with
  relevant technical and biological context before making decisions.
- No differential-expression tests, log fold changes, multiple-testing
  corrections, sample exclusions, or biological conclusions were produced.

## Input integrity

The compressed matrix SHA-256 before and after this run was identical:
`{source_sha256}`. The GEO matrix and sample manifest were read-only inputs.
All generated files were written under the project `reports` directory.
"""
    report_path = REPORT_DIR / "GSE65194_QC_report.md"
    report_path.write_text(report, encoding="utf-8")


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Project: {PROJECT_DIR}")
    print(f"Loading expression matrix: {MATRIX_PATH}")
    matrix, matrix_info = load_expression_matrix(MATRIX_PATH)
    matrix_sample_ids = matrix.columns.tolist()
    print(
        "Matrix loaded: "
        f"{matrix_info['probe_count']} probes × {matrix_info['sample_count']} arrays; "
        f"{matrix_info['missing_count']} missing and "
        f"{matrix_info['non_numeric_count']} non-numeric values."
    )

    print(f"Loading sample manifest: {MANIFEST_PATH}")
    manifest, manifest_info = load_and_align_manifest(MANIFEST_PATH, matrix_sample_ids)
    sample_info = validate_sample_composition(manifest)
    if sample_info["group_discrepancies"]:
        print(f"WARNING: sample-group count discrepancies: {sample_info['group_discrepancies']}")
    if sample_info["replicate_discrepancies"]:
        print("WARNING: " + "; ".join(sample_info["replicate_discrepancies"]))

    # The required replicate count is checked before any pair metrics are calculated.
    replicate_metrics = calculate_replicate_metrics(matrix, manifest, sample_info)
    sample_summary = build_sample_summary(matrix, manifest)

    summary_path = REPORT_DIR / "qc_sample_summary.tsv"
    replicate_path = REPORT_DIR / "qc_replicate_metrics.tsv"
    sample_summary.to_csv(summary_path, sep="\t", index=False, float_format="%.8g")
    replicate_metrics.to_csv(replicate_path, sep="\t", index=False, float_format="%.8g")
    save_expression_plots(matrix, manifest, sample_summary)
    save_replicate_plot(replicate_metrics)
    write_report(
        matrix,
        matrix_info,
        manifest_info,
        sample_info,
        sample_summary,
        replicate_metrics,
        str(matrix_info["sha256"]),
    )

    print("Sample groups:")
    for group, count in sample_info["group_counts"].items():
        print(f"  {group}: {count}")
    print(
        "Tumour patients: "
        f"{sample_info['patient_count']} total; "
        f"{sample_info['single_array_patients']} single-array; "
        f"{sample_info['replicate_pair_count']} replicate pairs."
    )
    print(
        "Replicate correlation range: "
        f"{replicate_metrics['pearson_correlation'].min():.6f} to "
        f"{replicate_metrics['pearson_correlation'].max():.6f}."
    )
    print(f"Matrix SHA-256 unchanged: {matrix_info['sha256']}")
    print("Generated files:")
    for filename in [
        "qc_sample_summary.tsv",
        "qc_replicate_metrics.tsv",
        "qc_expression_distributions.png",
        "qc_sample_medians.png",
        "qc_replicate_correlations.png",
        "GSE65194_QC_report.md",
    ]:
        print(f"  {REPORT_DIR / filename}")
    print("No differential-expression analysis was performed.")


if __name__ == "__main__":
    try:
        main()
    except QCValidationError as error:
        print(f"QC validation error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
