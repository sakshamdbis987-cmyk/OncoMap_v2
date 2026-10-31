"""Run targeted Welch tests on the validated patient-level OncoMap matrix.

This is a 38-probe targeted analysis, not a genome-wide analysis. It does not
modify expression values or any input files.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import ttest_ind


PROJECT_DIR = Path(__file__).resolve().parents[1]
MATRIX_PATH = PROJECT_DIR / "data" / "processed" / "GSE65194_patient_level_matrix.tsv.gz"
MANIFEST_PATH = PROJECT_DIR / "data" / "processed" / "GSE65194_patient_level_manifest.tsv"
COMPATIBILITY_REPORT_PATH = PROJECT_DIR / "reports" / "GSE65194_compatibility_audit.md"
AGGREGATION_REPORT_PATH = PROJECT_DIR / "reports" / "GSE65194_patient_level_aggregation_report.md"
REPORT_DIR = PROJECT_DIR / "reports"

TUMOUR_HEALTHY_PATH = REPORT_DIR / "GSE65194_tumor_vs_healthy_welch.tsv"
SUBTYPE_PATH = REPORT_DIR / "GSE65194_subtype_pairwise_welch.tsv"
GENE_SUMMARY_PATH = REPORT_DIR / "GSE65194_gene_direction_summary.tsv"
REPORT_PATH = REPORT_DIR / "GSE65194_differential_expression_report.md"

EXPECTED_PROBES = 54_673
EXPECTED_OBSERVATIONS = 141
EXPECTED_TUMOUR_PATIENTS = 130
EXPECTED_HEALTHY_CONTROLS = 11
EXPECTED_SUBTYPE_COUNTS = {
    "TNBC": 41,
    "Her2": 30,
    "Luminal A": 29,
    "Luminal B": 30,
}
SUBTYPES = list(EXPECTED_SUBTYPE_COUNTS)
ALPHA = 0.05

# These probe IDs and annotations are transcribed from the validated target
# table in the compatibility audit. Keep the ambiguous PTEN annotation intact.
TARGET_PROBES: dict[str, list[str]] = {
    "BRCA1": ["204531_s_at", "211851_x_at"],
    "BRCA2": ["208368_s_at", "214727_at"],
    "TP53": ["201746_at", "211300_s_at"],
    "ERBB2": ["210930_s_at", "216836_s_at", "234354_x_at"],
    "PIK3CA": ["204369_at", "231854_at", "235980_at"],
    "ESR1": [
        "205225_at",
        "211233_x_at",
        "211234_x_at",
        "211235_s_at",
        "211627_x_at",
        "215551_at",
        "215552_s_at",
        "217163_at",
        "217190_x_at",
    ],
    "PTEN": [
        "204053_x_at",
        "204054_at",
        "211711_s_at",
        "217492_s_at",
        "225363_at",
        "228006_at",
        "233254_x_at",
        "242622_x_at",
    ],
    "EGFR": [
        "1565483_at",
        "1565484_x_at",
        "201983_s_at",
        "201984_s_at",
        "210984_x_at",
        "211550_at",
        "211551_at",
        "211607_x_at",
        "224999_at",
    ],
}
PROBE_TO_GENE = {
    probe: gene for gene, probes in TARGET_PROBES.items() for probe in probes
}
PROBE_ANNOTATION = {
    probe: ("PTENP1 /// PTEN" if probe == "217492_s_at" else gene)
    for probe, gene in PROBE_TO_GENE.items()
}


class AnalysisValidationError(ValueError):
    """Raised when an analysis input fails a required validation."""


def refuse_existing_outputs() -> None:
    """Avoid replacing any requested analysis output."""
    outputs = [TUMOUR_HEALTHY_PATH, SUBTYPE_PATH, GENE_SUMMARY_PATH, REPORT_PATH]
    existing = [path for path in outputs if path.exists()]
    if existing:
        raise AnalysisValidationError(
            "Refusing to overwrite existing output(s): "
            + ", ".join(str(path) for path in existing)
        )


def load_and_validate_inputs() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Load the patient matrix and metadata, validating expected structure."""
    for path in [
        MATRIX_PATH,
        MANIFEST_PATH,
        COMPATIBILITY_REPORT_PATH,
        AGGREGATION_REPORT_PATH,
    ]:
        if not path.is_file():
            raise FileNotFoundError(f"Required input not found: {path}")

    compatibility_text = COMPATIBILITY_REPORT_PATH.read_text(encoding="utf-8")
    aggregation_text = AGGREGATION_REPORT_PATH.read_text(encoding="utf-8")
    if "217492_s_at" not in compatibility_text or "PTENP1///PTEN" not in compatibility_text:
        raise AnalysisValidationError(
            "Compatibility audit does not contain the expected ambiguous PTEN probe annotation."
        )
    expected_probe_mapping = {
        "BRCA1": ["204531_s_at", "211851_x_at"],
        "BRCA2": ["208368_s_at", "214727_at"],
        "TP53": ["201746_at", "211300_s_at"],
        "ERBB2": ["210930_s_at", "216836_s_at", "234354_x_at"],
        "PIK3CA": ["204369_at", "231854_at", "235980_at"],
        "ESR1": [
            "205225_at",
            "211233_x_at",
            "211234_x_at",
            "211235_s_at",
            "211627_x_at",
            "215551_at",
            "215552_s_at",
            "217163_at",
            "217190_x_at",
        ],
        "PTEN": [
            "204053_x_at",
            "204054_at",
            "211711_s_at",
            "217492_s_at",
            "225363_at",
            "228006_at",
            "233254_x_at",
            "242622_x_at",
        ],
        "EGFR": [
            "1565483_at",
            "1565484_x_at",
            "201983_s_at",
            "201984_s_at",
            "210984_x_at",
            "211550_at",
            "211551_at",
            "211607_x_at",
            "224999_at",
        ],
    }
    audit_probe_mapping: dict[str, list[str]] = {}
    for line in compatibility_text.splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) == 3 and cells[0] in expected_probe_mapping:
            audit_probe_mapping[cells[0]] = [
                probe.strip() for probe in cells[1].split(",")
            ]
    if audit_probe_mapping != expected_probe_mapping:
        raise AnalysisValidationError(
            "The 38-probe target mapping does not exactly match the compatibility audit."
        )
    if "explicit patient-level derived expression matrix" not in aggregation_text:
        raise AnalysisValidationError(
            "Aggregation report does not document the patient-level derived matrix."
        )

    # Read as strings first so conversion errors cannot be silently hidden.
    raw_matrix = pd.read_csv(
        MATRIX_PATH,
        sep="\t",
        compression="gzip",
        index_col=0,
        dtype=str,
        keep_default_na=False,
    )
    if raw_matrix.index.name != "ID_REF":
        raise AnalysisValidationError(
            f"Expected matrix probe index 'ID_REF', found {raw_matrix.index.name!r}."
        )
    if raw_matrix.shape != (EXPECTED_PROBES, EXPECTED_OBSERVATIONS):
        raise AnalysisValidationError(
            f"Expected {EXPECTED_PROBES} probes × {EXPECTED_OBSERVATIONS} observations; "
            f"found {raw_matrix.shape[0]} × {raw_matrix.shape[1]}."
        )
    if not raw_matrix.index.is_unique:
        raise AnalysisValidationError("Patient-level matrix contains duplicate probe IDs.")
    if not raw_matrix.columns.is_unique:
        raise AnalysisValidationError("Patient-level matrix contains duplicate sample IDs.")

    numeric_matrix = raw_matrix.apply(pd.to_numeric, errors="coerce")
    missing_count = int(numeric_matrix.isna().to_numpy().sum())
    non_numeric_count = int(
        (raw_matrix.ne("") & numeric_matrix.isna()).to_numpy().sum()
    )
    if missing_count or non_numeric_count:
        raise AnalysisValidationError(
            "Patient-level matrix contains invalid values: "
            f"{missing_count} missing/NaN, {non_numeric_count} non-numeric."
        )
    if not np.isfinite(numeric_matrix.to_numpy(dtype=np.float64)).all():
        raise AnalysisValidationError("Patient-level matrix contains infinite values.")

    manifest = pd.read_csv(MANIFEST_PATH, sep="\t", dtype=str, keep_default_na=False)
    required_columns = {
        "sample_id",
        "sample_type",
        "patient_id",
        "sample_group",
        "subtype",
        "source_gsm_ids",
        "source_array_count",
        "aggregation_method",
    }
    missing_columns = required_columns - set(manifest.columns)
    if missing_columns:
        raise AnalysisValidationError(
            f"Patient manifest is missing columns: {sorted(missing_columns)}."
        )
    if len(manifest) != EXPECTED_OBSERVATIONS:
        raise AnalysisValidationError(
            f"Expected {EXPECTED_OBSERVATIONS} manifest records, found {len(manifest)}."
        )
    if manifest["sample_id"].eq("").any() or manifest["sample_id"].duplicated().any():
        raise AnalysisValidationError("Patient-level manifest has blank or duplicate sample IDs.")
    if manifest["sample_id"].tolist() != numeric_matrix.columns.tolist():
        raise AnalysisValidationError(
            "Patient-level manifest sample IDs/order do not exactly match matrix columns."
        )

    tumour = manifest["sample_type"] == "tumour"
    healthy = manifest["sample_type"] == "healthy_control"
    if int(tumour.sum()) != EXPECTED_TUMOUR_PATIENTS:
        raise AnalysisValidationError(
            f"Expected {EXPECTED_TUMOUR_PATIENTS} tumour observations, found {int(tumour.sum())}."
        )
    if int(healthy.sum()) != EXPECTED_HEALTHY_CONTROLS:
        raise AnalysisValidationError(
            f"Expected {EXPECTED_HEALTHY_CONTROLS} healthy controls, found {int(healthy.sum())}."
        )
    if int((~(tumour | healthy)).sum()) != 0:
        raise AnalysisValidationError("Unexpected sample types found in patient-level manifest.")
    if manifest.loc[tumour, "patient_id"].eq("").any():
        raise AnalysisValidationError("A tumour observation has a missing patient ID.")
    if manifest.loc[tumour, "patient_id"].duplicated().any():
        raise AnalysisValidationError("Tumour patient IDs are not unique in the derived manifest.")

    subtype_counts = manifest.loc[tumour, "subtype"].value_counts().to_dict()
    if subtype_counts != EXPECTED_SUBTYPE_COUNTS:
        raise AnalysisValidationError(
            f"Tumour subtype counts differ: found {subtype_counts}; "
            f"expected {EXPECTED_SUBTYPE_COUNTS}."
        )
    if not manifest.loc[tumour, "sample_group"].equals(manifest.loc[tumour, "subtype"]):
        raise AnalysisValidationError(
            "Tumour sample_group labels and subtype labels do not match."
        )
    if manifest.loc[healthy, "sample_group"].ne("Healthy").any():
        raise AnalysisValidationError("Healthy control rows do not preserve sample_group='Healthy'.")

    expected_probe_ids = list(PROBE_TO_GENE)
    present_probes = set(numeric_matrix.index)
    missing_probes = sorted(set(expected_probe_ids) - present_probes)
    if missing_probes:
        raise AnalysisValidationError(
            f"Missing target probes from the matrix: {missing_probes}."
        )
    if len(expected_probe_ids) != 38 or len(set(expected_probe_ids)) != 38:
        raise AnalysisValidationError("The configured target probe list is not exactly 38 unique IDs.")

    selected = numeric_matrix.loc[expected_probe_ids]
    if not np.isfinite(selected.to_numpy(dtype=np.float64)).all():
        raise AnalysisValidationError("Selected target observations have missing/non-finite values.")

    return numeric_matrix, manifest, expected_probe_ids


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Adjust a complete family of p-values using BH step-up correction."""
    values = np.asarray(p_values, dtype=np.float64)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise AnalysisValidationError("BH adjustment requires a finite one-dimensional p-value array.")
    if ((values < 0) | (values > 1)).any():
        raise AnalysisValidationError("P-values must be between zero and one.")

    count = len(values)
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    ranks = np.arange(1, count + 1, dtype=np.float64)
    adjusted_sorted = sorted_values * count / ranks
    adjusted_sorted = np.minimum.accumulate(adjusted_sorted[::-1])[::-1]
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)
    adjusted = np.empty(count, dtype=np.float64)
    adjusted[order] = adjusted_sorted
    return adjusted


def welch_test(
    first_values: np.ndarray, second_values: np.ndarray
) -> tuple[float, float]:
    """Return Welch's t-statistic and two-sided p-value."""
    result = ttest_ind(
        first_values,
        second_values,
        equal_var=False,
        alternative="two-sided",
        nan_policy="raise",
    )
    statistic = float(result.statistic)
    p_value = float(result.pvalue)
    if not np.isfinite(statistic) or not np.isfinite(p_value):
        raise AnalysisValidationError(
            "Welch's t-test returned a non-finite statistic or p-value; "
            "inspect zero-variance or invalid group values."
        )
    return statistic, p_value


def run_tumour_vs_healthy(
    matrix: pd.DataFrame, manifest: pd.DataFrame, probe_ids: list[str]
) -> pd.DataFrame:
    """Test tumour against healthy controls and correct across 38 probes."""
    tumour_ids = manifest.loc[manifest["sample_type"] == "tumour", "sample_id"].tolist()
    healthy_ids = manifest.loc[
        manifest["sample_type"] == "healthy_control", "sample_id"
    ].tolist()
    rows: list[dict[str, object]] = []
    for probe_id in probe_ids:
        tumour_values = matrix.loc[probe_id, tumour_ids].to_numpy(dtype=np.float64)
        healthy_values = matrix.loc[probe_id, healthy_ids].to_numpy(dtype=np.float64)
        statistic, p_value = welch_test(tumour_values, healthy_values)
        tumour_mean = float(np.mean(tumour_values))
        healthy_mean = float(np.mean(healthy_values))
        rows.append(
            {
                "probe_id": probe_id,
                "gene_symbol": PROBE_TO_GENE[probe_id],
                "gene_annotation": PROBE_ANNOTATION[probe_id],
                "n_tumour": len(tumour_values),
                "n_healthy": len(healthy_values),
                "mean_tumour_log2_expression": tumour_mean,
                "mean_healthy_log2_expression": healthy_mean,
                "mean_log2_difference_tumour_minus_healthy": tumour_mean - healthy_mean,
                "t_statistic": statistic,
                "p_value_two_sided": p_value,
            }
        )
    results = pd.DataFrame(rows)
    results["bh_adjusted_p_value"] = benjamini_hochberg(
        results["p_value_two_sided"].to_numpy()
    )
    results["significant_bh_0_05"] = results["bh_adjusted_p_value"] < ALPHA
    return results


def run_subtype_comparisons(
    matrix: pd.DataFrame, manifest: pd.DataFrame, probe_ids: list[str]
) -> pd.DataFrame:
    """Run all six subtype contrasts and BH-correct all 228 tests together."""
    tumour_manifest = manifest.loc[manifest["sample_type"] == "tumour"]
    contrast_rows: list[dict[str, object]] = []

    for first_group, second_group in itertools.combinations(SUBTYPES, 2):
        first_ids = tumour_manifest.loc[
            tumour_manifest["subtype"] == first_group, "sample_id"
        ].tolist()
        second_ids = tumour_manifest.loc[
            tumour_manifest["subtype"] == second_group, "sample_id"
        ].tolist()
        if (
            len(first_ids) != EXPECTED_SUBTYPE_COUNTS[first_group]
            or len(second_ids) != EXPECTED_SUBTYPE_COUNTS[second_group]
        ):
            raise AnalysisValidationError(
                f"Unexpected group size in contrast {first_group} vs {second_group}."
            )

        for probe_id in probe_ids:
            first_values = matrix.loc[probe_id, first_ids].to_numpy(dtype=np.float64)
            second_values = matrix.loc[probe_id, second_ids].to_numpy(dtype=np.float64)
            statistic, p_value = welch_test(first_values, second_values)
            first_mean = float(np.mean(first_values))
            second_mean = float(np.mean(second_values))
            contrast_rows.append(
                {
                    "first_subtype": first_group,
                    "second_subtype": second_group,
                    "probe_id": probe_id,
                    "gene_symbol": PROBE_TO_GENE[probe_id],
                    "gene_annotation": PROBE_ANNOTATION[probe_id],
                    "n_first": len(first_values),
                    "n_second": len(second_values),
                    "mean_first_log2_expression": first_mean,
                    "mean_second_log2_expression": second_mean,
                    "mean_log2_difference_first_minus_second": first_mean - second_mean,
                    "t_statistic": statistic,
                    "p_value_two_sided": p_value,
                }
            )

    results = pd.DataFrame(contrast_rows)
    if len(results) != 228:
        raise AnalysisValidationError(
            f"Expected 228 subtype-probe tests, calculated {len(results)}."
        )
    results["bh_adjusted_p_value"] = benjamini_hochberg(
        results["p_value_two_sided"].to_numpy()
    )
    results["significant_bh_0_05"] = results["bh_adjusted_p_value"] < ALPHA
    return results


def make_gene_summary(
    tumour_healthy: pd.DataFrame,
    subtype_results: pd.DataFrame,
    matrix: pd.DataFrame,
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize probe direction without inventing gene-level expression."""
    tumour_manifest = manifest.loc[manifest["sample_type"] == "tumour"]
    subtype_means: dict[str, pd.Series] = {}
    for subtype in SUBTYPES:
        sample_ids = tumour_manifest.loc[
            tumour_manifest["subtype"] == subtype, "sample_id"
        ].tolist()
        subtype_means[subtype] = matrix.loc[list(PROBE_TO_GENE), sample_ids].mean(axis=1)

    rows: list[dict[str, object]] = []
    for gene, probes in TARGET_PROBES.items():
        gene_results = tumour_healthy.loc[tumour_healthy["gene_symbol"] == gene]
        differences = gene_results["mean_log2_difference_tumour_minus_healthy"]
        positive_count = int((differences > 0).sum())
        negative_count = int((differences < 0).sum())
        zero_count = int((differences == 0).sum())
        if positive_count and negative_count:
            direction = "mixed"
        elif positive_count:
            direction = "all_positive"
        elif negative_count:
            direction = "all_negative"
        else:
            direction = "all_zero"

        significant_probes = gene_results.loc[
            gene_results["significant_bh_0_05"], "probe_id"
        ].tolist()
        significant_subtype_tests = subtype_results.loc[
            (subtype_results["gene_symbol"] == gene)
            & subtype_results["significant_bh_0_05"]
        ]
        subtype_significant_labels = [
            f"{row.first_subtype}_vs_{row.second_subtype}:{row.probe_id}"
            for row in significant_subtype_tests.itertuples(index=False)
        ]

        if gene in {"ESR1", "ERBB2"}:
            per_probe_means = []
            for probe in probes:
                means = ", ".join(
                    f"{subtype}={float(subtype_means[subtype].loc[probe]):.6f}"
                    for subtype in SUBTYPES
                )
                per_probe_means.append(f"{probe}[{means}]")
            subtype_pattern = "; ".join(per_probe_means)
        else:
            subtype_pattern = ""

        rows.append(
            {
                "gene_symbol": gene,
                "probe_count": len(probes),
                "positive_probe_count_tumour_vs_healthy": positive_count,
                "negative_probe_count_tumour_vs_healthy": negative_count,
                "zero_difference_probe_count_tumour_vs_healthy": zero_count,
                "tumour_vs_healthy_direction_consistency": direction,
                "tumour_vs_healthy_bh_significant_probes": ";".join(significant_probes),
                "subtype_pairwise_bh_significant_probes": ";".join(
                    subtype_significant_labels
                ),
                "ESR1_ERBB2_probe_means_by_subtype": subtype_pattern,
                "probe_level_note": (
                    "217492_s_at annotation is PTENP1 /// PTEN; ambiguous"
                    if gene == "PTEN"
                    else "Probe measurements remain separate; no gene-level averaging."
                ),
            }
        )
    return pd.DataFrame(rows)


def dataframe_markdown(frame: pd.DataFrame, floatfmt: str = ".6g") -> str:
    """Format a DataFrame as a Markdown table without requiring tabulate."""
    headers = [str(column) for column in frame.columns]
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in frame.itertuples(index=False, name=None):
        values = []
        for value in row:
            if isinstance(value, (float, np.floating)):
                values.append(format(float(value), floatfmt))
            elif isinstance(value, (bool, np.bool_)):
                values.append("yes" if value else "no")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def build_report(
    tumour_healthy: pd.DataFrame,
    subtype_results: pd.DataFrame,
    gene_summary: pd.DataFrame,
    manifest: pd.DataFrame,
) -> str:
    """Create the methods/results report from the actual output tables."""
    significant_th = tumour_healthy.loc[tumour_healthy["significant_bh_0_05"]]
    significant_subtype = subtype_results.loc[subtype_results["significant_bh_0_05"]]
    contrast_counts = (
        significant_subtype.groupby(["first_subtype", "second_subtype"], sort=False)
        .size()
        .to_dict()
    )

    contrast_summary_rows = []
    for first, second in itertools.combinations(SUBTYPES, 2):
        contrast_summary_rows.append(
            f"| {first} - {second} | {EXPECTED_SUBTYPE_COUNTS[first]} | "
            f"{EXPECTED_SUBTYPE_COUNTS[second]} | "
            f"{contrast_counts.get((first, second), 0)} |"
        )

    esr1_erbb2_summary = gene_summary.loc[
        gene_summary["gene_symbol"].isin(["ESR1", "ERBB2"]),
        [
            "gene_symbol",
            "probe_count",
            "subtype_pairwise_bh_significant_probes",
            "ESR1_ERBB2_probe_means_by_subtype",
        ],
    ]

    direction_columns = [
        "gene_symbol",
        "probe_count",
        "positive_probe_count_tumour_vs_healthy",
        "negative_probe_count_tumour_vs_healthy",
        "tumour_vs_healthy_direction_consistency",
        "tumour_vs_healthy_bh_significant_probes",
    ]

    return f"""# GSE65194 targeted differential-expression report

## Scope and caution

This is a **targeted analysis of 38 validated GPL570 probe sets** representing
eight genes. It is not genome-wide differential-expression analysis. Results
are probe-level; multiple probes for one gene are not averaged into an invented
gene-expression value.

The analysis is exploratory and does not establish causation, diagnostic
performance, or independent validation. A statistically significant p-value
does not by itself validate a biological conclusion.

No input matrix, manifest, compressed GEO source, or prior report was modified.
No samples were removed. The patient-level matrix already represents the
documented aggregation of replicate pairs; this script did not re-average
replicates.

## Inputs and validation

- Patient-level matrix: 54,673 probes × 141 observations.
- Tumour patients: 130; healthy controls: 11.
- Tumour subtype counts: TNBC 41, Her2 30, Luminal A 29, Luminal B 30.
- All 141 sample identifiers matched exactly between matrix columns and
  manifest rows, with no duplicates.
- All 38 target probes were present.
- No missing, non-numeric, or infinite expression values were found.
- SciPy version: {scipy.__version__}.

## Methods

For each probe, two-sided Welch independent-samples t-tests were run using
`scipy.stats.ttest_ind(equal_var=False)`. The expression values are the
processed log2-scale values as stored in the derived patient-level matrix;
there was no additional transformation or normalization.

### Tumour versus healthy

All 130 tumour patients were compared with all 11 healthy controls. The
reported mean log2 difference is tumour mean minus healthy-control mean.
Positive values indicate higher mean expression in tumour; negative values
indicate lower mean expression in tumour. Benjamini-Hochberg (BH) correction
was applied across the complete family of **38 probe tests**. Significance is
defined as BH-adjusted p < 0.05.

### Pairwise subtype comparisons

All six two-sided contrasts were tested: TNBC−Her2, TNBC−Luminal A,
TNBC−Luminal B, Her2−Luminal A, Her2−Luminal B, and Luminal A−Luminal B.
The reported difference is first subtype mean minus second subtype mean.
Benjamini-Hochberg correction was applied once across the complete family of
**228 tests (6 contrasts × 38 probes)**, not separately within each contrast.
Significance is defined as adjusted p < 0.05.

The recorded subtypes are analyzed as listed; this analysis does not infer
ER-positive status from Luminal A or Luminal B labels.

## Tumour versus healthy results

- Probes passing BH-adjusted p < 0.05: **{len(significant_th)} of 38**.
- Positive tumour-minus-healthy differences: **{int((tumour_healthy['mean_log2_difference_tumour_minus_healthy'] > 0).sum())}**.
- Negative tumour-minus-healthy differences: **{int((tumour_healthy['mean_log2_difference_tumour_minus_healthy'] < 0).sum())}**.

Significant probe results:

{dataframe_markdown(significant_th[['probe_id', 'gene_symbol', 'n_tumour', 'n_healthy', 'mean_log2_difference_tumour_minus_healthy', 'p_value_two_sided', 'bh_adjusted_p_value']]) if len(significant_th) else 'No probes met the BH-adjusted p < 0.05 threshold.'}

Complete probe-level results are in
`GSE65194_tumor_vs_healthy_welch.tsv`.

## Pairwise subtype results

| Contrast (first minus second) | First n | Second n | BH-significant probes in contrast |
|---|---:|---:|---:|
{chr(10).join(contrast_summary_rows)}

There are **{len(significant_subtype)} of 228** subtype-probe tests with
BH-adjusted p < 0.05 across the full family. Complete results, including
non-significant probes, are in
`GSE65194_subtype_pairwise_welch.tsv`.

## Probe directions by gene

The following direction summary is based on tumour-minus-healthy mean
differences. The sign counts refer to separate probes, not a combined gene
measurement.

{dataframe_markdown(gene_summary[direction_columns])}

`217492_s_at` retains the annotation `PTENP1 /// PTEN`; because the annotation
is ambiguous, it should not be attributed uniquely to PTEN. All individual
probes remain separately visible in both result tables.

## ESR1 and ERBB2 across the recorded subtypes

For ESR1 and ERBB2, the gene summary records each probe's mean expression
within TNBC, Her2, Luminal A, and Luminal B separately. These values are
probe-specific; no gene-level average is calculated. Adjusted-significant
pairwise subtype probes are listed below.

{dataframe_markdown(esr1_erbb2_summary)}

Use the full subtype result table for each individual probe's six contrasts,
mean difference, and family-wide adjusted p-value. A subtype mean pattern is
descriptive and does not by itself establish receptor status or biological
mechanism.

## Limitations

- This tests only the 38 selected probes; it is not genome-wide screening.
- Multiple probes may measure different transcript regions or have ambiguous
  annotations. Probe-level direction disagreement is retained rather than
  concealed by averaging.
- The healthy-control sample size is 11, compared with 130 tumour patients;
  group sizes are reported for transparency.
- The upstream audit notes that the exact CDF and explicit log2 transformation
  history are not fully resolved. The analysis uses the processed values as
  supplied and does not renormalize them.
- GEO's processing text mentions averaging technical replicates even though
  separate GSM columns are present in its downloadable matrix. The patient-
  level matrix is the project's explicit documented aggregation; this targeted
  analysis does not resolve GEO's original processing provenance.
- Luminal labels were not converted to ER-positive labels. Clinical receptor
  status is not inferred by this subtype analysis.
- Statistical significance is not proof of biological validation, causality,
  or diagnostic performance.

## Outputs

- `GSE65194_tumor_vs_healthy_welch.tsv`
- `GSE65194_subtype_pairwise_welch.tsv`
- `GSE65194_gene_direction_summary.tsv`
- `GSE65194_differential_expression_report.md`
"""


def save_outputs(
    tumour_healthy: pd.DataFrame,
    subtype_results: pd.DataFrame,
    gene_summary: pd.DataFrame,
    report: str,
) -> None:
    """Write all analysis outputs only after checking no targets exist."""
    refuse_existing_outputs()
    temporary_paths = [
        TUMOUR_HEALTHY_PATH.with_suffix(TUMOUR_HEALTHY_PATH.suffix + ".tmp"),
        SUBTYPE_PATH.with_suffix(SUBTYPE_PATH.suffix + ".tmp"),
        GENE_SUMMARY_PATH.with_suffix(GENE_SUMMARY_PATH.suffix + ".tmp"),
        REPORT_PATH.with_suffix(REPORT_PATH.suffix + ".tmp"),
    ]
    existing_temporary = [path for path in temporary_paths if path.exists()]
    if existing_temporary:
        raise AnalysisValidationError(
            "Refusing to overwrite temporary output(s): "
            + ", ".join(str(path) for path in existing_temporary)
        )

    tumour_healthy.to_csv(temporary_paths[0], sep="\t", index=False, float_format="%.12g")
    subtype_results.to_csv(temporary_paths[1], sep="\t", index=False, float_format="%.12g")
    gene_summary.to_csv(temporary_paths[2], sep="\t", index=False, lineterminator="\n")
    temporary_paths[3].write_text(report, encoding="utf-8")

    final_paths = [TUMOUR_HEALTHY_PATH, SUBTYPE_PATH, GENE_SUMMARY_PATH, REPORT_PATH]
    for temporary_path, final_path in zip(temporary_paths, final_paths):
        temporary_path.replace(final_path)


def main() -> None:
    refuse_existing_outputs()
    print(f"Project: {PROJECT_DIR}")
    print(f"SciPy: {scipy.__version__}")
    print(f"Loading patient matrix: {MATRIX_PATH}")
    matrix, manifest, probe_ids = load_and_validate_inputs()
    print(
        f"Input validation passed: {matrix.shape[0]} probes × {matrix.shape[1]} observations; "
        f"{len(probe_ids)} targets; 130 tumour patients and 11 healthy controls."
    )
    print(f"Manifest columns: {', '.join(manifest.columns)}")

    tumour_healthy = run_tumour_vs_healthy(matrix, manifest, probe_ids)
    subtype_results = run_subtype_comparisons(matrix, manifest, probe_ids)
    gene_summary = make_gene_summary(tumour_healthy, subtype_results, matrix, manifest)
    report = build_report(tumour_healthy, subtype_results, gene_summary, manifest)
    save_outputs(tumour_healthy, subtype_results, gene_summary, report)

    print(
        "Tumour vs healthy: "
        f"{int(tumour_healthy['significant_bh_0_05'].sum())}/38 probes "
        "passed BH-adjusted p < 0.05."
    )
    print(
        "Subtype contrasts: "
        f"{int(subtype_results['significant_bh_0_05'].sum())}/228 tests "
        "passed family-wide BH-adjusted p < 0.05."
    )
    print("Subtype contrast group sizes:")
    for first, second in itertools.combinations(SUBTYPES, 2):
        count = int(
            (
                (subtype_results["first_subtype"] == first)
                & (subtype_results["second_subtype"] == second)
            ).sum()
        )
        print(
            f"  {first} - {second}: "
            f"n={EXPECTED_SUBTYPE_COUNTS[first]} vs {EXPECTED_SUBTYPE_COUNTS[second]}, "
            f"{int(subtype_results.loc[(subtype_results['first_subtype'] == first) & (subtype_results['second_subtype'] == second), 'significant_bh_0_05'].sum())}/{count} significant"
        )
    print("Generated files:")
    for path in [TUMOUR_HEALTHY_PATH, SUBTYPE_PATH, GENE_SUMMARY_PATH, REPORT_PATH]:
        print(f"  {path}")
    print("No samples were removed; no additional replicate averaging was performed.")


if __name__ == "__main__":
    try:
        main()
    except (AnalysisValidationError, FileNotFoundError) as error:
        print(f"Analysis validation error: {error}", file=sys.stderr)
        raise SystemExit(1) from error
