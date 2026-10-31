"""Build the OncoMap v2 final report and self-contained results dashboard."""

from __future__ import annotations

import html
import gzip
import json
import math
import re
import sys
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"
PROCESSED = ROOT / "data" / "processed"
OUTPUT_REPORT = REPORTS / "OncoMap_v2_Final_Research_Report.md"
OUTPUT_DASHBOARD = ROOT / "dashboard" / "index.html"
ALPHA = 0.05

INPUTS = {
    "source_matrix": ROOT / "data" / "raw" / "GSE65194_series_matrix.txt.gz",
    "source_manifest": ROOT / "data" / "metadata" / "GSE65194_sample_manifest.tsv",
    "patient_matrix": PROCESSED / "GSE65194_patient_level_matrix.tsv.gz",
    "patient_manifest": PROCESSED / "GSE65194_patient_level_manifest.tsv",
    "compatibility": REPORTS / "GSE65194_compatibility_audit.md",
    "qc": REPORTS / "GSE65194_QC_report.md",
    "replicate_audit": REPORTS / "GSE65194_replicate_resolution_audit.md",
    "aggregation": REPORTS / "GSE65194_patient_level_aggregation_report.md",
    "tumor_welch": REPORTS / "GSE65194_tumor_vs_healthy_welch.tsv",
    "tumor_sensitivity": REPORTS / "GSE65194_limma_repA_vs_repB.tsv",
    "tumor_limma_report": REPORTS / "GSE65194_limma_vs_welch_report.md",
    "subtype_welch": REPORTS / "GSE65194_subtype_pairwise_welch.tsv",
    "subtype_limma": REPORTS / "GSE65194_limma_subtype_pairwise.tsv",
    "subtype_limma_report": REPORTS / "GSE65194_limma_subtype_vs_welch_report.md",
    "replication": REPORTS / "GSE42568_vs_GSE65194_replication.tsv",
    "replication_report": REPORTS / "GSE42568_vs_GSE65194_replication_report.md",
    "matched_panel": REPORTS / "GSE42568_vs_GSE65194_matched_panel_BH.tsv",
    "matched_panel_report": REPORTS / "GSE42568_vs_GSE65194_matched_panel_BH_report.md",
    "limma_rep_a_full": PROCESSED / "GSE65194_limma_repA_full_table.tsv",
    "limma_rep_b_full": PROCESSED / "GSE65194_limma_repB_full_table.tsv",
    "target_method_script": ROOT / "scripts" / "03_differential_expression.py",
    "replicate_limma_script": ROOT / "scripts" / "04_limma_sensitivity_analysis.R",
    "subtype_limma_script": ROOT / "scripts" / "05_limma_subtype_analysis.R",
    "qc_script": ROOT / "scripts" / "01_load_and_qc.py",
    "aggregation_script": ROOT / "scripts" / "02_build_patient_level_matrix.py",
}

SUBTYPE_ORDER = ["TNBC", "Her2", "Luminal A", "Luminal B"]
EXPECTED_SUBTYPE_COUNTS = {"TNBC": 41, "Her2": 30, "Luminal A": 29, "Luminal B": 30}
EXPECTED_CONTRASTS = [
    "TNBC - Her2",
    "TNBC - Luminal A",
    "TNBC - Luminal B",
    "Her2 - Luminal A",
    "Her2 - Luminal B",
    "Luminal A - Luminal B",
]


class BuildValidationError(RuntimeError):
    """Raised when an existing analysis input cannot support a safe build."""


def fail(message: str) -> None:
    raise BuildValidationError(message)


def require_columns(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        fail(f"{label} is missing required column(s): {', '.join(missing)}")


def read_tsv(key: str) -> pd.DataFrame:
    path = INPUTS[key]
    try:
        return pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as exc:
        fail(f"Could not read {path.relative_to(ROOT)}: {exc}")


def read_text(key: str) -> str:
    path = INPUTS[key]
    try:
        return path.read_text(encoding="utf-8")
    except Exception as exc:
        fail(f"Could not read {path.relative_to(ROOT)}: {exc}")


def verify_unique(frame: pd.DataFrame, column: str, label: str) -> None:
    if frame[column].isna().any() or frame[column].astype(str).duplicated().any():
        fail(f"{label} must contain unique, non-missing values in {column}.")


def as_bool(series: pd.Series, label: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    mapped = series.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )
    if mapped.isna().any():
        fail(f"{label} contains values that are not valid booleans.")
    return mapped.astype(bool)


def finite_column(frame: pd.DataFrame, name: str, label: str) -> pd.Series:
    values = pd.to_numeric(frame[name], errors="coerce")
    if values.isna().any() or not values.map(math.isfinite).all():
        fail(f"{label}.{name} contains missing or non-finite numeric values.")
    return values


def format_p(value: float) -> str:
    if pd.isna(value):
        return "—"
    if value == 0:
        return "0"
    if value < 0.0001:
        return f"{value:.2e}"
    return f"{value:.4g}"


def format_num(value: float, digits: int = 4) -> str:
    if pd.isna(value):
        return "—"
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[Any]], alignments: list[str] | None = None) -> str:
    header = "| " + " | ".join(headers) + " |"
    if alignments is None:
        alignments = ["---"] * len(headers)
    separator = "| " + " | ".join(alignments) + " |"
    body = [
        "| " + " | ".join(str(item).replace("|", "\\|").replace("\n", " ") for item in row) + " |"
        for row in rows
    ]
    return "\n".join([header, separator, *body])


def parse_target_mapping(compatibility_text: str) -> tuple[list[str], dict[str, str]]:
    mapping: dict[str, str] = {}
    for line in compatibility_text.splitlines():
        match = re.match(r"^\|\s*(BRCA1|BRCA2|TP53|ERBB2|PIK3CA|ESR1|PTEN|EGFR)\s*\|\s*(.*?)\s*\|\s*\d+\s*\|$", line)
        if match:
            gene, probes_text = match.groups()
            mapping[gene] = [probe.strip() for probe in probes_text.split(",")]
    expected_order = ["BRCA1", "BRCA2", "TP53", "ERBB2", "PIK3CA", "ESR1", "PTEN", "EGFR"]
    if list(mapping) != expected_order:
        fail("Could not parse the eight-gene target mapping from the compatibility audit.")
    probe_to_gene = {probe: gene for gene, probes in mapping.items() for probe in probes}
    if len(probe_to_gene) != 38:
        fail(f"Compatibility audit target mapping contains {len(probe_to_gene)} unique probes, expected 38.")
    return list(probe_to_gene), probe_to_gene


def assert_mentions(text: str, patterns: list[str], label: str) -> None:
    for pattern in patterns:
        if not re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL):
            fail(f"{label} does not substantiate expected fact matching: {pattern}")


def load_and_validate() -> dict[str, Any]:
    existing_outputs = [path for path in [OUTPUT_REPORT, OUTPUT_DASHBOARD] if path.exists()]
    if existing_outputs:
        fail(
            "Refusing to overwrite existing output(s): "
            + ", ".join(str(path.relative_to(ROOT)) for path in existing_outputs)
        )

    missing = [path for path in INPUTS.values() if not path.is_file()]
    if missing:
        fail("Missing required input(s): " + ", ".join(str(path.relative_to(ROOT)) for path in missing))

    compatibility = read_text("compatibility")
    qc = read_text("qc")
    replicate_audit = read_text("replicate_audit")
    aggregation = read_text("aggregation")
    tumor_limma_report = read_text("tumor_limma_report")
    subtype_limma_report = read_text("subtype_limma_report")
    replication_report = read_text("replication_report")
    matched_report = read_text("matched_panel_report")
    target_method = INPUTS["target_method_script"].read_text(encoding="utf-8")
    replicate_method = INPUTS["replicate_limma_script"].read_text(encoding="utf-8")
    subtype_method = INPUTS["subtype_limma_script"].read_text(encoding="utf-8")

    assert_mentions(
        compatibility,
        [
            r"54,?673\s+(?:unique\s+)?probe rows\s*[×x]\s*178\s+sample columns",
            r"GPL570",
            r"technical replicates were subsequently averaged",
            r"PTENP1///PTEN",
        ],
        "Compatibility audit",
    )
    assert_mentions(
        replication_report,
        [r"104 cancer samples and 17 normal samples", r"exact-probe matching", r"38 exact shared probes"],
        "Cross-cohort replication report",
    )
    assert_mentions(qc, [r"54673\s*×\s*178", r"23\s*\(expected 23\)", r"130\s*\(expected 130\)"], "QC report")
    assert_mentions(
        replicate_audit,
        [r"exactly \*\*23\*\*", r"130 patients", r"178 GSM columns"],
        "Replicate audit",
    )
    assert_mentions(
        aggregation,
        [
            r"Deposited expression matrix\s*\|\s*54,?673\s*\|\s*178 GSM arrays",
            r"derived patient-level matrix",
            r"arithmetic mean",
            r"14 cell-line arrays",
            r"TNBC\s*\|\s*41",
            r"Her2\s*\|\s*30",
            r"Luminal A\s*\|\s*29",
            r"Luminal B\s*\|\s*30",
        ],
        "Patient-level aggregation report",
    )

    probe_ids, probe_to_gene = parse_target_mapping(compatibility)
    if any(not re.fullmatch(r"[0-9]+_[a-z0-9]+", probe) for probe in probe_ids):
        fail("Target annotations contain an unexpected probe identifier format.")
    gene_annotations = dict(probe_to_gene)
    gene_annotations["217492_s_at"] = "PTENP1 /// PTEN"
    assert_mentions(target_method, [r"TARGET_PROBES", r"ttest_ind", r"benjamini_hochberg"], "Welch analysis script")
    assert_mentions(replicate_method, [r"limma::lmFit", r"limma::eBayes", r"adjust.method\s*=\s*\"BH\""], "Tumour/healthy limma script")
    assert_mentions(subtype_method, [r"limma::makeContrasts", r"limma::eBayes", r"p.adjust\(limma_results\$p_value,\s*method\s*=\s*\"BH\"\)"], "Subtype limma script")

    replication = read_tsv("replication")
    require_columns(
        replication,
        [
            "probe_id", "gene_symbol", "GSE42568_gene_annotation", "GSE42568_logFC",
            "GSE42568_adjusted_p_value", "GSE42568_cancer_n", "GSE42568_normal_n",
            "GSE65194_tumour_minus_healthy_mean_difference", "GSE65194_adjusted_p_value",
            "GSE65194_tumour_n", "GSE65194_healthy_n", "direction_comparison",
            "GSE42568_significant_adj_p_lt_0_05", "GSE65194_significant_adj_p_lt_0_05",
            "significance_pattern",
        ],
        "Cross-cohort replication TSV",
    )
    if len(replication) != 38:
        fail(f"Cross-cohort replication TSV has {len(replication)} rows, expected 38.")
    verify_unique(replication, "probe_id", "Cross-cohort replication TSV")
    if set(replication["probe_id"]) != set(probe_ids):
        fail("Cross-cohort replication probe IDs do not exactly match the compatibility audit.")
    if not (replication["GSE42568_cancer_n"] == 104).all() or not (replication["GSE42568_normal_n"] == 17).all():
        fail("GSE42568 sample counts in replication TSV are not 104 cancer and 17 normal.")
    if not (replication["GSE65194_tumour_n"] == 130).all() or not (replication["GSE65194_healthy_n"] == 11).all():
        fail("GSE65194 sample counts in replication TSV are not 130 tumour and 11 healthy.")
    for column in [
        "GSE42568_logFC", "GSE42568_adjusted_p_value",
        "GSE65194_tumour_minus_healthy_mean_difference", "GSE65194_adjusted_p_value",
    ]:
        finite_column(replication, column, "Cross-cohort replication TSV")
    expected_direction = (
        replication["GSE42568_logFC"].gt(0)
        == replication["GSE65194_tumour_minus_healthy_mean_difference"].gt(0)
    )
    expected_direction &= replication["GSE42568_logFC"].ne(0)
    expected_direction &= replication["GSE65194_tumour_minus_healthy_mean_difference"].ne(0)
    reported_direction = replication["direction_comparison"].eq("concordant")
    if not expected_direction.equals(reported_direction):
        fail("Cross-cohort direction classifications do not match the observed effect signs.")

    matched = read_tsv("matched_panel")
    require_columns(
        matched,
        [
            "probe_id", "GSE42568_adjusted_p_value_original",
            "GSE42568_adjusted_p_value_matched_38_BH",
            "GSE65194_adjusted_p_value_original",
            "GSE65194_adjusted_p_value_matched_38_BH",
            "GSE42568_significant_matched_38_BH_0_05",
            "GSE65194_significant_matched_38_BH_0_05",
        ],
        "Matched-panel BH TSV",
    )
    if len(matched) != 38 or set(matched["probe_id"]) != set(probe_ids):
        fail("Matched-panel BH TSV must contain the same 38 exact target probes.")
    verify_unique(matched, "probe_id", "Matched-panel BH TSV")
    for cohort in ["GSE42568", "GSE65194"]:
        original_sig = as_bool(matched[f"{cohort}_significant_adj_p_lt_0_05"], f"{cohort} original significance")
        matched_sig = as_bool(matched[f"{cohort}_significant_matched_38_BH_0_05"], f"{cohort} matched-panel significance")
        matched_p = finite_column(matched, f"{cohort}_adjusted_p_value_matched_38_BH", "Matched-panel BH TSV")
        if not matched_sig.equals(matched_p.lt(ALPHA)):
            fail(f"{cohort} matched-panel significance flags disagree with adjusted p-values.")
        if int(matched_sig.sum()) != int(original_sig.sum()):
            fail(f"{cohort} significance counts differ between original and matched-panel BH values.")
    if not matched["probe_id"].equals(replication["probe_id"]):
        matched = matched.set_index("probe_id").loc[replication["probe_id"]].reset_index()

    tumor_welch = read_tsv("tumor_welch")
    tumor_sensitivity = read_tsv("tumor_sensitivity")
    require_columns(
        tumor_welch,
        ["probe_id", "gene_symbol", "gene_annotation", "mean_log2_difference_tumour_minus_healthy",
         "p_value_two_sided", "bh_adjusted_p_value", "significant_bh_0_05"],
        "Tumour-versus-healthy Welch TSV",
    )
    require_columns(
        tumor_sensitivity,
        ["probe_id", "gene_symbol", "logFC_repA", "logFC_repB",
         "logFC_difference_repA_minus_repB", "p_value_repA", "p_value_repB",
         "adj_p_value_repA", "adj_p_value_repB", "significant_repA_adj_p_lt_0_05",
         "significant_repB_adj_p_lt_0_05"],
        "Tumour-versus-healthy limma sensitivity TSV",
    )
    for frame, label in [(tumor_welch, "Welch tumor-versus-healthy"), (tumor_sensitivity, "Limma repA/repB")]:
        if len(frame) != 38 or set(frame["probe_id"]) != set(probe_ids):
            fail(f"{label} results must contain the same 38 target probes.")
        verify_unique(frame, "probe_id", label)

    limma_a = read_tsv("limma_rep_a_full")
    limma_b = read_tsv("limma_rep_b_full")
    for frame, label in [(limma_a, "Full limma repA"), (limma_b, "Full limma repB")]:
        require_columns(frame, ["probe_id", "logFC", "P.Value", "adj.P.Val"], label)
        if len(frame) != 54673:
            fail(f"{label} contains {len(frame)} probes, expected 54,673.")
        verify_unique(frame, "probe_id", label)
    for label, full, logfc_col, raw_col, adj_col in [
        ("repA", limma_a, "logFC_repA", "p_value_repA", "adj_p_value_repA"),
        ("repB", limma_b, "logFC_repB", "p_value_repB", "adj_p_value_repB"),
    ]:
        target_full = full.set_index("probe_id").loc[probe_ids]
        sensitivity_indexed = tumor_sensitivity.set_index("probe_id").loc[probe_ids]
        for field, column in [("logFC", logfc_col), ("P.Value", raw_col), ("adj.P.Val", adj_col)]:
            if not (target_full[field].to_numpy() == sensitivity_indexed[column].to_numpy()).all():
                if not (
                    (target_full[field] - sensitivity_indexed[column]).abs() < 1e-12
                ).all():
                    fail(f"Full limma {label} {field} values disagree with the target comparison TSV.")

    subtype_welch = read_tsv("subtype_welch")
    subtype_limma = read_tsv("subtype_limma")
    require_columns(
        subtype_welch,
        ["first_subtype", "second_subtype", "probe_id", "gene_symbol",
         "mean_log2_difference_first_minus_second", "p_value_two_sided",
         "bh_adjusted_p_value", "significant_bh_0_05"],
        "Subtype Welch TSV",
    )
    require_columns(
        subtype_limma,
        ["first_subtype", "second_subtype", "contrast", "probe_id", "gene_symbol",
         "logFC", "welch_log2_difference_first_minus_second",
         "logFC_difference_limma_minus_welch", "limma_p_value", "welch_p_value",
         "limma_bh_adjusted_p_value", "welch_bh_adjusted_p_value",
         "direction_agrees_with_welch", "limma_significant_bh_0_05",
         "welch_significant_bh_0_05", "significance_comparison"],
        "Subtype limma comparison TSV",
    )
    if len(subtype_welch) != 228 or len(subtype_limma) != 228:
        fail("Both subtype comparison TSVs must contain all 228 probe-by-contrast rows.")
    key_cols = ["first_subtype", "second_subtype", "probe_id"]
    if subtype_welch.duplicated(key_cols).any() or subtype_limma.duplicated(key_cols).any():
        fail("Subtype comparison TSVs contain duplicate contrast/probe rows.")
    subtype_welch["contrast"] = subtype_welch["first_subtype"] + " - " + subtype_welch["second_subtype"]
    for frame, label in [(subtype_welch, "Subtype Welch"), (subtype_limma, "Subtype limma")]:
        if set(frame["probe_id"]) != set(probe_ids) or set(frame["contrast"]) != set(EXPECTED_CONTRASTS):
            fail(f"{label} results do not contain all 38 probes for each expected contrast.")
        counts = frame.groupby("contrast").size()
        if not (counts == 38).all():
            fail(f"{label} results do not have exactly 38 rows per contrast.")
    expected_keys = set(
        zip(subtype_welch["first_subtype"], subtype_welch["second_subtype"], subtype_welch["probe_id"])
    )
    limma_keys = set(
        zip(subtype_limma["first_subtype"], subtype_limma["second_subtype"], subtype_limma["probe_id"])
    )
    if expected_keys != limma_keys:
        fail("Subtype limma and Welch files do not match by exact contrast and probe ID.")
    if set(subtype_limma["gene_symbol"]) != set(probe_to_gene.values()):
        fail("Subtype results do not represent all eight target genes.")
    if not subtype_limma["direction_agrees_with_welch"].map(lambda x: str(x).lower() == "true").all():
        fail("Subtype limma direction-agreement flags contain false values.")

    for frame, name in [(tumor_welch, "Welch tumor-versus-healthy"), (tumor_sensitivity, "Limma tumor sensitivity")]:
        if set(frame["gene_symbol"]) != set(probe_to_gene.values()):
            fail(f"{name} does not contain all eight target genes.")

    if not subtype_limma["contrast"].equals(subtype_welch["contrast"]):
        order = ["first_subtype", "second_subtype", "probe_id"]
        subtype_limma = subtype_limma.sort_values(order).reset_index(drop=True)
        subtype_welch = subtype_welch.sort_values(order).reset_index(drop=True)

    try:
        with gzip.open(INPUTS["source_matrix"], "rb") as source:
            if source.read(2) != b"!S":
                fail("The original GSE65194 matrix is not a readable gzip series-matrix file.")
        with gzip.open(INPUTS["patient_matrix"], "rb") as patient:
            if not patient.readline().startswith(b"ID_REF\t"):
                fail("The derived patient matrix is not a readable gzip TSV with ID_REF.")
    except (OSError, EOFError) as exc:
        fail(f"Could not validate compressed matrix input: {exc}")
    manifest = pd.read_csv(INPUTS["patient_manifest"], sep="\t")
    require_columns(manifest, ["sample_id", "sample_type", "subtype"], "Patient-level manifest")
    if len(manifest) != 141 or int((manifest["sample_type"] == "tumour").sum()) != 130:
        fail("Patient-level manifest does not contain 130 tumour patients and 11 controls.")
    manifest_subtypes = manifest.loc[manifest["sample_type"] == "tumour", "subtype"].value_counts().to_dict()
    if manifest_subtypes != EXPECTED_SUBTYPE_COUNTS:
        fail(f"Patient-level manifest subtype counts differ from expected: {manifest_subtypes}.")
    source_manifest = pd.read_csv(INPUTS["source_manifest"], sep="\t", keep_default_na=False)
    require_columns(
        source_manifest,
        ["gsm_id", "tumor_id", "sample_group", "is_cell_line"],
        "Source-array manifest",
    )
    if len(source_manifest) != 178 or source_manifest["gsm_id"].duplicated().any():
        fail("Source-array manifest must contain exactly 178 unique GSM records.")
    source_groups = source_manifest["sample_group"].value_counts().to_dict()
    if (
        source_groups.get("Healthy") != 11
        or source_groups.get("CellLine") != 14
        or sum(source_groups.get(group, 0) for group in EXPECTED_SUBTYPE_COUNTS) != 153
    ):
        fail(f"Source-array manifest composition is inconsistent: {source_groups}.")
    tumour_source = source_manifest[
        source_manifest["sample_group"].isin(EXPECTED_SUBTYPE_COUNTS)
    ]
    patient_array_counts = tumour_source.groupby("tumor_id").size()
    if (
        len(patient_array_counts) != 130
        or int((patient_array_counts == 2).sum()) != 23
        or int((patient_array_counts == 1).sum()) != 107
        or int((patient_array_counts > 2).sum()) != 0
    ):
        fail("Source-array manifest does not confirm 130 patients, 23 pairs and 107 singleton arrays.")
    source_subtype_patient_counts = {
        group: int(tumour_source.loc[tumour_source["sample_group"] == group, "tumor_id"].nunique())
        for group in EXPECTED_SUBTYPE_COUNTS
    }
    if source_subtype_patient_counts != EXPECTED_SUBTYPE_COUNTS:
        fail(f"Source-array patient counts by subtype differ from expected: {source_subtype_patient_counts}.")

    return {
        "probe_ids": probe_ids,
        "probe_to_gene": probe_to_gene,
        "gene_annotations": gene_annotations,
        "composition": {
            "subtypes": source_subtype_patient_counts,
            "healthy": int((source_manifest["sample_group"] == "Healthy").sum()),
            "cell_lines": int((source_manifest["sample_group"] == "CellLine").sum()),
            "tumour_patients": int(tumour_source["tumor_id"].nunique()),
            "tumour_arrays": len(tumour_source),
            "total_arrays": len(source_manifest),
            "replicate_pairs": int((patient_array_counts == 2).sum()),
            "single_array_patients": int((patient_array_counts == 1).sum()),
        },
        "replication": replication,
        "matched": matched,
        "tumor_welch": tumor_welch,
        "tumor_sensitivity": tumor_sensitivity,
        "subtype_welch": subtype_welch,
        "subtype_limma": subtype_limma,
        "limma_a": limma_a,
        "limma_b": limma_b,
        "texts": {
            "compatibility": compatibility,
            "qc": qc,
            "replicate_audit": replicate_audit,
            "aggregation": aggregation,
            "tumor_limma_report": tumor_limma_report,
            "subtype_limma_report": subtype_limma_report,
            "replication_report": replication_report,
            "matched_report": matched_report,
            "target_method_script": target_method,
            "replicate_method_script": replicate_method,
            "subtype_method_script": subtype_method,
        },
    }


def make_summaries(data: dict[str, Any]) -> dict[str, Any]:
    replication = data["replication"]
    matched = data["matched"]
    tumor_welch = data["tumor_welch"]
    sensitivity = data["tumor_sensitivity"]
    subtype = data["subtype_limma"]
    subtype_welch = data["subtype_welch"]
    replication_pattern = replication["significance_pattern"].value_counts()
    direction_counts = replication["direction_comparison"].value_counts()

    matched_a = as_bool(matched["GSE42568_significant_matched_38_BH_0_05"], "GSE42568 matched BH")
    matched_b = as_bool(matched["GSE65194_significant_matched_38_BH_0_05"], "GSE65194 matched BH")
    matched_categories = {
        "both": int((matched_a & matched_b).sum()),
        "GSE42568 only": int((matched_a & ~matched_b).sum()),
        "GSE65194 only": int((~matched_a & matched_b).sum()),
        "neither": int((~matched_a & ~matched_b).sum()),
    }
    subtype_summary: list[dict[str, Any]] = []
    for contrast in EXPECTED_CONTRASTS:
        rows = subtype[subtype["contrast"] == contrast]
        subtype_summary.append(
            {
                "contrast": contrast,
                "limma_sig": int(as_bool(rows["limma_significant_bh_0_05"], contrast).sum()),
                "welch_sig": int(as_bool(rows["welch_significant_bh_0_05"], contrast).sum()),
                "both": int((rows["significance_comparison"] == "Both").sum()),
                "limma_only": int((rows["significance_comparison"] == "Limma only").sum()),
                "welch_only": int((rows["significance_comparison"] == "Welch only").sum()),
                "neither": int((rows["significance_comparison"] == "Neither").sum()),
                "direction_agreement": int(as_bool(rows["direction_agrees_with_welch"], contrast).sum()),
            }
        )

    gene_summary = []
    for gene, probes in data["probe_to_gene"].items():
        rows = tumor_welch[tumor_welch["gene_symbol"] == gene]
        gene_summary.append(
            {
                "gene": gene,
                "probe_count": len(probes),
                "tumour_higher": int((rows["mean_log2_difference_tumour_minus_healthy"] > 0).sum()),
                "tumour_lower": int((rows["mean_log2_difference_tumour_minus_healthy"] < 0).sum()),
                "welch_significant": int(as_bool(rows["significant_bh_0_05"], gene).sum()),
                "limma_repA_significant": int(
                    as_bool(sensitivity[sensitivity["gene_symbol"] == gene]["significant_repA_adj_p_lt_0_05"], gene).sum()
                ),
                "limma_repB_significant": int(
                    as_bool(sensitivity[sensitivity["gene_symbol"] == gene]["significant_repB_adj_p_lt_0_05"], gene).sum()
                ),
            }
        )

    subtype_total = {
        "limma_sig": int(as_bool(subtype["limma_significant_bh_0_05"], "limma subtype").sum()),
        "welch_sig": int(as_bool(subtype["welch_significant_bh_0_05"], "Welch subtype").sum()),
        "both": int((subtype["significance_comparison"] == "Both").sum()),
        "limma_only": int((subtype["significance_comparison"] == "Limma only").sum()),
        "welch_only": int((subtype["significance_comparison"] == "Welch only").sum()),
        "neither": int((subtype["significance_comparison"] == "Neither").sum()),
        "direction_agreement": int(as_bool(subtype["direction_agrees_with_welch"], "direction agreement").sum()),
    }
    subtype_welch_sig = int(as_bool(subtype_welch["significant_bh_0_05"], "Welch subtype").sum())
    if subtype_welch_sig != subtype_total["welch_sig"]:
        fail("Welch adjusted-significance count differs between the source and comparison TSV.")

    return {
        "direction_counts": {
            "concordant": int(direction_counts.get("concordant", 0)),
            "discordant": int(direction_counts.get("discordant", 0)),
        },
        "replication_pattern": {str(k): int(v) for k, v in replication_pattern.items()},
        "matched_sig": {
            "GSE42568": int(matched_a.sum()),
            "GSE65194": int(matched_b.sum()),
        },
        "matched_categories": matched_categories,
        "tumor_welch_sig": int(as_bool(tumor_welch["significant_bh_0_05"], "Welch tumour healthy").sum()),
        "tumor_repA_sig": int(as_bool(sensitivity["significant_repA_adj_p_lt_0_05"], "limma repA").sum()),
        "tumor_repB_sig": int(as_bool(sensitivity["significant_repB_adj_p_lt_0_05"], "limma repB").sum()),
        "tumor_full_repA_sig": int((data["limma_a"]["adj.P.Val"] < ALPHA).sum()),
        "tumor_full_repB_sig": int((data["limma_b"]["adj.P.Val"] < ALPHA).sum()),
        "tumor_sensitivity": sensitivity,
        "subtype_summary": subtype_summary,
        "subtype_total": subtype_total,
        "gene_summary": gene_summary,
    }


def row_table_replication(data: dict[str, Any]) -> str:
    rows = []
    for _, row in data["replication"].iterrows():
        rows.append(
            [
                row["probe_id"], row["gene_symbol"], row["GSE42568_gene_annotation"],
                format_num(row["GSE42568_logFC"], 3), format_p(row["GSE42568_adjusted_p_value"]),
                format_num(row["GSE65194_tumour_minus_healthy_mean_difference"], 4),
                format_p(row["GSE65194_adjusted_p_value"]), row["direction_comparison"],
                "Yes" if row["significance_pattern"] == "significant_in_both" else "No",
            ]
        )
    return markdown_table(
        ["Probe", "Gene", "GSE42568 annotation", "GSE42568 logFC", "GSE42568 original BH p",
         "GSE65194 tumour−healthy difference", "GSE65194 BH p", "Direction", "Significant in both"],
        rows,
    )


def row_table_matched(data: dict[str, Any]) -> str:
    rows = []
    for _, row in data["matched"].iterrows():
        rows.append(
            [
                row["probe_id"], row["gene_symbol"],
                format_p(row["GSE42568_adjusted_p_value_matched_38_BH"]),
                format_p(row["GSE65194_adjusted_p_value_matched_38_BH"]),
                "Yes" if bool(row["GSE42568_significant_matched_38_BH_0_05"]) else "No",
                "Yes" if bool(row["GSE65194_significant_matched_38_BH_0_05"]) else "No",
            ]
        )
    return markdown_table(
        ["Probe", "Gene", "GSE42568 matched-38 BH p", "GSE65194 matched-38 BH p",
         "GSE42568 significant", "GSE65194 significant"],
        rows,
    )


def row_table_tumor_health(data: dict[str, Any]) -> str:
    welch = data["tumor_welch"].set_index("probe_id")
    sensitivity = data["tumor_sensitivity"].set_index("probe_id")
    rows = []
    for probe in data["probe_ids"]:
        w = welch.loc[probe]
        s = sensitivity.loc[probe]
        rows.append(
            [
                probe, w["gene_symbol"], format_num(w["mean_log2_difference_tumour_minus_healthy"], 4),
                format_p(w["bh_adjusted_p_value"]), format_num(s["logFC_repA"], 4),
                format_p(s["adj_p_value_repA"]), format_num(s["logFC_repB"], 4),
                format_p(s["adj_p_value_repB"]),
            ]
        )
    return markdown_table(
        ["Probe", "Gene", "Welch difference", "Welch BH p (38)", "limma repA logFC",
         "repA BH p (54,673)", "limma repB logFC", "repB BH p (54,673)"],
        rows,
    )


def row_table_subtypes(data: dict[str, Any], select_genes: set[str] | None = None) -> str:
    rows = []
    frame = data["subtype_limma"]
    if select_genes:
        frame = frame[frame["gene_symbol"].isin(select_genes)]
    for _, row in frame.iterrows():
        rows.append(
            [
                row["contrast"], row["probe_id"], row["gene_symbol"],
                format_num(row["logFC"], 4),
                format_num(row["welch_log2_difference_first_minus_second"], 4),
                format_p(row["limma_bh_adjusted_p_value"]), format_p(row["welch_bh_adjusted_p_value"]),
                "Yes" if bool(row["direction_agrees_with_welch"]) else "No",
                row["significance_comparison"],
            ]
        )
    return markdown_table(
        ["Contrast", "Probe", "Gene", "limma logFC", "Welch difference", "limma BH p (228)",
         "Welch BH p (228)", "Direction agrees", "Significance"],
        rows,
    )


def create_report(data: dict[str, Any], summary: dict[str, Any]) -> str:
    replication = data["replication"]
    concordant = summary["direction_counts"]["concordant"]
    discordant = summary["direction_counts"]["discordant"]
    pattern = summary["replication_pattern"]
    st = summary["subtype_total"]
    composition = data["composition"]
    full_limma_a = summary["tumor_full_repA_sig"]
    full_limma_b = summary["tumor_full_repB_sig"]
    gse42568_cancer_n = int(replication["GSE42568_cancer_n"].iloc[0])
    gse42568_normal_n = int(replication["GSE42568_normal_n"].iloc[0])
    gse65194_healthy_n = int(replication["GSE65194_healthy_n"].iloc[0])
    matrix_probe_count = len(data["limma_a"])
    limma_version_match = re.search(r"limma version:\s*([0-9.]+)", data["texts"]["tumor_limma_report"], re.I)
    limma_version = limma_version_match.group(1) if limma_version_match else "recorded in analysis report"

    abstract = (
        f"This exploratory analysis compares probe-level expression for eight breast-cancer-related target genes "
        f"({len(data['probe_ids'])} GPL570 probes) across GSE42568 and GSE65194, and examines tumour–healthy and subtype differences "
        f"in GSE65194. The source cohort results include {gse42568_cancer_n} cancer and {gse42568_normal_n} normal samples for GSE42568 and "
        f"{composition['tumour_patients']} tumour patients plus {gse65194_healthy_n} healthy controls in GSE65194. Across exact-probe cohort matches, "
        f"{concordant}/38 directions were concordant and {discordant}/38 discordant. In the matched-38 BH "
        f"comparison, {summary['matched_categories']['both']} probes were significant in both cohorts. "
        f"Patient-level subtype comparisons showed {st['direction_agreement']}/{len(data['subtype_limma'])} directional agreement "
        f"between Welch and limma; {st['limma_sig']} and {st['welch_sig']} tests, respectively, met adjusted "
        f"p < 0.05. These results are computational and exploratory, not clinical or diagnostic validation."
    )

    subtype_rows = [
        [
            row["contrast"], row["limma_sig"], row["welch_sig"], row["both"],
            row["limma_only"], row["welch_only"], row["neither"],
            f"{row['direction_agreement']}/38",
        ]
        for row in summary["subtype_summary"]
    ]
    subtype_summary_table = markdown_table(
        ["Contrast", "limma significant", "Welch significant", "Both", "limma only",
         "Welch only", "Neither", "Direction agrees"],
        subtype_rows,
    )
    gene_rows = [
        [
            row["gene"], row["probe_count"], row["tumour_higher"], row["tumour_lower"],
            row["welch_significant"], row["limma_repA_significant"], row["limma_repB_significant"],
        ]
        for row in summary["gene_summary"]
    ]
    gene_table = markdown_table(
        ["Target gene", "Probes", "Tumour higher (Welch)", "Tumour lower (Welch)",
         "Welch significant", "limma repA significant", "limma repB significant"],
        gene_rows,
    )

    sensitivity_rows = []
    for _, row in summary["tumor_sensitivity"].iterrows():
        sensitivity_rows.append(
            [
                row["probe_id"], row["gene_symbol"], format_num(row["logFC_repA"], 4),
                format_num(row["logFC_repB"], 4), format_num(row["logFC_difference_repA_minus_repB"], 4),
                format_p(row["adj_p_value_repA"]), format_p(row["adj_p_value_repB"]),
                "Yes" if bool(row["significant_repA_adj_p_lt_0_05"]) else "No",
                "Yes" if bool(row["significant_repB_adj_p_lt_0_05"]) else "No",
            ]
        )
    sensitivity_table = markdown_table(
        ["Probe", "Gene", "repA logFC", "repB logFC", "repA−repB", "repA BH p",
         "repB BH p", "repA significant", "repB significant"],
        sensitivity_rows,
    )

    subtype_direction_rows = []
    for gene in ["ESR1", "ERBB2"]:
        gene_data = data["subtype_limma"][data["subtype_limma"]["gene_symbol"] == gene]
        for contrast in EXPECTED_CONTRASTS:
            rows = gene_data[gene_data["contrast"] == contrast]
            subtype_direction_rows.append(
                [
                    gene, contrast,
                    f"{int((rows['logFC'] > 0).sum())} higher in first / {int((rows['logFC'] < 0).sum())} higher in second",
                    f"{int(as_bool(rows['limma_significant_bh_0_05'], gene).sum())}/"
                    f"{int(as_bool(rows['welch_significant_bh_0_05'], gene).sum())}",
                    f"{int(as_bool(rows['direction_agrees_with_welch'], gene).sum())}/{len(rows)}",
                ]
            )
    esr_erbb_summary = markdown_table(
        ["Gene", "Contrast", "Probe directions", "limma/Welch significant probes", "Direction agreement"],
        subtype_direction_rows,
    )
    erbb2_rows = data["subtype_limma"][data["subtype_limma"]["gene_symbol"] == "ERBB2"]
    erbb2_direction_text = "; ".join(
        f"{contrast}: {int((rows['logFC'] > 0).sum())}/3 higher in the first subtype"
        for contrast in EXPECTED_CONTRASTS
        for rows in [erbb2_rows[erbb2_rows["contrast"] == contrast]]
    )
    esr1_luminal = data["subtype_limma"][
        (data["subtype_limma"]["gene_symbol"] == "ESR1")
        & (data["subtype_limma"]["contrast"] == "Luminal A - Luminal B")
    ]
    esr1_luminal_limma_sig = int(as_bool(esr1_luminal["limma_significant_bh_0_05"], "ESR1 limma") .sum())
    esr1_luminal_welch_sig = int(as_bool(esr1_luminal["welch_significant_bh_0_05"], "ESR1 Welch").sum())

    target_probe_table = row_table_tumor_health(data)
    subtype_full_table = row_table_subtypes(data)
    result_files = [
        "GSE42568_vs_GSE65194_replication.tsv",
        "GSE42568_vs_GSE65194_matched_panel_BH.tsv",
        "GSE65194_tumor_vs_healthy_welch.tsv",
        "GSE65194_limma_repA_vs_repB.tsv",
        "GSE65194_subtype_pairwise_welch.tsv",
        "GSE65194_limma_subtype_pairwise.tsv",
        "../data/processed/GSE65194_limma_repA_full_table.tsv",
        "../data/processed/GSE65194_limma_repB_full_table.tsv",
    ]
    result_links = ", ".join(f"[`{Path(item).name}`]({item})" for item in result_files)

    return f"""# OncoMap v2 — Breast Cancer Expression and Molecular Subtype Explorer

**Final exploratory research report**<br>
**Data source:** GEO series GSE42568 and GSE65194<br>
**Scope:** Probe-level computational analysis of a prespecified eight-gene, 38-probe panel

## 1. Abstract

{abstract}

## 2. Introduction

Breast tumours show heterogeneous transcriptional profiles. Comparing expression between tumour and non-tumour breast tissue, and across clinically annotated molecular subtype groups, can help describe which probe-level patterns are present in public expression datasets. This project focuses on eight target genes—BRCA1, BRCA2, TP53, ERBB2, PIK3CA, ESR1, PTEN and EGFR—selected in the existing OncoMap target panel. Multiple microarray probes can represent a gene and may behave differently; therefore, all results remain at the exact-probe level rather than averaging probes into a synthetic gene-expression value.

An independent-cohort comparison is useful because a pattern that recurs in another dataset is less likely to depend entirely on one cohort's sample composition or analysis choices. It is not, by itself, biological validation: preprocessing, recruitment, sample sizes and statistical methods differ between cohorts, and the panel is targeted rather than genome-wide.

## 3. Research objectives

1. Describe how the 38 target probes differ between tumour and healthy breast tissue in GSE65194.
2. Describe probe-level differences across TNBC, Her2, Luminal A and Luminal B patient groups.
3. Assess whether tumour-versus-normal expression directions are consistent for exact shared probes in the independent GSE42568 cohort.

## 4. Dataset description

### GSE42568

The existing GSE42568 result source is a probe-level GEO2R/limma table from the GPL570 platform. The documented comparison is cancer tissue minus normal breast tissue, with **{gse42568_cancer_n} cancer and {gse42568_normal_n} normal samples**. Positive logFC indicates higher expression in cancer. Its source GEO2R adjusted p-values come from a broader probe-testing family. This report uses the existing result and does not reconstruct the cohort analysis.

### GSE65194

The downloaded GEO series matrix contains **{matrix_probe_count:,} probes × {composition['total_arrays']} GSM arrays** on GPL570. Array composition in the audit is {composition['tumour_arrays']} tumour arrays representing {composition['tumour_patients']} patients, {composition['healthy']} healthy breast controls and {composition['cell_lines']} cell-line arrays. The analysed patient-level cohort has **{composition['tumour_patients']} tumour patients** and {composition['healthy']} controls; subtype counts are {', '.join(f'{name} {count}' for name, count in composition['subtypes'].items())}. The cell-line arrays are excluded from primary tumour/healthy and patient-level subtype comparisons.

Both series identify GPL570 and report GC-RMA-related processing. GSE65194 also describes fixed-effect correction for array and hybridization batch; equivalent batch correction is not documented for GSE42568. GSE65194's exact CDF choice (standard AffyCDF or Brainarray) and explicit transformation history remain unresolved. The project analyses existing processed values as log2-scale expression and does not repeat transformation or normalization. The GEO records are available at [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568) and [GSE65194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65194).

## 5. Data preparation and quality control

The source-array matrix and manifest were checked for unique probe/GSM identifiers, numeric expression values and exact metadata matching. Preliminary QC found 54,673 unique probe rows and 178 unique GSM columns, no missing or non-numeric measurements, and no sample-median screen that was used to remove arrays. The screen was descriptive only; flagged arrays were retained.

The source manifest identifies {composition['tumour_arrays']} tumour arrays for {composition['tumour_patients']} tumour patients: **{composition['replicate_pairs']} technical replicate pairs** and {composition['single_array_patients']} singleton arrays. Pair-level correlations across {matrix_probe_count:,} probes were high, but correlation does not establish sample provenance or prove prior averaging. The deposited series matrix retains separate GSM columns for both members of all {composition['replicate_pairs']} pairs even though the GEO processing annotation states, “Samples with technical replicates were subsequently averaged.” The available records do not resolve where or whether that documented averaging occurred. Accordingly, the patient-level matrix is explicitly a **project-derived dataset**: for each validated pair, the two deposited processed values were arithmetically averaged probe-by-probe; single-array patients and healthy controls were retained unchanged. The original {composition['total_arrays']}-array matrix was not rewritten. Replicate-source uncertainty remains a limitation.

The derived matrix has {matrix_probe_count:,} probes × {composition['tumour_patients'] + composition['healthy']} observations. Patient-level tumour/healthy and subtype analyses use one tumour observation per patient; the {composition['healthy']} controls are used only for tumour-versus-healthy comparisons. The {composition['cell_lines']} cell-line arrays are not included in the derived matrix.

## 6. Statistical methodology

### Welch tumour-versus-healthy analysis

The existing patient-level Welch results use two-sided independent-samples Welch t-tests, without assuming equal group variances. Effect is the mean tumour minus healthy log2 expression; positive values indicate higher tumour expression. Benjamini–Hochberg (BH) adjustment was applied over the **38 target-probe tests**.

### limma tumour-versus-healthy sensitivity

Two local limma fits select one array per paired patient using alternative assignments: repA or repB, while retaining the same 107 singleton tumour samples and 11 healthy controls. These are **sensitivity analyses over alternative technical-replicate selections, not independent cohorts**. Both use the contrast Tumour − Healthy; positive logFC means higher tumour expression. Full limma tables retain all **54,673 probes**, and BH adjustment for each fit spans that full set. Hence adjusted-p significance for the 38 targets need not match Welch's 38-test family. The local fits use the downloaded matrix and planned sample assignments; they are not claimed to reproduce the GEO2R website exactly.

### Patient-level subtype analyses

The subtype analyses use the existing averaged patient-level matrix and {composition['tumour_patients']} tumour patients only. For each of {len(EXPECTED_CONTRASTS)} contrasts—TNBC − Her2, TNBC − Luminal A, TNBC − Luminal B, Her2 − Luminal A, Her2 − Luminal B and Luminal A − Luminal B—Welch tests and empirical-Bayes limma estimates are available for the same {len(data['probe_ids'])} probes. Positive effect means higher expression in the first-named subtype. BH adjustment is across the complete **{len(data['probe_ids']) * len(EXPECTED_CONTRASTS)}-test family ({len(data['probe_ids'])} probes × {len(EXPECTED_CONTRASTS)} contrasts)** for each method, not separately within each contrast. Limma uses a four-group design, `lmFit`, `makeContrasts`, `contrasts.fit` and `eBayes`.

### Independent-cohort matching and interpretation of effects

Cross-cohort comparisons use exact GPL570 probe IDs; gene labels are descriptive, not the join key. GSE42568 logFC and GSE65194 tumour-minus-healthy mean differences have the same stated orientation. Original adjusted p-values reflect the original method-specific families. A separate matched-panel file recalculates BH within the same **38 exact probes separately in each cohort**, which permits a like-for-like panel-size correction but does not harmonize preprocessing or study design. No cohort matrices were pooled.

## 7. Results

### Independent-cohort replication

Among the 38 exact shared probes, **{concordant}** had concordant and **{discordant}** had discordant tumour/cancer-minus-normal/healthy directions. Under the original cohort-specific adjustments, {pattern.get('significant_in_both', 0)} were significant in both, {pattern.get('GSE42568_only', 0)} only in GSE42568, {pattern.get('GSE65194_only', 0)} only in GSE65194, and {pattern.get('significant_in_neither', 0)} in neither. These are descriptive replication counts, not a pooled meta-analysis.

{row_table_replication(data)}

### Matched-panel BH comparison

Recalculation across exactly 38 probes yielded {summary['matched_sig']['GSE42568']} significant probes in GSE42568 and {summary['matched_sig']['GSE65194']} in GSE65194. Categories were {summary['matched_categories']['both']} significant in both, {summary['matched_categories']['GSE42568 only']} only in GSE42568, {summary['matched_categories']['GSE65194 only']} only in GSE65194 and {summary['matched_categories']['neither']} in neither. The original significance counts were unchanged by this matched-panel recalculation in these outputs. The two datasets still differ in cohort and preprocessing.

{row_table_matched(data)}

### Tumour versus healthy in GSE65194

The {len(data['probe_ids'])}-probe Welch analysis identified **{summary['tumor_welch_sig']}/{len(data['probe_ids'])}** target probes at BH-adjusted p < 0.05. Each full-probe limma fit identified **{full_limma_a:,}** (repA) and **{full_limma_b:,}** (repB) probes of {matrix_probe_count:,} at BH-adjusted p < 0.05; within the target panel, **{summary['tumor_repA_sig']}/{len(data['probe_ids'])}** and **{summary['tumor_repB_sig']}/{len(data['probe_ids'])}**, respectively, were significant. The limma–Welch sensitivity report records 38/38 directional agreement for each assignment. The different BH families ({len(data['probe_ids'])} versus {matrix_probe_count:,}) matter when comparing significance.

The 38 target probes by gene are listed below. The per-gene counts summarize individual probe directions/significance; they are not gene-level expression estimates.

{markdown_table(['Target gene', 'Probes', 'Tumour higher (Welch)', 'Tumour lower (Welch)', 'Welch significant', 'limma repA significant', 'limma repB significant'], gene_rows)}

{target_probe_table}

### limma repA/repB sensitivity

repA and repB are alternative selections from 23 paired patients. The probe-level comparison includes each logFC, difference, raw p-value and full-matrix BH p-value.

{sensitivity_table}

### Patient-level subtype comparisons

The patient-level Welch and limma files each contain all 228 probe-by-contrast comparisons. Overall, limma identified **{st['limma_sig']}/228** significant tests and Welch **{st['welch_sig']}/228**. **{st['both']}** were significant by both, **{st['limma_only']}** by limma only, **{st['welch_only']}** by Welch only and **{st['neither']}** by neither. Direction agreed for **{st['direction_agreement']}/228** tests. The table gives per-contrast counts.

{subtype_summary_table}

### ESR1 and ERBB2 subtype patterns

Across the six contrasts, the number of ERBB2 probes higher in the first-named subtype is: {erbb2_direction_text}. The exact-probe results and adjusted p-values below retain each probe separately. These group-level differences do not establish receptor status for individual samples.

ESR1 has nine probes with heterogeneous results rather than one interchangeable gene-level value. In the Luminal A − Luminal B contrast, {esr1_luminal_limma_sig} limma and {esr1_luminal_welch_sig} Welch ESR1 probes meet adjusted p < 0.05; for other contrasts, significance varies by probe. Exact signs and test counts are shown below without inferring clinical ER-positive status from Luminal labels.

{esr_erbb_summary}

#### ESR1 and ERBB2 probe-level subtype results

{row_table_subtypes(data, {'ESR1', 'ERBB2'})}

#### Complete 228-row subtype result

All target genes and probes, all six contrast directions, both methods' adjusted p-values, and significance categories are preserved in this table.

{subtype_full_table}

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

This analysis describes exact-probe tumour/healthy and subtype patterns in GSE65194 and compares tumour/normal directions with an independent GPL570 cohort. The outputs show {concordant}/38 cross-cohort direction concordance, {st['direction_agreement']}/228 subtype-direction agreement between Welch and limma, and clear probe-level heterogeneity that warrants preserving probe identities. These findings support exploratory follow-up only. Replication in additional cohorts, clarification of processing provenance, clinical annotation, and laboratory validation remain necessary before making biological or clinical claims.

## 11. Reproducibility

### Inputs and software

- GEO series: [GSE42568](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568) and [GSE65194](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65194).
- GSE65194 input: `data/raw/GSE65194_series_matrix.txt.gz`; derived inputs: `data/processed/GSE65194_patient_level_matrix.tsv.gz` and `data/processed/GSE65194_patient_level_manifest.tsv`.
- Existing cohort-level GSE42568 GEO2R/limma output and report are external to this v2 workspace; the v2 replication and matched-panel TSVs retain the probe-level values used here.
- Software recorded in source reports: R 4.6.2, limma {limma_version}; the final builder uses the existing `.venv` Python {sys.version.split()[0]} and pandas {pd.__version__}. No packages were installed and no scientific analysis was rerun for this report/dashboard build.

### Analysis scripts and commands

- `scripts/01_load_and_qc.py` — initial matrix/metadata QC.
- `scripts/02_build_patient_level_matrix.py` — derived patient-level aggregation.
- `scripts/03_differential_expression.py` — targeted Welch analyses and report.
- `scripts/04_limma_sensitivity_analysis.R` — alternative repA/repB full-probe limma fits.
- `scripts/05_limma_subtype_analysis.R` — patient-level limma subtype contrasts.
- Builder command: `\\.venv\\Scripts\\python.exe scripts\\06_build_final_report_dashboard.py`
- Limma commands previously used: `Rscript scripts/04_limma_sensitivity_analysis.R` and `Rscript scripts/05_limma_subtype_analysis.R` with R 4.6.2 and limma 3.68.5.

### Existing result files used

{result_links}

The standalone dashboard is [`../dashboard/index.html`](../dashboard/index.html). Its visualization data are embedded from the validated result tables; it uses no external scripts, fonts, APIs, CDNs or runtime data fetches.
"""


def json_records(frame: pd.DataFrame, columns: list[str] | None = None) -> list[dict[str, Any]]:
    if columns:
        frame = frame[columns]
    return json.loads(frame.to_json(orient="records", double_precision=15))


def dashboard_data(data: dict[str, Any], summary: dict[str, Any]) -> dict[str, Any]:
    replication = data["replication"].copy()
    matched = data["matched"].copy()
    sensitivity = data["tumor_sensitivity"].copy()
    tumor_welch = data["tumor_welch"].copy()
    subtype = data["subtype_limma"].copy()
    composition = data["composition"].copy()
    composition["subtypes"] = [
        {"name": name, "count": count} for name, count in composition["subtypes"].items()
    ]
    composition["gse42568_cancer"] = int(replication["GSE42568_cancer_n"].iloc[0])
    composition["gse42568_normal"] = int(replication["GSE42568_normal_n"].iloc[0])
    composition["gse65194_healthy"] = int(replication["GSE65194_healthy_n"].iloc[0])
    composition["probe_count"] = len(data["limma_a"])
    composition["gene_count"] = len(set(data["probe_to_gene"].values()))
    composition["target_probe_count"] = len(data["probe_ids"])
    composition["contrast_count"] = len(EXPECTED_CONTRASTS)
    composition["subtype_test_count"] = len(subtype)

    replicate_results = replication.merge(
        matched[
            [
                "probe_id", "GSE42568_adjusted_p_value_matched_38_BH",
                "GSE65194_adjusted_p_value_matched_38_BH",
                "GSE42568_significant_matched_38_BH_0_05",
                "GSE65194_significant_matched_38_BH_0_05",
            ]
        ],
        on="probe_id",
        validate="one_to_one",
    )
    tumor = tumor_welch.merge(sensitivity, on=["probe_id", "gene_symbol"], validate="one_to_one")
    subtype_records = subtype.copy()
    return {
        "genes": list(dict.fromkeys(data["probe_to_gene"].values())),
        "subtypes": SUBTYPE_ORDER,
        "composition": composition,
        "summary": {
            "direction": summary["direction_counts"],
            "replication_significance": summary["replication_pattern"],
            "matched_significance": summary["matched_categories"],
            "matched_sig_counts": summary["matched_sig"],
            "tumor_welch_sig": summary["tumor_welch_sig"],
            "tumor_repA_sig": summary["tumor_repA_sig"],
            "tumor_repB_sig": summary["tumor_repB_sig"],
            "full_limma_repA_sig": summary["tumor_full_repA_sig"],
            "full_limma_repB_sig": summary["tumor_full_repB_sig"],
            "subtype": summary["subtype_total"],
            "subtype_contrasts": summary["subtype_summary"],
        },
        "replication": json_records(
            replicate_results,
            [
                "probe_id", "gene_symbol", "GSE42568_gene_annotation",
                "GSE42568_logFC", "GSE42568_adjusted_p_value",
                "GSE65194_tumour_minus_healthy_mean_difference",
                "GSE65194_adjusted_p_value", "direction_comparison",
                "GSE42568_significant_adj_p_lt_0_05", "GSE65194_significant_adj_p_lt_0_05",
                "GSE42568_adjusted_p_value_matched_38_BH",
                "GSE65194_adjusted_p_value_matched_38_BH",
                "GSE42568_significant_matched_38_BH_0_05",
                "GSE65194_significant_matched_38_BH_0_05",
            ],
        ),
        "tumorHealthy": json_records(
            tumor,
            [
                "probe_id", "gene_symbol", "gene_annotation",
                "mean_log2_difference_tumour_minus_healthy", "p_value_two_sided",
                "bh_adjusted_p_value", "significant_bh_0_05",
                "logFC_repA", "logFC_repB", "logFC_difference_repA_minus_repB",
                "p_value_repA", "p_value_repB", "adj_p_value_repA", "adj_p_value_repB",
                "significant_repA_adj_p_lt_0_05", "significant_repB_adj_p_lt_0_05",
            ],
        ),
        "subtype": json_records(
            subtype_records,
            [
                "first_subtype", "second_subtype", "contrast", "probe_id", "gene_symbol",
                "gene_annotation", "logFC", "welch_log2_difference_first_minus_second",
                "logFC_difference_limma_minus_welch", "limma_p_value", "welch_p_value",
                "limma_bh_adjusted_p_value", "welch_bh_adjusted_p_value",
                "direction_agrees_with_welch", "limma_significant_bh_0_05",
                "welch_significant_bh_0_05", "significance_comparison",
            ],
        ),
        "links": {
            "replication": "../reports/GSE42568_vs_GSE65194_replication.tsv",
            "matched": "../reports/GSE42568_vs_GSE65194_matched_panel_BH.tsv",
            "tumorWelch": "../reports/GSE65194_tumor_vs_healthy_welch.tsv",
            "tumorLimma": "../reports/GSE65194_limma_repA_vs_repB.tsv",
            "subtypeWelch": "../reports/GSE65194_subtype_pairwise_welch.tsv",
            "subtypeLimma": "../reports/GSE65194_limma_subtype_pairwise.tsv",
            "fullRepA": "../data/processed/GSE65194_limma_repA_full_table.tsv",
            "fullRepB": "../data/processed/GSE65194_limma_repB_full_table.tsv",
        },
    }


def create_dashboard(dashboard: dict[str, Any]) -> str:
    payload = json.dumps(dashboard, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    payload = payload.replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    payload = payload.replace("\u2028", "\\u2028").replace("\u2029", "\\u2029")
    return r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="theme-color" content="#f5f8fa">
  <meta name="description" content="Exploratory, probe-level breast cancer expression research dashboard for GSE42568 and GSE65194.">
  <title>OncoMap v2 | Breast Cancer Expression Explorer</title>
  <style>
    :root{--navy:#102a43;--ink:#243b53;--muted:#627d98;--teal:#087f8c;--teal2:#dff5f3;--blue:#2878b5;--blue2:#e7f1f8;--line:#d9e2ec;--paper:#fff;--bg:#f4f7fa;--gold:#e5a839;--rose:#c85363;--green:#218567;--shadow:0 12px 32px rgba(16,42,67,.08);--radius:18px}
    *{box-sizing:border-box}html{scroll-behavior:smooth;scroll-padding-top:82px}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.58 "Segoe UI",Inter,system-ui,-apple-system,BlinkMacSystemFont,Arial,sans-serif}
    a{color:var(--teal);text-decoration:none}a:hover{text-decoration:underline}button,input,select{font:inherit}button{cursor:pointer}
    .topbar{position:sticky;top:0;z-index:5;background:rgba(255,255,255,.94);backdrop-filter:blur(14px);border-bottom:1px solid var(--line)}
    .top-inner{max-width:1440px;margin:auto;padding:12px 28px;display:flex;align-items:center;justify-content:space-between;gap:24px}
    .brand{font-weight:800;letter-spacing:-.03em;color:var(--navy);font-size:18px;white-space:nowrap}.brand span{color:var(--teal)}
    nav{display:flex;gap:18px;align-items:center;overflow:auto;white-space:nowrap}nav a{font-size:13px;font-weight:650;color:#486581}
    .page{max-width:1440px;margin:0 auto;padding:28px;display:grid;grid-template-columns:220px minmax(0,1fr);gap:28px}
    .rail{position:sticky;top:78px;align-self:start;padding:18px 0}.rail .eyebrow{font-size:11px;text-transform:uppercase;letter-spacing:.14em;color:var(--muted);font-weight:750;margin:8px 12px}
    .rail a{display:block;color:#486581;padding:9px 12px;border-radius:10px;font-size:13px;font-weight:650}.rail a:hover,.rail a.active{background:#e8f4f4;color:var(--teal);text-decoration:none}
    main{min-width:0}.hero{background:linear-gradient(116deg,#102a43 0%,#164e63 62%,#087f8c 100%);color:#fff;border-radius:24px;padding:42px 44px;box-shadow:var(--shadow);position:relative;overflow:hidden}
    .hero:after{content:"";position:absolute;width:300px;height:300px;border:1px solid rgba(255,255,255,.16);border-radius:50%;right:-55px;top:-155px;box-shadow:0 0 0 42px rgba(255,255,255,.04),0 0 0 88px rgba(255,255,255,.03)}
    .hero .kicker{color:#a7e7e1;text-transform:uppercase;letter-spacing:.15em;font-size:11px;font-weight:800}.hero h1{max-width:820px;font-size:clamp(32px,5vw,54px);line-height:1.04;letter-spacing:-.045em;margin:14px 0}.hero p{max-width:750px;color:#d9e8ef;font-size:17px}
    .hero-actions{display:flex;flex-wrap:wrap;gap:12px;margin-top:24px}.pill{display:inline-flex;align-items:center;border:1px solid rgba(255,255,255,.32);border-radius:99px;color:#fff;padding:8px 13px;font-size:13px;font-weight:650;background:rgba(255,255,255,.08)}.pill:hover{background:rgba(255,255,255,.15);text-decoration:none}
    .metrics{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px;margin:18px 0 28px}.metric,.card,.section-card{background:var(--paper);border:1px solid var(--line);border-radius:var(--radius);box-shadow:0 3px 14px rgba(16,42,67,.035)}
    .metric{padding:18px 20px}.metric strong{display:block;color:var(--navy);font-size:29px;letter-spacing:-.04em}.metric span{color:var(--muted);font-size:12px;font-weight:700}
    .section-card{padding:28px;margin:24px 0;scroll-margin-top:84px}.section-head{display:flex;justify-content:space-between;align-items:flex-start;gap:16px;margin-bottom:20px}.section-head h2{color:var(--navy);font-size:24px;letter-spacing:-.03em;margin:0}.section-head p{margin:6px 0 0;color:var(--muted);max-width:760px}.tag{display:inline-block;background:var(--teal2);color:#08717b;padding:5px 9px;border-radius:99px;font-size:11px;font-weight:800;white-space:nowrap}
    .grid-2{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:16px}.grid-3{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:14px}
    .card{padding:20px}.card h3{margin:0 0 7px;font-size:15px;color:var(--navy)}.card p{margin:6px 0;color:var(--muted);font-size:13px}.callout{border-left:4px solid var(--gold);background:#fff8e9;padding:14px 17px;border-radius:0 12px 12px 0;color:#72531b;font-size:13px;margin:16px 0}
    .chart{width:100%;min-height:180px}.chart svg{width:100%;height:auto;overflow:visible}.chart text{font-family:inherit;fill:#486581;font-size:12px}.chart .value{fill:var(--navy);font-size:13px;font-weight:800}
    .legend{display:flex;flex-wrap:wrap;gap:12px;color:var(--muted);font-size:12px}.legend i{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:5px}
    .table-tools{display:flex;flex-wrap:wrap;gap:10px;margin:14px 0}.table-tools input,.table-tools select{border:1px solid var(--line);border-radius:10px;background:#fff;padding:9px 11px;color:var(--ink);min-width:160px}.table-tools input{flex:1}
    .table-wrap{overflow:auto;border:1px solid var(--line);border-radius:12px;max-height:560px;background:#fff}.data-table{width:100%;border-collapse:collapse;font-size:12px;white-space:nowrap}.data-table th{position:sticky;top:0;background:#edf3f7;color:#486581;text-align:left;padding:10px 11px;font-size:11px;z-index:1}.data-table td{padding:9px 11px;border-top:1px solid #edf2f6}.data-table tbody tr:hover{background:#f3fbfa}.num{text-align:right;font-variant-numeric:tabular-nums}.badge{display:inline-block;border-radius:99px;padding:3px 8px;font-size:10px;font-weight:800}.yes{background:#e4f4ee;color:#176b53}.no{background:#f0f3f6;color:#627d98}.warn{background:#fff1d2;color:#8a5a00}.toolbar{display:flex;flex-wrap:wrap;gap:8px}.tab-button{border:1px solid var(--line);background:#fff;color:#486581;border-radius:99px;padding:8px 14px;font-weight:750;font-size:12px}.tab-button[aria-pressed="true"]{background:var(--navy);border-color:var(--navy);color:#fff}
    .result-links{display:flex;flex-wrap:wrap;gap:8px;margin-top:14px}.result-links a{border:1px solid var(--line);border-radius:9px;padding:7px 10px;font-size:11px;font-weight:700}
    .gene-cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px}.gene-card{border:1px solid var(--line);border-radius:12px;padding:12px;background:#fbfdfe}.gene-card strong{color:var(--navy)}.gene-card span{display:block;color:var(--muted);font-size:11px;margin-top:3px}
    .notice{color:var(--muted);font-size:12px}.footer{padding:18px 4px 38px;color:var(--muted);font-size:12px}.footer strong{color:var(--navy)}
    [hidden]{display:none!important}
    @media(max-width:1050px){.page{grid-template-columns:1fr;padding:20px}.rail{display:none}.metrics{grid-template-columns:repeat(2,minmax(0,1fr))}.top-inner{padding:10px 18px}.top-inner nav{display:none}}
    @media(max-width:650px){.page{padding:12px}.hero{padding:29px 23px;border-radius:18px}.section-card{padding:19px 15px}.grid-2,.grid-3{grid-template-columns:1fr}.metrics{gap:9px}.metric{padding:13px}.metric strong{font-size:24px}.gene-cards{grid-template-columns:repeat(2,minmax(0,1fr))}.section-head{display:block}.tag{margin-top:10px}.table-tools input,.table-tools select{min-width:100%}}
  </style>
</head>
<body>
<header class="topbar"><div class="top-inner">
  <a class="brand" href="#overview" aria-label="OncoMap home">Onco<span>Map</span> v2</a>
  <nav aria-label="Primary navigation">
    <a href="#overview">Overview</a><a href="#composition">Cohorts</a><a href="#replication">Replication</a>
    <a href="#tumor-healthy">Tumour / healthy</a><a href="#subtypes">Subtypes</a><a href="#methods">Methods</a>
  </nav>
</div></header>
<div class="page">
  <aside class="rail" aria-label="Section navigation">
    <div class="eyebrow">Explore</div>
    <a href="#overview">Overview</a><a href="#composition">Dataset composition</a><a href="#replication">Cross-cohort replication</a>
    <a href="#tumor-healthy">Tumour versus healthy</a><a href="#subtypes">Subtype explorer</a><a href="#methods">Methodology &amp; limits</a>
  </aside>
  <main>
    <section class="hero" id="overview">
      <div class="kicker">Exploratory public-data research · GPL570</div>
      <h1>Breast cancer expression &amp; molecular subtype explorer</h1>
      <p id="cohortSummary">Probe-level comparison of eight target genes across two GEO cohorts, with patient-balanced tumour/healthy and subtype analyses in GSE65194.</p>
      <div class="hero-actions">
        <a class="pill" href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568" target="_blank" rel="noopener">GSE42568 ↗</a>
        <a class="pill" href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65194" target="_blank" rel="noopener">GSE65194 ↗</a>
        <a class="pill" href="#replication">Explore probe results ↓</a>
      </div>
    </section>
    <div class="metrics" id="overviewMetrics"></div>
    <section class="section-card">
      <div class="section-head"><div><h2>Key findings</h2><p>All counts are derived from the validated project result tables; probe signals are not collapsed to gene-level values.</p></div><span class="tag">Exploratory, non-diagnostic</span></div>
      <div class="grid-3" id="keyFindings"></div>
    </section>
    <section class="section-card" id="composition">
      <div class="section-head"><div><h2>Dataset composition</h2><p>GSE65194 contains array-level source data and a project-derived patient-level matrix for primary tumour analyses.</p></div><span class="tag" id="patientCountTag"></span></div>
      <div class="grid-2">
        <div class="card"><h3>Patient-level tumour subtypes</h3><div class="chart" id="subtypeChart" role="img" aria-label="Bar chart of patient counts by tumour subtype"></div><div class="legend"><span>Each bar is one patient count from the patient-level manifest.</span></div></div>
        <div class="card"><h3>Source-array accounting</h3><div class="chart" id="compositionChart" role="img" aria-label="Bar chart of sample and replicate counts"></div><p id="replicateCompositionText"></p></div>
      </div>
      <div class="callout"><strong>Provenance note.</strong> GEO processing text says technical replicates were averaged, yet the downloaded matrix retains separate GSM columns for all 23 pairs. Where or whether that earlier averaging occurred is unresolved; this dashboard treats the project matrix as derived from the deposited array-level values.</div>
    </section>
    <section class="section-card" id="replication">
      <div class="section-head"><div><h2>Cross-cohort replication</h2><p>Exact GPL570 probe IDs join the cancer-minus-normal GSE42568 results to tumour-minus-healthy GSE65194 results.</p></div><span class="tag" id="probeCountTag"></span></div>
      <div class="grid-2">
        <div class="card"><h3>Effect direction</h3><div class="chart" id="directionChart" role="img" aria-label="Concordant and discordant probe directions"></div><p>Both cohort effects are oriented cancer/tumour minus normal/healthy.</p></div>
        <div class="card"><h3>Matched-panel BH significance</h3><div class="chart" id="matchedChart" role="img" aria-label="Matched 38-probe significance categories"></div><p>BH is recalculated independently within the same 38 probes for each cohort.</p></div>
      </div>
      <div class="result-links"><a id="linkReplication" href="../reports/GSE42568_vs_GSE65194_replication.tsv">Download cross-cohort TSV</a><a id="linkMatched" href="../reports/GSE42568_vs_GSE65194_matched_panel_BH.tsv">Download matched-panel TSV</a></div>
      <h3 style="margin-top:24px">Probe-level cohort comparison</h3>
      <div class="table-tools"><input id="repSearch" type="search" placeholder="Search gene or probe ID" aria-label="Search cross-cohort probe table"><select id="repGene" aria-label="Filter by target gene"><option value="">All genes</option></select><select id="repSort" aria-label="Sort cross-cohort table"><option value="probe_id">Sort: probe ID</option><option value="GSE65194_tumour_minus_healthy_mean_difference">Sort: GSE65194 effect</option><option value="GSE42568_logFC">Sort: GSE42568 effect</option><option value="GSE65194_adjusted_p_value">Sort: GSE65194 adjusted p</option><option value="GSE42568_adjusted_p_value">Sort: GSE42568 adjusted p</option></select></div>
      <div class="table-wrap"><table class="data-table"><thead><tr><th>Probe / gene</th><th>GSE42568 effect</th><th>GSE42568 original BH p</th><th>GSE65194 effect</th><th>GSE65194 BH p</th><th>Direction</th><th>Matched-panel status</th></tr></thead><tbody id="repBody"></tbody></table></div>
      <p class="notice">Original adjusted p-values use different testing families. The matched-panel chart and table use BH across the same 38 probes in each cohort.</p>
    </section>
    <section class="section-card" id="tumor-healthy">
      <div class="section-head"><div><h2>Tumour versus healthy</h2><p>Patient-level Welch results are compared with local limma repA/repB array-selection sensitivity fits.</p></div><span class="tag">38 target probes</span></div>
      <div class="grid-3" id="tumorSummary"></div>
      <div class="card" style="margin-top:16px"><h3>Interactive probe explorer</h3>
        <div class="toolbar" role="group" aria-label="Select analysis"><button class="tab-button" data-mode="tumor" aria-pressed="true">Tumour versus healthy</button><button class="tab-button" data-mode="subtype" aria-pressed="false">Subtype comparisons</button></div>
        <p id="explorerDescription" class="notice">Welch BH family: 38 probes. limma repA and repB BH family: all 54,673 probes in each full fit.</p>
        <div class="table-tools"><select id="explorerGene" aria-label="Filter explorer by gene"><option value="">All genes</option></select><input id="explorerSearch" type="search" placeholder="Search gene, probe or contrast" aria-label="Search interactive results"><select id="explorerSort" aria-label="Sort interactive results"><option value="probe_id">Sort: probe ID</option><option value="effect">Sort: effect</option><option value="adj">Sort: adjusted p-value</option></select><select id="contrastFilter" aria-label="Filter subtype contrast" hidden><option value="">All subtype contrasts</option></select></div>
        <div class="table-wrap"><table class="data-table"><thead id="explorerHead"></thead><tbody id="explorerBody"></tbody></table></div>
        <div class="result-links" id="explorerLinks"></div>
      </div>
    </section>
    <section class="section-card" id="subtypes">
      <div class="section-head"><div><h2>Subtype explorer</h2><p>Patient-level contrasts use one observation per tumour patient; sign is first subtype minus second subtype.</p></div><span class="tag" id="subtypeTestTag"></span></div>
      <div class="grid-2">
        <div class="card"><h3>Welch versus limma significant probes</h3><div class="chart" id="subtypeSigChart" role="img" aria-label="Comparison of significant probe counts by subtype contrast"></div><div class="legend"><span><i style="background:#087f8c"></i>limma</span><span><i style="background:#2878b5"></i>Welch</span></div></div>
        <div class="card"><h3>ESR1 and ERBB2 probe patterns</h3><p>Probe-level summaries across all six contrasts; sign is first subtype minus second subtype.</p><div id="genePatternCards" class="gene-cards"></div></div>
      </div>
      <h3 style="margin-top:24px">Featured ESR1 / ERBB2 probe results</h3>
      <p class="notice" id="featuredGeneNote"></p>
      <div class="table-wrap"><table class="data-table"><thead><tr><th>Contrast</th><th>Probe / gene</th><th>limma logFC</th><th>Welch difference</th><th>limma BH p</th><th>Welch BH p</th><th>Direction</th><th>Significance</th></tr></thead><tbody id="featuredGeneBody"></tbody></table></div>
      <div class="callout">Subtype names are preserved from the sample manifest. Luminal labels are not used to infer explicit clinical ER-positive status, and Her2 group labels are not a patient-level receptor assay.</div>
      <div class="result-links"><a href="../reports/GSE65194_subtype_pairwise_welch.tsv">Download subtype Welch TSV</a><a href="../reports/GSE65194_limma_subtype_pairwise.tsv">Download subtype limma/Welch TSV</a></div>
    </section>
    <section class="section-card" id="methods">
      <div class="section-head"><div><h2>Methodology and limitations</h2><p>Interpret the results as exploratory computational evidence from public microarray datasets.</p></div><span class="tag">No diagnostic use</span></div>
      <div class="grid-2">
        <div class="card"><h3>Pipeline</h3><p>Sample metadata were matched by GSM. Technical replicate pairs were verified and averaged probe-by-probe for the derived patient-level matrix. Welch independent-sample tests and empirical-Bayes limma fits use existing processed expression; no additional normalization or log transformation was performed.</p><p id="methodStatsText"></p></div>
        <div class="card"><h3>Limitations</h3><p id="limitationStatsText"></p><p>No laboratory, prospective clinical, treatment-response or diagnostic validation. Expression association does not establish causation, receptor status or clinical utility.</p></div>
      </div>
      <div class="result-links">
        <a href="../reports/GSE65194_limma_repA_vs_repB.tsv">Tumour/healthy limma sensitivity TSV</a>
        <a href="../data/processed/GSE65194_limma_repA_full_table.tsv">Full limma repA table</a>
        <a href="../data/processed/GSE65194_limma_repB_full_table.tsv">Full limma repB table</a>
      </div>
    </section>
    <footer class="footer"><strong>OncoMap v2.</strong> Exploratory research dashboard; not medical advice, diagnostic software, or a clinically validated biomarker. Data summaries are embedded from project TSV outputs. Open-source data: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568" target="_blank" rel="noopener">GSE42568</a> · <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65194" target="_blank" rel="noopener">GSE65194</a>.</footer>
  </main>
</div>
<script id="dashboard-data" type="application/json">__DATA__</script>
<script>
"use strict";
const D=JSON.parse(document.getElementById("dashboard-data").textContent);
const $=id=>document.getElementById(id);
const num=(x,d=3)=>Number.isFinite(Number(x))?Number(x).toFixed(d):"—";
const pval=x=>{x=Number(x);if(!Number.isFinite(x))return"—";if(x===0)return"0";return x<.0001?x.toExponential(2):x.toPrecision(3)};
const directionLabel=x=>Number(x)>0?"First subtype higher":Number(x)<0?"Second subtype higher":"No change";
const badge=(yes,label)=>`<span class="badge ${yes?"yes":"no"}">${label}</span>`;
function fillSelect(id){const e=$(id);D.genes.forEach(g=>{const o=document.createElement("option");o.value=g;o.textContent=g;e.appendChild(o)})}
["repGene","explorerGene"].forEach(fillSelect);
function svgBarChart(target,items,opts={}){const host=$(target),w=650,h=Math.max(150,items.length*38+32),left=opts.left||150,right=52,plot=w-left-right,row=34,max=Math.max(1,...items.map(x=>Math.max(0,Number(x.value)||0))),svg=document.createElementNS("http://www.w3.org/2000/svg","svg");svg.setAttribute("viewBox",`0 0 ${w} ${h}`);svg.setAttribute("role","img");items.forEach((item,i)=>{const y=16+i*row+6,barW=Math.max(1,plot*(Number(item.value)||0)/max),text=document.createElementNS(svg.namespaceURI,"text");text.setAttribute("x","0");text.setAttribute("y",String(y+13));text.textContent=item.label;svg.appendChild(text);const rect=document.createElementNS(svg.namespaceURI,"rect");rect.setAttribute("x",left);rect.setAttribute("y",String(y));rect.setAttribute("width",String(barW));rect.setAttribute("height","18");rect.setAttribute("rx","5");rect.setAttribute("fill",item.color||"#087f8c");svg.appendChild(rect);const val=document.createElementNS(svg.namespaceURI,"text");val.setAttribute("x",String(left+barW+8));val.setAttribute("y",String(y+14));val.setAttribute("class","value");val.textContent=String(item.value);svg.appendChild(val)});host.replaceChildren(svg)}
const dir=D.summary.direction,mc=D.summary.matched_significance,st=D.summary.subtype;
$("cohortSummary").textContent=`Probe-level comparison of ${D.composition.gene_count} target genes across GSE42568 (${D.composition.gse42568_cancer} cancer, ${D.composition.gse42568_normal} normal) and GSE65194 (${D.composition.tumour_patients} tumour patients, ${D.composition.gse65194_healthy} healthy controls), with patient-balanced subtype analyses.`;
$("patientCountTag").textContent=`${D.composition.tumour_patients} tumour patients`;
$("probeCountTag").textContent=`${D.composition.target_probe_count} exact probes`;
$("subtypeTestTag").textContent=`${D.composition.contrast_count} contrasts · ${D.composition.subtype_test_count} tests`;
$("replicateCompositionText").textContent=`${D.composition.tumour_arrays} tumour arrays represent ${D.composition.tumour_patients} patients; ${D.composition.replicate_pairs} validated pairs were arithmetically averaged probe-by-probe in the derived patient matrix. ${D.composition.single_array_patients} singleton patients retain their observed array.`;
$("featuredGeneNote").textContent=`All ${D.tumorHealthy.filter(x=>x.gene_symbol==="ESR1").length/6} ESR1 probes and ${D.tumorHealthy.filter(x=>x.gene_symbol==="ERBB2").length/6} ERBB2 probes appear separately for every contrast. Adjusted p-values use the ${D.composition.subtype_test_count}-test family.`;
$("methodStatsText").textContent=`Subtype BH adjustment spans ${D.composition.subtype_test_count} tests per method. Patient-level Welch tumour/healthy adjustment spans ${D.composition.target_probe_count} probes; each full limma tumour/healthy fit adjusts across ${D.composition.probe_count.toLocaleString()} probes. The separate matched-panel cohort comparison adjusts each cohort over ${D.composition.target_probe_count} exact probes.`;
$("limitationStatsText").textContent=`One independent validation cohort; targeted ${D.composition.target_probe_count}-probe panel; small GSE65194 healthy group (n=${D.composition.healthy}); different sample sizes, preprocessing histories and statistical methods; unresolved replicate-averaging provenance; and ambiguous annotation for 217492_s_at (PTENP1 /// PTEN).`;
$("overviewMetrics").innerHTML=[
  ["2","independent GEO cohorts"],[String(D.composition.gene_count),"target genes"],[String(D.composition.target_probe_count),"exact target probes"],[String(D.composition.tumour_patients),"GSE65194 tumour patients"]
].map(x=>`<div class="metric"><strong>${x[0]}</strong><span>${x[1]}</span></div>`).join("");
$("keyFindings").innerHTML=[
  [`${dir.concordant}/38`,`cross-cohort directions concordant`],
  [`${st.direction_agreement}/228`,`subtype-test directions agree: limma vs Welch`],
  [`${st.both}`,`subtype probe tests significant in both methods`]
].map(x=>`<div class="card"><h3>${x[0]}</h3><p>${x[1]}</p></div>`).join("");
svgBarChart("subtypeChart",D.composition.subtypes.map((x,i)=>({label:x.name,value:x.count,color:["#087f8c","#2878b5","#3b91a8","#69a8b5"][i]})));
svgBarChart("compositionChart",[
 {label:"Tumour patients",value:D.composition.tumour_patients,color:"#087f8c"},
 {label:"Tumour arrays",value:D.composition.tumour_arrays,color:"#2878b5"},
 {label:"Healthy controls",value:D.composition.healthy,color:"#69a8b5"},
 {label:"Cell-line arrays",value:D.composition.cell_lines,color:"#e5a839"},
 {label:"Replicate pairs",value:D.composition.replicate_pairs,color:"#9a7bc1"}
],{left:145});
svgBarChart("directionChart",[
 {label:"Concordant",value:dir.concordant,color:"#087f8c"},
 {label:"Discordant",value:dir.discordant,color:"#c85363"}
]);
svgBarChart("matchedChart",Object.entries(mc).map(([label,value],i)=>({label,value,color:["#218567","#2878b5","#e5a839","#aab8c5"][i]})),{left:140});
const repSig=D.summary.replication_significance;
const tv=D.summary;
$("tumorSummary").innerHTML=[
 ["Welch target probes",`${tv.tumor_welch_sig}/38 significant`],
 ["limma repA targets",`${tv.tumor_repA_sig}/38 significant`],
 ["limma repB targets",`${tv.tumor_repB_sig}/38 significant`]
].map(x=>`<div class="card"><h3>${x[0]}</h3><p><strong style="font-size:22px;color:var(--navy)">${x[1]}</strong></p><p>limma full-probe significant counts: repA ${tv.full_limma_repA_sig.toLocaleString()} / repB ${tv.full_limma_repB_sig.toLocaleString()} of 54,673.</p></div>`).join("");
function currentMode(){return document.querySelector('.tab-button[aria-pressed="true"]').dataset.mode}
function setMode(mode){document.querySelectorAll(".tab-button").forEach(b=>b.setAttribute("aria-pressed",String(b.dataset.mode===mode)));$("contrastFilter").hidden=mode!=="subtype";$("explorerDescription").textContent=mode==="tumor"?`Welch BH family: ${D.composition.target_probe_count} probes. limma repA and repB BH family: all ${D.composition.probe_count.toLocaleString()} probes in each full fit.`:`Patient-level subtype results: BH adjustment across the complete ${D.composition.subtype_test_count}-test family for each method.`;renderExplorer()}
document.querySelectorAll(".tab-button").forEach(b=>b.addEventListener("click",()=>setMode(b.dataset.mode)));
const contrastSelect=$("contrastFilter");D.subtypes;D.summary.subtype_contrasts.forEach(x=>{const o=document.createElement("option");o.value=x.contrast;o.textContent=x.contrast;contrastSelect.appendChild(o)});
function renderExplorer(){const mode=currentMode(),gene=$("explorerGene").value,search=$("explorerSearch").value.trim().toLowerCase(),sort=$("explorerSort").value;let rows,heads;
 if(mode==="tumor"){rows=D.tumorHealthy.filter(x=>(!gene||x.gene_symbol===gene)&&(!search||`${x.gene_symbol} ${x.probe_id}`.toLowerCase().includes(search)));heads=["Probe / gene","Welch difference","Welch BH p (38)","limma repA logFC","repA BH p (54,673)","limma repB logFC","repB BH p (54,673)"];rows.sort((a,b)=>sort==="effect"?b.mean_log2_difference_tumour_minus_healthy-a.mean_log2_difference_tumour_minus_healthy:sort==="adj"?a.bh_adjusted_p_value-b.bh_adjusted_p_value:a.probe_id.localeCompare(b.probe_id));$("explorerLinks").innerHTML=`<a href="${D.links.tumorWelch}">Welch TSV</a><a href="${D.links.tumorLimma}">limma sensitivity TSV</a><a href="${D.links.fullRepA}">full repA TSV</a><a href="${D.links.fullRepB}">full repB TSV</a>`;$("explorerHead").innerHTML=`<tr>${heads.map(x=>`<th>${x}</th>`).join("")}</tr>`;$("explorerBody").innerHTML=rows.map(x=>`<tr><td><strong>${x.probe_id}</strong><br><span class="notice">${x.gene_symbol}</span></td><td class="num">${num(x.mean_log2_difference_tumour_minus_healthy,4)}</td><td class="num">${pval(x.bh_adjusted_p_value)} ${badge(x.significant_bh_0_05,"Welch")}</td><td class="num">${num(x.logFC_repA,4)}</td><td class="num">${pval(x.adj_p_value_repA)} ${badge(x.significant_repA_adj_p_lt_0_05,"repA")}</td><td class="num">${num(x.logFC_repB,4)}</td><td class="num">${pval(x.adj_p_value_repB)} ${badge(x.significant_repB_adj_p_lt_0_05,"repB")}</td></tr>`).join("")}
 else{const contrast=$("contrastFilter").value;rows=D.subtype.filter(x=>(!gene||x.gene_symbol===gene)&&(!contrast||x.contrast===contrast)&&(!search||`${x.gene_symbol} ${x.probe_id} ${x.contrast}`.toLowerCase().includes(search)));rows.sort((a,b)=>sort==="effect"?b.logFC-a.logFC:sort==="adj"?a.limma_bh_adjusted_p_value-b.limma_bh_adjusted_p_value:a.probe_id.localeCompare(b.probe_id));heads=["Contrast","Probe / gene","limma logFC","Welch difference","limma BH p","Welch BH p","Direction","Significance"];$("explorerLinks").innerHTML=`<a href="${D.links.subtypeWelch}">Welch TSV</a><a href="${D.links.subtypeLimma}">limma subtype TSV</a>`;$("explorerHead").innerHTML=`<tr>${heads.map(x=>`<th>${x}</th>`).join("")}</tr>`;$("explorerBody").innerHTML=rows.map(x=>`<tr><td>${x.contrast}</td><td><strong>${x.probe_id}</strong><br><span class="notice">${x.gene_symbol}</span></td><td class="num">${num(x.logFC,4)}</td><td class="num">${num(x.welch_log2_difference_first_minus_second,4)}</td><td class="num">${pval(x.limma_bh_adjusted_p_value)} ${badge(x.limma_significant_bh_0_05,"limma")}</td><td class="num">${pval(x.welch_bh_adjusted_p_value)} ${badge(x.welch_significant_bh_0_05,"Welch")}</td><td>${directionLabel(x.logFC)}</td><td>${x.significance_comparison}</td></tr>`).join("")}
 if(!rows.length)$("explorerBody").innerHTML=`<tr><td colspan="${heads.length}">No probes match these filters.</td></tr>`;
}
["explorerGene","explorerSearch","explorerSort","contrastFilter"].forEach(id=>$(id).addEventListener("input",renderExplorer));
function renderReplication(){let gene=$("repGene").value,search=$("repSearch").value.trim().toLowerCase(),sort=$("repSort").value,rows=D.replication.filter(x=>(!gene||x.gene_symbol===gene)&&(!search||`${x.gene_symbol} ${x.probe_id}`.toLowerCase().includes(search)));rows.sort((a,b)=>typeof a[sort]==="string"?a[sort].localeCompare(b[sort]):Number(a[sort])-Number(b[sort]));$("repBody").innerHTML=rows.map(x=>{const matched=x.GSE42568_significant_matched_38_BH_0_05&&x.GSE65194_significant_matched_38_BH_0_05?"Both":x.GSE42568_significant_matched_38_BH_0_05?"GSE42568 only":x.GSE65194_significant_matched_38_BH_0_05?"GSE65194 only":"Neither";return `<tr><td><strong>${x.probe_id}</strong><br><span class="notice">${x.gene_symbol}</span></td><td class="num" title="Positive means higher cancer expression">${num(x.GSE42568_logFC,4)}</td><td class="num">${pval(x.GSE42568_adjusted_p_value)}</td><td class="num" title="Positive means higher tumour expression">${num(x.GSE65194_tumour_minus_healthy_mean_difference,4)}</td><td class="num">${pval(x.GSE65194_adjusted_p_value)}</td><td>${x.direction_comparison==="concordant"?badge(true,"Concordant"):badge(false,"Discordant")}</td><td>${matched}</td></tr>`}).join("")||`<tr><td colspan="7">No probes match these filters.</td></tr>`}
["repGene","repSearch","repSort"].forEach(id=>$(id).addEventListener("input",renderReplication));
svgBarChart("subtypeSigChart",D.summary.subtype_contrasts.flatMap((x,i)=>[{label:`${x.contrast} · limma`,value:x.limma_sig,color:"#087f8c"},{label:`${x.contrast} · Welch`,value:x.welch_sig,color:"#2878b5"}]),{left:180});
const esr=D.subtype.filter(x=>x.gene_symbol==="ESR1"),erbb=D.subtype.filter(x=>x.gene_symbol==="ERBB2");
$("genePatternCards").innerHTML=[["ERBB2",erbb],["ESR1",esr]].map(([gene,rows])=>{const sig=rows.filter(x=>x.limma_significant_bh_0_05).length;return `<div class="gene-card"><strong>${gene}</strong><span>${rows.length/6} probes · ${sig} limma-significant probe tests / ${rows.length} contrasts</span></div>`}).join("");
const featured=D.subtype.filter(x=>x.gene_symbol==="ESR1"||x.gene_symbol==="ERBB2");
$("featuredGeneBody").innerHTML=featured.map(x=>`<tr><td>${x.contrast}</td><td><strong>${x.probe_id}</strong><br><span class="notice">${x.gene_symbol}</span></td><td class="num">${num(x.logFC,4)}</td><td class="num">${num(x.welch_log2_difference_first_minus_second,4)}</td><td class="num">${pval(x.limma_bh_adjusted_p_value)}</td><td class="num">${pval(x.welch_bh_adjusted_p_value)}</td><td>${directionLabel(x.logFC)}</td><td>${x.significance_comparison}</td></tr>`).join("");
renderReplication();renderExplorer();
document.querySelectorAll(".rail a").forEach(a=>a.addEventListener("click",()=>{document.querySelectorAll(".rail a").forEach(x=>x.classList.remove("active"));a.classList.add("active")}));
</script>
</body>
</html>""".replace("__DATA__", payload)


def verify_dashboard(html_text: str, data: dict[str, Any]) -> None:
    if not html_text.startswith("<!doctype html>"):
        fail("Generated dashboard is not an HTML5 document.")
    if re.search(r"<(?:script|link)[^>]+(?:src|href)=[\"']https?://", html_text, re.I):
        fail("Dashboard contains an external script or stylesheet dependency.")
    if re.search(r"\bfetch\s*\(", html_text) or re.search(r"\bXMLHttpRequest\b", html_text):
        fail("Dashboard contains a network data-fetch dependency.")
    if html_text.count('id="dashboard-data"') != 1 or "__DATA__" in html_text:
        fail("Embedded dashboard data is missing or not safely substituted.")
    expected_probe_count = len(data["probe_ids"])
    if expected_probe_count != 38:
        fail("Dashboard target-probe source is not exactly 38 probes.")
    if len(data["subtype_limma"]) != 228:
        fail("Dashboard subtype data is not exactly 228 comparisons.")
    for link in [
        "../reports/GSE42568_vs_GSE65194_replication.tsv",
        "../reports/GSE65194_limma_subtype_pairwise.tsv",
        "../data/processed/GSE65194_limma_repA_full_table.tsv",
        "../data/processed/GSE65194_limma_repB_full_table.tsv",
    ]:
        if f'href="{link}"' not in html_text:
            fail(f"Dashboard is missing expected relative result link: {link}")
    if "https://fonts.googleapis.com" in html_text or "cdn" in html_text.lower():
        fail("Dashboard unexpectedly references a CDN or external font.")
    class DashboardParser(HTMLParser):
        def __init__(self) -> None:
            super().__init__()
            self.ids: list[str] = []
            self.hrefs: list[str] = []
            self.sources: list[str] = []

        def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
            attributes = dict(attrs)
            if attributes.get("id"):
                self.ids.append(str(attributes["id"]))
            if attributes.get("href"):
                self.hrefs.append(str(attributes["href"]))
            if attributes.get("src"):
                self.sources.append(str(attributes["src"]))

    parser = DashboardParser()
    parser.feed(html_text)
    parser.close()
    if len(parser.ids) != len(set(parser.ids)):
        fail("Dashboard contains duplicate HTML IDs.")
    ids = set(parser.ids)
    for href in parser.hrefs:
        if href.startswith("#"):
            if href[1:] not in ids:
                fail(f"Dashboard contains a broken in-page reference: {href}")
        elif not href.startswith(("http://", "https://", "mailto:")):
            target = (OUTPUT_DASHBOARD.parent / href.split("?", 1)[0].split("#", 1)[0]).resolve()
            if not target.is_file():
                fail(f"Dashboard contains a broken local link: {href}")
    if parser.sources:
        fail("Dashboard contains external or local runtime source dependencies.")


def main() -> int:
    try:
        data = load_and_validate()
        summary = make_summaries(data)
        report = create_report(data, summary)
        dashboard = create_dashboard(dashboard_data(data, summary))
        verify_dashboard(dashboard, data)
        OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT_DASHBOARD.parent.mkdir(parents=True, exist_ok=True)
        if OUTPUT_REPORT.exists() or OUTPUT_DASHBOARD.exists():
            fail("An output appeared during generation; refusing to overwrite it.")
        OUTPUT_REPORT.write_text(report, encoding="utf-8", newline="\n")
        OUTPUT_DASHBOARD.write_text(dashboard, encoding="utf-8", newline="\n")
        print("Input validation passed.")
        print(f"Validated 38 exact target probes across 8 genes.")
        print(f"Validated 54,673 probes in each full limma table.")
        print(f"Validated 228 subtype rows across all six contrasts.")
        print(
            "Cross-cohort direction: "
            f"{summary['direction_counts']['concordant']} concordant, "
            f"{summary['direction_counts']['discordant']} discordant."
        )
        print(
            "Subtype BH significance: "
            f"limma {summary['subtype_total']['limma_sig']}/228; "
            f"Welch {summary['subtype_total']['welch_sig']}/228; "
            f"direction agreement {summary['subtype_total']['direction_agreement']}/228."
        )
        for path in [OUTPUT_REPORT, OUTPUT_DASHBOARD]:
            print(f"Created: {path} ({path.stat().st_size:,} bytes)")
        return 0
    except BuildValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except Exception as exc:
        print(f"ERROR: final report/dashboard build failed: {exc}", file=sys.stderr)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
