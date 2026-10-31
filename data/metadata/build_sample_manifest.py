from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent.parent
MATRIX_PATH = BASE_DIR / "raw" / "GSE65194_series_matrix.txt.gz"
CLINICAL_PATH = BASE_DIR / "metadata" / "GSE65194_clinical_data_update.txt.gz"
MANIFEST_PATH = BASE_DIR / "metadata" / "GSE65194_sample_manifest.tsv"


def read_matrix_annotations() -> tuple[list[str], list[str], list[str], dict[str, dict[str, str]]]:
    with gzip.open(MATRIX_PATH, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    def get_values(prefix: str) -> list[str]:
        for line in lines:
            if line.startswith(prefix):
                return [value.strip('"') for value in line.rstrip("\n").split("\t")[1:]]
        raise ValueError(f"Missing annotation block: {prefix!r}")

    sample_titles = get_values("!Sample_title\t")
    sample_accessions = get_values("!Sample_geo_accession\t")
    source_names = get_values("!Sample_source_name_ch1\t")

    sample_attributes: dict[str, dict[str, str]] = {}
    for line in lines:
        if not line.startswith("!Sample_characteristics_ch1"):
            continue
        parts = [part.strip('"') for part in line.rstrip("\n").split("\t")[1:]]
        for idx, value in enumerate(parts):
            if not value or ": " not in value:
                sample_attributes.setdefault(sample_accessions[idx], {})
                continue
            key, val = value.split(": ", 1)
            sample_attributes.setdefault(sample_accessions[idx], {})[key.strip()] = val.strip()

    return sample_titles, sample_accessions, source_names, sample_attributes


def read_clinical_records() -> tuple[set[str], dict[str, str]]:
    clinical_gsm_ids: set[str] = set()
    patient_to_clinical_gsm: dict[str, str] = {}

    with gzip.open(CLINICAL_PATH, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gsm_id = (row.get("GSE65194_ID") or "").strip()
            if not gsm_id:
                continue
            clinical_gsm_ids.add(gsm_id)

    with gzip.open(MATRIX_PATH, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    matrix_gsms: dict[str, str] = {}
    for line in lines:
        if line.startswith("!Sample_geo_accession\t"):
            accession_values = [value.strip('"') for value in line.rstrip("\n").split("\t")[1:]]
            for gsm_id in accession_values:
                matrix_gsms[gsm_id] = gsm_id
            break

    for line in lines:
        if line.startswith("!Sample_title\t"):
            title_values = [value.strip('"') for value in line.rstrip("\n").split("\t")[1:]]
            for gsm_id, title in zip(matrix_gsms.keys(), title_values):
                tumor_match = re.search(r"TUM\d+", title)
                if tumor_match and gsm_id in clinical_gsm_ids:
                    tumor_id = tumor_match.group(0)
                    patient_to_clinical_gsm.setdefault(tumor_id, gsm_id)

    return clinical_gsm_ids, patient_to_clinical_gsm


def main() -> None:
    sample_titles, sample_accessions, source_names, sample_attributes = read_matrix_annotations()
    clinical_gsm_ids, patient_to_clinical_gsm = read_clinical_records()

    rows: list[dict[str, str]] = []
    for idx, gsm_id in enumerate(sample_accessions):
        title = sample_titles[idx]
        source_name = source_names[idx]
        attributes = sample_attributes.get(gsm_id, {})

        tumor_match = re.search(r"TUM\d+", title)
        tumor_id = tumor_match.group(0) if tumor_match else ""
        sample_group = attributes.get("sample_group", "")

        replicate = ""
        if "_repA" in title:
            replicate = "repA"
        elif "_repB" in title:
            replicate = "repB"

        if tumor_id and sample_group not in {"Healthy", "CellLine"}:
            clinical_record_gsm_id = patient_to_clinical_gsm.get(tumor_id, "")
            patient_clinical_record_matched = "yes" if tumor_id in patient_to_clinical_gsm else "no"
            direct_gsm_clinical_match = "yes" if gsm_id in clinical_gsm_ids else "no"
        else:
            clinical_record_gsm_id = ""
            patient_clinical_record_matched = "no"
            direct_gsm_clinical_match = "no"

        row = {
            "gsm_id": gsm_id,
            "sample_title": title,
            "tumor_id": tumor_id,
            "replicate_label": replicate,
            "sample_group": sample_group,
            "array_batch": attributes.get("array_batch", ""),
            "hybridation_batch": attributes.get("hybridation batch", ""),
            "source_name": source_name,
            "cell_line": attributes.get("cell line", ""),
            "is_cell_line": "yes" if attributes.get("cell line") else "no",
            "platform": "GPL570",
            "clinical_record_gsm_id": clinical_record_gsm_id,
            "patient_clinical_record_matched": patient_clinical_record_matched,
            "direct_gsm_clinical_match": direct_gsm_clinical_match,
        }
        rows.append(row)

    fieldnames = [
        "gsm_id",
        "sample_title",
        "tumor_id",
        "replicate_label",
        "sample_group",
        "array_batch",
        "hybridation_batch",
        "source_name",
        "cell_line",
        "is_cell_line",
        "platform",
        "clinical_record_gsm_id",
        "patient_clinical_record_matched",
        "direct_gsm_clinical_match",
    ]

    with MANIFEST_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
