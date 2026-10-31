project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
local_library <- normalizePath(file.path(project_root, ".Rlib"), winslash = "/", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("limma is not available in .Rlib or the configured R libraries.", call. = FALSE)
}

limma_version <- as.character(utils::packageVersion("limma"))
r_version <- paste(R.version$major, R.version$minor, sep = ".")

matrix_file <- file.path(project_root, "data", "processed", "GSE65194_patient_level_matrix.tsv.gz")
manifest_file <- file.path(project_root, "data", "processed", "GSE65194_patient_level_manifest.tsv")
welch_file <- file.path(project_root, "reports", "GSE65194_subtype_pairwise_welch.tsv")
compatibility_file <- file.path(project_root, "reports", "GSE65194_compatibility_audit.md")
output_files <- c(
  file.path(project_root, "reports", "GSE65194_limma_subtype_pairwise.tsv"),
  file.path(project_root, "reports", "GSE65194_limma_subtype_vs_welch_report.md")
)

existing_outputs <- output_files[file.exists(output_files)]
if (length(existing_outputs) > 0L) {
  stop(
    "Refusing to overwrite existing output(s): ",
    paste(existing_outputs, collapse = ", "),
    call. = FALSE
  )
}
inputs <- c(matrix_file, manifest_file, welch_file, compatibility_file)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs) > 0L) {
  stop("Missing required input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

subtypes <- c("TNBC", "Her2", "Luminal A", "Luminal B")
expected_counts <- c("TNBC" = 41L, "Her2" = 30L, "Luminal A" = 29L, "Luminal B" = 30L)
expected_pairs <- data.frame(
  first_subtype = c("TNBC", "TNBC", "TNBC", "Her2", "Her2", "Luminal A"),
  second_subtype = c("Her2", "Luminal A", "Luminal B", "Luminal A", "Luminal B", "Luminal B"),
  stringsAsFactors = FALSE
)

probe_to_gene <- c(
  "204531_s_at" = "BRCA1",
  "211851_x_at" = "BRCA1",
  "208368_s_at" = "BRCA2",
  "214727_at" = "BRCA2",
  "201746_at" = "TP53",
  "211300_s_at" = "TP53",
  "210930_s_at" = "ERBB2",
  "216836_s_at" = "ERBB2",
  "234354_x_at" = "ERBB2",
  "204369_at" = "PIK3CA",
  "231854_at" = "PIK3CA",
  "235980_at" = "PIK3CA",
  "205225_at" = "ESR1",
  "211233_x_at" = "ESR1",
  "211234_x_at" = "ESR1",
  "211235_s_at" = "ESR1",
  "211627_x_at" = "ESR1",
  "215551_at" = "ESR1",
  "215552_s_at" = "ESR1",
  "217163_at" = "ESR1",
  "217190_x_at" = "ESR1",
  "204053_x_at" = "PTEN",
  "204054_at" = "PTEN",
  "211711_s_at" = "PTEN",
  "217492_s_at" = "PTEN",
  "225363_at" = "PTEN",
  "228006_at" = "PTEN",
  "233254_x_at" = "PTEN",
  "242622_x_at" = "PTEN",
  "1565483_at" = "EGFR",
  "1565484_x_at" = "EGFR",
  "201983_s_at" = "EGFR",
  "201984_s_at" = "EGFR",
  "210984_x_at" = "EGFR",
  "211550_at" = "EGFR",
  "211551_at" = "EGFR",
  "211607_x_at" = "EGFR",
  "224999_at" = "EGFR"
)
probe_ids <- names(probe_to_gene)
probe_annotation <- probe_to_gene
probe_annotation["217492_s_at"] <- "PTENP1 /// PTEN"

compatibility_text <- paste(readLines(compatibility_file, warn = FALSE), collapse = "\n")
audit_probe_mapping <- list()
for (line in strsplit(compatibility_text, "\n", fixed = TRUE)[[1L]]) {
  if (grepl("^\\| (BRCA1|BRCA2|TP53|ERBB2|PIK3CA|ESR1|PTEN|EGFR) \\|", line)) {
    cells <- trimws(strsplit(line, "\\|", perl = TRUE)[[1L]])
    if (length(cells) >= 4L) {
      audit_probe_mapping[[cells[2L]]] <- trimws(strsplit(cells[3L], ",", fixed = TRUE)[[1L]])
    }
  }
}
expected_gene_order <- c("BRCA1", "BRCA2", "TP53", "ERBB2", "PIK3CA", "ESR1", "PTEN", "EGFR")
for (gene in expected_gene_order) {
  expected <- names(probe_to_gene)[probe_to_gene == gene]
  if (!identical(audit_probe_mapping[[gene]], expected)) {
    stop("Target probe mapping for ", gene, " does not match the compatibility audit.", call. = FALSE)
  }
}
if (!grepl("217492_s_at", compatibility_text, fixed = TRUE) ||
    !grepl("PTENP1///PTEN", compatibility_text, fixed = TRUE) ||
    length(probe_ids) != 38L || anyDuplicated(probe_ids)) {
  stop("The validated 38-probe target list or ambiguous PTEN annotation failed validation.", call. = FALSE)
}

read_gzip_tsv <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  utils::read.delim(
    connection,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = character()
  )
}
matrix_table <- read_gzip_tsv(matrix_file)
if (!identical(names(matrix_table)[1L], "ID_REF") ||
    nrow(matrix_table) != 54673L ||
    ncol(matrix_table) != 142L) {
  stop(
    "Expected patient-level matrix dimensions 54,673 probes × 141 samples; found ",
    nrow(matrix_table), " rows and ", ncol(matrix_table) - 1L, " samples.",
    call. = FALSE
  )
}
matrix_probe_ids <- as.character(matrix_table[[1L]])
if (anyNA(matrix_probe_ids) || any(!nzchar(matrix_probe_ids)) || anyDuplicated(matrix_probe_ids)) {
  stop("Matrix probe IDs must be present and unique.", call. = FALSE)
}
if (!all(probe_ids %in% matrix_probe_ids)) {
  stop(
    "Target probe IDs missing from the patient matrix: ",
    paste(setdiff(probe_ids, matrix_probe_ids), collapse = ", "),
    call. = FALSE
  )
}

manifest <- utils::read.delim(
  manifest_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = character()
)
required_manifest_columns <- c(
  "sample_id", "sample_type", "patient_id", "sample_group", "subtype"
)
if (!all(required_manifest_columns %in% names(manifest)) ||
    nrow(manifest) != 141L ||
    anyNA(manifest$sample_id) || any(!nzchar(manifest$sample_id)) ||
    anyDuplicated(manifest$sample_id) ||
    !identical(manifest$sample_id, names(matrix_table)[-1L])) {
  stop("Manifest dimensions, required fields, unique IDs, or matrix column order failed validation.", call. = FALSE)
}
if (sum(manifest$sample_type == "tumour") != 130L ||
    sum(manifest$sample_type == "healthy_control") != 11L ||
    any(!manifest$sample_type %in% c("tumour", "healthy_control"))) {
  stop("Expected 130 tumour patients and 11 healthy controls in the patient-level manifest.", call. = FALSE)
}
tumour_manifest <- manifest[manifest$sample_type == "tumour", , drop = FALSE]
if (anyNA(tumour_manifest$patient_id) || any(!nzchar(tumour_manifest$patient_id)) ||
    anyDuplicated(tumour_manifest$patient_id)) {
  stop("Tumour patient IDs must be complete and unique.", call. = FALSE)
}
observed_counts <- table(factor(tumour_manifest$subtype, levels = subtypes))
if (!identical(as.integer(observed_counts), unname(expected_counts)) ||
    any(tumour_manifest$sample_group != tumour_manifest$subtype)) {
  stop("Tumour subtype counts or subtype labels do not match the expected values.", call. = FALSE)
}
numeric_values <- lapply(matrix_table[-1L], function(column) {
  suppressWarnings(as.numeric(column))
})
invalid_columns <- vapply(numeric_values, function(column) any(!is.finite(column)), logical(1))
if (any(invalid_columns)) {
  stop(
    "Missing, non-numeric, or non-finite patient matrix values in sample(s): ",
    paste(names(invalid_columns)[invalid_columns], collapse = ", "),
    call. = FALSE
  )
}
expression <- do.call(cbind, numeric_values)
rownames(expression) <- matrix_probe_ids
colnames(expression) <- names(matrix_table)[-1L]
rm(matrix_table, numeric_values)

welch <- utils::read.delim(
  welch_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = character()
)
required_welch_columns <- c(
  "first_subtype", "second_subtype", "probe_id", "gene_symbol",
  "mean_log2_difference_first_minus_second", "p_value_two_sided",
  "bh_adjusted_p_value"
)
if (!all(required_welch_columns %in% names(welch)) ||
    nrow(welch) != 228L ||
    anyDuplicated(paste(welch$first_subtype, welch$second_subtype, welch$probe_id, sep = "\r"))) {
  stop("Welch results must contain 228 unique probe-by-contrast rows and required fields.", call. = FALSE)
}
expected_key <- unlist(lapply(seq_len(nrow(expected_pairs)), function(i) {
  paste(
    expected_pairs$first_subtype[i],
    expected_pairs$second_subtype[i],
    probe_ids,
    sep = "\r"
  )
}), use.names = FALSE)
welch_key <- paste(welch$first_subtype, welch$second_subtype, welch$probe_id, sep = "\r")
if (!setequal(welch_key, expected_key) ||
    any(welch$gene_symbol != unname(probe_to_gene[welch$probe_id]))) {
  stop("Welch results do not exactly cover the expected contrasts and 38 probe IDs.", call. = FALSE)
}

selected_samples <- tumour_manifest$sample_id
selected_expression <- expression[, selected_samples, drop = FALSE]
group <- factor(tumour_manifest$subtype, levels = subtypes)
design <- stats::model.matrix(~ 0 + group)
colnames(design) <- make.names(subtypes)
contrast_matrix <- limma::makeContrasts(
  TNBC - Her2,
  TNBC - Luminal.A,
  TNBC - Luminal.B,
  Her2 - Luminal.A,
  Her2 - Luminal.B,
  Luminal.A - Luminal.B,
  levels = design
)
contrast_names <- c(
  "TNBC - Her2",
  "TNBC - Luminal A",
  "TNBC - Luminal B",
  "Her2 - Luminal A",
  "Her2 - Luminal B",
  "Luminal A - Luminal B"
)
colnames(contrast_matrix) <- contrast_names
fit <- limma::lmFit(selected_expression, design)
fit <- limma::contrasts.fit(fit, contrast_matrix)
fit <- limma::eBayes(fit)

results_by_contrast <- lapply(seq_along(contrast_names), function(i) {
  top <- limma::topTable(
    fit,
    coef = i,
    number = Inf,
    adjust.method = "none",
    sort.by = "none"
  )
  if (nrow(top) != nrow(selected_expression) ||
      !identical(rownames(top), rownames(selected_expression))) {
    stop("limma did not return every probe for contrast ", contrast_names[i], ".", call. = FALSE)
  }
  selected <- top[probe_ids, , drop = FALSE]
  pair <- expected_pairs[i, ]
  data.frame(
    first_subtype = pair$first_subtype,
    second_subtype = pair$second_subtype,
    contrast = contrast_names[i],
    probe_id = probe_ids,
    gene_symbol = unname(probe_to_gene[probe_ids]),
    gene_annotation = unname(probe_annotation[probe_ids]),
    logFC = selected$logFC,
    AveExpr = selected$AveExpr,
    t = selected$t,
    p_value = selected$P.Value,
    stringsAsFactors = FALSE
  )
})
limma_results <- do.call(rbind, results_by_contrast)
if (nrow(limma_results) != 228L || any(!is.finite(limma_results$p_value))) {
  stop("Expected 228 finite raw p-values from limma.", call. = FALSE)
}
limma_results$bh_adjusted_p_value <- stats::p.adjust(limma_results$p_value, method = "BH")
limma_results$significant_bh_0_05 <- limma_results$bh_adjusted_p_value < 0.05
limma_results$direction <- ifelse(
  limma_results$logFC > 0,
  "First subtype higher",
  ifelse(limma_results$logFC < 0, "Second subtype higher", "No change")
)
names(limma_results)[names(limma_results) == "significant_bh_0_05"] <-
  "limma_significant_bh_0_05"

welch$contrast <- paste(welch$first_subtype, "-", welch$second_subtype)
names(welch)[names(welch) == "gene_symbol"] <- "gene_symbol_welch"
names(welch)[names(welch) == "p_value_two_sided"] <- "welch_p_value"
names(welch)[names(welch) == "bh_adjusted_p_value"] <- "welch_bh_adjusted_p_value"
names(welch)[names(welch) == "mean_log2_difference_first_minus_second"] <-
  "welch_log2_difference"
welch <- welch[, c(
  "first_subtype", "second_subtype", "contrast", "probe_id",
  "gene_symbol_welch", "welch_p_value", "welch_bh_adjusted_p_value",
  "welch_log2_difference"
)]
names(limma_results)[names(limma_results) == "gene_symbol"] <- "gene_symbol_limma"
names(limma_results)[names(limma_results) == "p_value"] <- "limma_p_value"
names(limma_results)[names(limma_results) == "bh_adjusted_p_value"] <-
  "limma_bh_adjusted_p_value"
comparison <- merge(
  limma_results,
  welch,
  by = c("first_subtype", "second_subtype", "contrast", "probe_id"),
  all = FALSE,
  sort = FALSE,
  suffixes = c("_limma", "_welch")
)
if (nrow(comparison) != 228L) {
  stop("Exact probe ID and contrast join did not produce all 228 comparisons.", call. = FALSE)
}
comparison$direction_agrees <- sign(comparison$logFC) ==
  sign(comparison$welch_log2_difference)
comparison$logFC_difference_limma_minus_welch <- comparison$logFC -
  comparison$welch_log2_difference
comparison$adjusted_p_difference_limma_minus_welch <-
  comparison$limma_bh_adjusted_p_value - comparison$welch_bh_adjusted_p_value
comparison$limma_significant <- comparison$limma_significant_bh_0_05
comparison$welch_significant <- comparison$welch_bh_adjusted_p_value < 0.05
comparison$significance_comparison <- ifelse(
  comparison$limma_significant & comparison$welch_significant,
  "Both",
  ifelse(
    comparison$limma_significant,
    "Limma only",
    ifelse(comparison$welch_significant, "Welch only", "Neither")
  )
)
comparison <- comparison[
  match(expected_key, paste(
    comparison$first_subtype,
    comparison$second_subtype,
    comparison$probe_id,
    sep = "\r"
  )),
  ,
  drop = FALSE
]

output_table <- comparison[, c(
  "first_subtype", "second_subtype", "contrast", "probe_id", "gene_symbol_limma",
  "gene_annotation", "logFC", "welch_log2_difference",
  "logFC_difference_limma_minus_welch", "limma_p_value", "welch_p_value",
  "limma_bh_adjusted_p_value", "welch_bh_adjusted_p_value",
  "adjusted_p_difference_limma_minus_welch", "direction", "direction_agrees",
  "limma_significant", "welch_significant", "significance_comparison"
)]
names(output_table)[names(output_table) == "gene_symbol_limma"] <- "gene_symbol"
names(output_table)[names(output_table) == "welch_log2_difference"] <-
  "welch_log2_difference_first_minus_second"
names(output_table)[names(output_table) == "direction_agrees"] <- "direction_agrees_with_welch"
names(output_table)[names(output_table) == "limma_significant"] <- "limma_significant_bh_0_05"
names(output_table)[names(output_table) == "welch_significant"] <- "welch_significant_bh_0_05"

utils::write.table(
  output_table,
  file = output_files[1L],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

summaries <- lapply(seq_len(nrow(expected_pairs)), function(i) {
  rows <- output_table[
    output_table$first_subtype == expected_pairs$first_subtype[i] &
      output_table$second_subtype == expected_pairs$second_subtype[i],
    ,
    drop = FALSE
  ]
  data.frame(
    contrast = contrast_names[i],
    limma_significant = sum(rows$limma_significant_bh_0_05),
    welch_significant = sum(rows$welch_significant_bh_0_05),
    significant_both = sum(rows$significance_comparison == "Both"),
    limma_only = sum(rows$significance_comparison == "Limma only"),
    welch_only = sum(rows$significance_comparison == "Welch only"),
    neither = sum(rows$significance_comparison == "Neither"),
    directional_agreement = sum(rows$direction_agrees_with_welch),
    stringsAsFactors = FALSE
  )
})
contrast_summary <- do.call(rbind, summaries)
total_summary <- data.frame(
  limma_significant = sum(output_table$limma_significant_bh_0_05),
  welch_significant = sum(output_table$welch_significant_bh_0_05),
  significant_both = sum(output_table$significance_comparison == "Both"),
  limma_only = sum(output_table$significance_comparison == "Limma only"),
  welch_only = sum(output_table$significance_comparison == "Welch only"),
  neither = sum(output_table$significance_comparison == "Neither"),
  directional_agreement = sum(output_table$direction_agrees_with_welch)
)
format_summary_row <- function(row) {
  paste0(
    "| ", row$contrast, " | ", row$limma_significant, " | ", row$welch_significant,
    " | ", row$significant_both, " | ", row$limma_only, " | ", row$welch_only,
    " | ", row$neither, " | ", row$directional_agreement, "/38 |"
  )
}
summary_lines <- vapply(seq_len(nrow(contrast_summary)), function(i) {
  format_summary_row(contrast_summary[i, ])
}, character(1))
number <- function(x) formatC(x, digits = 5L, format = "g")
probe_lines <- vapply(seq_len(nrow(output_table)), function(i) {
  row <- output_table[i, ]
  paste0(
    "| ", row$contrast, " | ", row$probe_id, " | ", row$gene_annotation,
    " | ", number(row$logFC), " | ", number(row$welch_log2_difference_first_minus_second),
    " | ", number(row$logFC_difference_limma_minus_welch),
    " | ", number(row$limma_p_value), " | ", number(row$welch_p_value),
    " | ", number(row$limma_bh_adjusted_p_value),
    " | ", number(row$welch_bh_adjusted_p_value),
    " | ", ifelse(row$direction_agrees_with_welch, "Yes", "No"),
    " | ", row$significance_comparison, " |"
  )
}, character(1))
report <- c(
  "# GSE65194 patient-level limma subtype analysis",
  "",
  "## Validation and reproducibility",
  "",
  paste0("- R version: ", r_version, " (", R.version$platform, ")."),
  paste0("- limma version: ", limma_version, "."),
  "- Expression input: existing patient-level averaged matrix; no re-averaging, log transformation, normalization, or batch correction was performed.",
  "- The patient-level matrix contains 54,673 probes × 141 observations. The analysis selected the 130 tumour patients only; 11 healthy controls were excluded. The 14 cell-line arrays are not in the patient-level matrix.",
  paste0(
    "- Tumour subtype counts: TNBC ", expected_counts["TNBC"],
    ", Her2 ", expected_counts["Her2"],
    ", Luminal A ", expected_counts["Luminal A"],
    ", Luminal B ", expected_counts["Luminal B"], "."
  ),
  "- All 38 target probes were validated against the compatibility audit and were present in the matrix. Probe measurements remain separate; no gene-level averaging was performed.",
  "- Design group order: TNBC, Her2, Luminal A, Luminal B.",
  "- Positive logFC means higher expression in the first-named subtype.",
  "",
  "## Model and multiple testing",
  "",
  "The model used `limma::lmFit()`, `limma::makeContrasts()`, `limma::contrasts.fit()`, and `limma::eBayes()`. Six pairwise contrasts were extracted for each of the 38 target probes. Benjamini-Hochberg adjustment was applied once across all 228 limma raw p-values; the Welch adjusted p-values use the existing complete 228-test family.",
  "",
  "## Comparison with Welch",
  "",
  paste0(
    "- Across all 228 probe-by-contrast results, limma has ",
    total_summary$limma_significant, " significant probes and Welch has ",
    total_summary$welch_significant, "."
  ),
  paste0(
    "- Significant in both methods: ", total_summary$significant_both,
    "; limma only: ", total_summary$limma_only,
    "; Welch only: ", total_summary$welch_only,
    "; neither: ", total_summary$neither, "."
  ),
  paste0(
    "- Directional agreement: ", total_summary$directional_agreement,
    "/228 probe-by-contrast results."
  ),
  "",
  "| Contrast | Limma significant | Welch significant | Both | Limma only | Welch only | Neither | Direction agrees |",
  "|---|---:|---:|---:|---:|---:|---:|---:|",
  summary_lines,
  "",
  "All comparisons below are probe-level and joined using exact probe IDs and contrast names. LogFC difference is limma minus Welch. Adjusted-p differences compare the same 228-test Welch family against limma BH correction across all 228 results.",
  "",
  "| Contrast | Probe ID | Annotation | limma logFC | Welch difference | logFC difference | limma raw p | Welch raw p | limma BH p | Welch BH p | Direction agrees | Significance |",
  "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|",
  probe_lines,
  ""
)
writeLines(report, output_files[2L], useBytes = TRUE)

cat("R version:", r_version, "; limma version:", limma_version, "\n")
cat("Validation passed: patient-level matrix 54,673 x 141; 130 tumour patients selected.\n")
cat("Subtype counts:", paste(names(expected_counts), expected_counts, sep = "=", collapse = ", "), "\n")
cat("Target probes validated:", length(probe_ids), "; limma/Welch joined results:", nrow(output_table), "\n")
cat(
  "Family-wide BH significant (of 228): limma =", total_summary$limma_significant,
  "; Welch =", total_summary$welch_significant,
  "; both =", total_summary$significant_both,
  "; limma only =", total_summary$limma_only,
  "; Welch only =", total_summary$welch_only,
  "; neither =", total_summary$neither, "\n"
)
cat("Directional agreement:", total_summary$directional_agreement, "/228\n")
cat("Per-contrast results:\n")
for (i in seq_len(nrow(contrast_summary))) {
  row <- contrast_summary[i, ]
  cat(
    row$contrast, ": limma=", row$limma_significant,
    ", Welch=", row$welch_significant,
    ", both=", row$significant_both,
    ", limma-only=", row$limma_only,
    ", Welch-only=", row$welch_only,
    ", neither=", row$neither,
    ", direction agreement=", row$directional_agreement, "/38\n",
    sep = ""
  )
}
cat("Created:\n", paste(output_files, collapse = "\n"), "\n")
