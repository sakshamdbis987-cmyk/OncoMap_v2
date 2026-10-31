project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
local_library <- normalizePath(file.path(project_root, ".Rlib"), winslash = "/", mustWork = TRUE)
.libPaths(c(local_library, .libPaths()))

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("limma is not available in .Rlib or the configured R libraries.", call. = FALSE)
}

limma_version <- as.character(utils::packageVersion("limma"))
r_version <- paste(R.version$major, R.version$minor, sep = ".")

matrix_file <- file.path(project_root, "data", "raw", "GSE65194_series_matrix.txt.gz")
rep_a_file <- file.path(project_root, "reports", "GSE65194_GEO2R_repA_group_assignments.tsv")
rep_b_file <- file.path(project_root, "reports", "GSE65194_GEO2R_repB_group_assignments.tsv")
welch_file <- file.path(project_root, "reports", "GSE65194_tumor_vs_healthy_welch.tsv")
context_files <- c(
  file.path(project_root, "reports", "GSE65194_differential_expression_report.md"),
  file.path(project_root, "reports", "GSE65194_patient_level_aggregation_report.md")
)

output_files <- c(
  file.path(project_root, "data", "processed", "GSE65194_limma_repA_full_table.tsv"),
  file.path(project_root, "data", "processed", "GSE65194_limma_repB_full_table.tsv"),
  file.path(project_root, "reports", "GSE65194_limma_repA_vs_repB.tsv"),
  file.path(project_root, "reports", "GSE65194_limma_vs_welch_report.md")
)

if (any(file.exists(output_files))) {
  stop(
    "Refusing to overwrite existing output(s): ",
    paste(output_files[file.exists(output_files)], collapse = ", "),
    call. = FALSE
  )
}

input_files <- c(matrix_file, rep_a_file, rep_b_file, welch_file, context_files)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "), call. = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = character()
  )
}

read_assignment <- function(path, label) {
  assignment <- read_tsv(path)
  required_columns <- c(
    "gsm_id", "patient_id", "sample_title", "original_subtype",
    "geo2r_group", "selection_reason"
  )
  if (!all(required_columns %in% names(assignment))) {
    stop(label, " assignment sheet is missing required columns.", call. = FALSE)
  }
  if (nrow(assignment) != 141L || anyDuplicated(assignment$gsm_id)) {
    stop(label, " must contain exactly 141 unique GSMs.", call. = FALSE)
  }
  if (sum(assignment$geo2r_group == "Tumour") != 130L ||
      sum(assignment$geo2r_group == "Healthy") != 11L ||
      any(!assignment$geo2r_group %in% c("Tumour", "Healthy"))) {
    stop(label, " must contain exactly 130 Tumour and 11 Healthy assignments.", call. = FALSE)
  }
  tumour <- assignment[assignment$geo2r_group == "Tumour", , drop = FALSE]
  healthy <- assignment[assignment$geo2r_group == "Healthy", , drop = FALSE]
  if (anyNA(tumour$patient_id) || any(!nzchar(tumour$patient_id)) ||
      length(unique(tumour$patient_id)) != 130L) {
    stop(label, " must represent 130 unique tumour patients.", call. = FALSE)
  }
  healthy_ids <- sub("^Healthy,\\s*(H[0-9]+).*$", "\\1", healthy$sample_title)
  if (any(healthy_ids == healthy$sample_title) ||
      length(unique(healthy_ids)) != 11L) {
    stop(label, " must represent 11 uniquely identified healthy controls.", call. = FALSE)
  }
  assignment
}

cat("limma version:", limma_version, "\n")
cat("R version:", r_version, "\n")
cat("Validating inputs and reading the GEO series matrix...\n")

rep_a <- read_assignment(rep_a_file, "repA")
rep_b <- read_assignment(rep_b_file, "repB")

readLinesGz <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection))
  readLines(connection, warn = FALSE)
}

matrix_lines <- readLinesGz(matrix_file)
begin_marker <- which(trimws(matrix_lines) == "!series_matrix_table_begin")
end_marker <- which(trimws(matrix_lines) == "!series_matrix_table_end")
if (length(begin_marker) != 1L || length(end_marker) != 1L ||
    begin_marker >= end_marker) {
  stop("Could not uniquely locate the series-matrix table markers.", call. = FALSE)
}

table_text <- paste(matrix_lines[seq.int(begin_marker + 1L, end_marker - 1L)], collapse = "\n")
matrix_table <- utils::read.delim(
  text = table_text,
  header = TRUE,
  sep = "\t",
  quote = "\"",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = character()
)
rm(matrix_lines, table_text)

if (!identical(names(matrix_table)[1L], "ID_REF") ||
    nrow(matrix_table) != 54673L ||
    ncol(matrix_table) != 179L) {
  stop(
    "Expected an ID_REF column plus 178 GSM columns and exactly 54,673 probes; found ",
    nrow(matrix_table), " rows and ", ncol(matrix_table) - 1L, " source columns.",
    call. = FALSE
  )
}
gsm_columns <- names(matrix_table)[-1L]
if (anyDuplicated(gsm_columns) || any(!grepl("^GSM[0-9]+$", gsm_columns))) {
  stop("The source matrix must have 178 unique GSM column names.", call. = FALSE)
}
probe_ids <- as.character(matrix_table[[1L]])
if (anyNA(probe_ids) || any(!nzchar(probe_ids)) || anyDuplicated(probe_ids)) {
  stop("Probe IDs must be present and unique.", call. = FALSE)
}

numeric_columns <- lapply(matrix_table[-1L], function(column) {
  suppressWarnings(as.numeric(column))
})
invalid_values <- vapply(numeric_columns, function(column) any(!is.finite(column)), logical(1))
if (any(invalid_values)) {
  stop(
    "Non-numeric, missing, or non-finite values found in source column(s): ",
    paste(gsm_columns[invalid_values], collapse = ", "),
    call. = FALSE
  )
}
expression <- do.call(cbind, numeric_columns)
dimnames(expression) <- list(probe_ids, gsm_columns)
rm(matrix_table, numeric_columns)

if (!all(rep_a$gsm_id %in% gsm_columns) || !all(rep_b$gsm_id %in% gsm_columns)) {
  stop("An assignment GSM is absent from the source matrix.", call. = FALSE)
}

tumour_a <- rep_a[rep_a$geo2r_group == "Tumour", , drop = FALSE]
tumour_b <- rep_b[rep_b$geo2r_group == "Tumour", , drop = FALSE]
healthy_a <- rep_a[rep_a$geo2r_group == "Healthy", , drop = FALSE]
healthy_b <- rep_b[rep_b$geo2r_group == "Healthy", , drop = FALSE]
if (!setequal(tumour_a$patient_id, tumour_b$patient_id) ||
    !setequal(healthy_a$gsm_id, healthy_b$gsm_id)) {
  stop("repA and repB do not represent the same tumour patients and healthy controls.", call. = FALSE)
}

tumour_a <- tumour_a[match(sort(tumour_a$patient_id), tumour_a$patient_id), , drop = FALSE]
tumour_b <- tumour_b[match(tumour_a$patient_id, tumour_b$patient_id), , drop = FALSE]
if (!identical(tumour_a$original_subtype, tumour_b$original_subtype)) {
  stop("Tumour patient metadata differs between repA and repB.", call. = FALSE)
}

changed_patients <- tumour_a$patient_id[tumour_a$gsm_id != tumour_b$gsm_id]
same_patients <- tumour_a$patient_id[tumour_a$gsm_id == tumour_b$gsm_id]
if (length(changed_patients) != 23L ||
    any(tumour_a$selection_reason[match(changed_patients, tumour_a$patient_id)] != "repA") ||
    any(tumour_b$selection_reason[match(changed_patients, tumour_b$patient_id)] != "repB") ||
    any(tumour_a$selection_reason[match(same_patients, tumour_a$patient_id)] != "singleton") ||
    any(tumour_b$selection_reason[match(same_patients, tumour_b$patient_id)] != "singleton")) {
  stop("Assignment differences are not limited to the expected paired-patient arrays.", call. = FALSE)
}
normalized_a_titles <- sub("_repA$", "", tumour_a$sample_title)
normalized_b_titles <- sub("_repB$", "", tumour_b$sample_title)
if (!identical(normalized_a_titles, normalized_b_titles) ||
    !setequal(healthy_a$gsm_id, healthy_b$gsm_id) ||
    !identical(
      healthy_a$sample_title[match(sort(healthy_a$gsm_id), healthy_a$gsm_id)],
      healthy_b$sample_title[match(sort(healthy_a$gsm_id), healthy_b$gsm_id)]
    )) {
  stop("repA and repB differ outside their selected paired tumour arrays.", call. = FALSE)
}

if (!all(c("probe_id", "gene_symbol", "mean_log2_difference_tumour_minus_healthy",
          "p_value_two_sided", "bh_adjusted_p_value") %in% names(read_tsv(welch_file)))) {
  stop("The Welch results file is missing required comparison columns.", call. = FALSE)
}
welch <- read_tsv(welch_file)
if (nrow(welch) != 38L || anyDuplicated(welch$probe_id) ||
    !all(welch$probe_id %in% probe_ids)) {
  stop("Welch results must contain 38 unique probe IDs present in the source matrix.", call. = FALSE)
}
target_probes <- welch$probe_id

run_limma <- function(assignment, label) {
  selected <- expression[, assignment$gsm_id, drop = FALSE]
  group <- factor(assignment$geo2r_group, levels = c("Tumour", "Healthy"))
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  contrast <- limma::makeContrasts(Tumour - Healthy, levels = design)
  fit <- limma::lmFit(selected, design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(fit)
  result <- limma::topTable(
    fit,
    coef = 1L,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )
  if (nrow(result) != 54673L || !identical(rownames(result), rownames(expression))) {
    stop(label, " limma output is not the full unfiltered probe set.", call. = FALSE)
  }
  data.frame(probe_id = rownames(result), result, check.names = FALSE, row.names = NULL)
}

cat("Running repA and repB patient-balanced limma fits...\n")
result_a <- run_limma(rep_a, "repA")
result_b <- run_limma(rep_b, "repB")

write_tsv <- function(data, path) {
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

write_tsv(result_a, output_files[1L])
write_tsv(result_b, output_files[2L])

direction <- function(values) {
  ifelse(values > 0, "Tumour higher", ifelse(values < 0, "Tumour lower", "No change"))
}
significant <- function(values) values < 0.05

target_a <- result_a[match(target_probes, result_a$probe_id), , drop = FALSE]
target_b <- result_b[match(target_probes, result_b$probe_id), , drop = FALSE]
replicate_comparison <- data.frame(
  probe_id = target_probes,
  gene_symbol = welch$gene_symbol,
  logFC_repA = target_a$logFC,
  direction_repA = direction(target_a$logFC),
  logFC_repB = target_b$logFC,
  direction_repB = direction(target_b$logFC),
  logFC_difference_repA_minus_repB = target_a$logFC - target_b$logFC,
  p_value_repA = target_a$P.Value,
  p_value_repB = target_b$P.Value,
  adj_p_value_repA = target_a$adj.P.Val,
  adj_p_value_repB = target_b$adj.P.Val,
  significant_repA_adj_p_lt_0_05 = significant(target_a$adj.P.Val),
  significant_repB_adj_p_lt_0_05 = significant(target_b$adj.P.Val),
  stringsAsFactors = FALSE
)
write_tsv(replicate_comparison, output_files[3L])

welch_comparison <- data.frame(
  probe_id = target_probes,
  gene_symbol = welch$gene_symbol,
  welch_log2_difference = welch$mean_log2_difference_tumour_minus_healthy,
  limma_logFC_repA = target_a$logFC,
  limma_logFC_repB = target_b$logFC,
  welch_p_value = welch$p_value_two_sided,
  limma_p_value_repA = target_a$P.Value,
  limma_p_value_repB = target_b$P.Value,
  welch_adj_p_value = welch$bh_adjusted_p_value,
  limma_adj_p_value_repA = target_a$adj.P.Val,
  limma_adj_p_value_repB = target_b$adj.P.Val,
  direction_agrees_repA = sign(welch$mean_log2_difference_tumour_minus_healthy) == sign(target_a$logFC),
  direction_agrees_repB = sign(welch$mean_log2_difference_tumour_minus_healthy) == sign(target_b$logFC),
  significant_welch = significant(welch$bh_adjusted_p_value),
  significant_limma_repA = significant(target_a$adj.P.Val),
  significant_limma_repB = significant(target_b$adj.P.Val),
  stringsAsFactors = FALSE
)

count_significance <- function(limma_significant) {
  c(
    both = sum(significant(welch$bh_adjusted_p_value) & limma_significant),
    welch_only = sum(significant(welch$bh_adjusted_p_value) & !limma_significant),
    limma_only = sum(!significant(welch$bh_adjusted_p_value) & limma_significant),
    neither = sum(!significant(welch$bh_adjusted_p_value) & !limma_significant)
  )
}

format_value <- function(value) {
  ifelse(is.na(value), "NA", formatC(value, digits = 4L, format = "g"))
}
table_rows <- vapply(seq_len(nrow(welch_comparison)), function(i) {
  row <- welch_comparison[i, ]
  paste0(
    "| ", row$probe_id, " | ", row$gene_symbol,
    " | ", format_value(row$welch_log2_difference),
    " | ", format_value(row$limma_logFC_repA),
    " | ", format_value(row$limma_logFC_repB),
    " | ", format_value(row$welch_adj_p_value),
    " | ", format_value(row$limma_adj_p_value_repA),
    " | ", format_value(row$limma_adj_p_value_repB),
    " | ", ifelse(row$direction_agrees_repA, "Yes", "No"),
    " | ", ifelse(row$direction_agrees_repB, "Yes", "No"),
    " | ", ifelse(row$significant_welch, "Yes", "No"),
    " | ", ifelse(row$significant_limma_repA, "Yes", "No"),
    " | ", ifelse(row$significant_limma_repB, "Yes", "No"), " |"
  )
}, character(1))

sig_a <- count_significance(welch_comparison$significant_limma_repA)
sig_b <- count_significance(welch_comparison$significant_limma_repB)
report_lines <- c(
  "# GSE65194 local limma sensitivity analysis",
  "",
  "## Validation and reproducibility",
  "",
  paste0("- R version: ", r_version, " (", R.version$platform, ")."),
  paste0("- limma version: ", limma_version, "."),
  "- Input GEO matrix: 54,673 probes and 178 source GSM columns; the expression table was read between the GEO series-matrix table markers.",
  "- Assignment validation: each sheet contains 141 unique GSMs (130 Tumour and 11 Healthy), representing 130 unique tumour patients and 11 uniquely identified healthy controls.",
  paste0("- repA and repB differ only in the selected arrays for ", length(changed_patients), " paired tumour patients; 107 singleton tumour selections and all 11 healthy controls are unchanged."),
  "- Group sizes per fit: Tumour = 130; Healthy = 11.",
  "- Contrast: Tumour - Healthy. Positive logFC means higher expression in Tumour.",
  "- The matrix contains processed log2 expression values; no further log transformation, normalization, batch correction, sample removal, or probe filtering was applied.",
  "- repA and repB are alternative technical-replicate selections, not independent cohorts.",
  "",
  "## Methods",
  "",
  "Both full-probe analyses used `limma::lmFit()`, `limma::contrasts.fit()`, `limma::eBayes()`, and `limma::topTable()`, with BH adjustment across all 54,673 probes. The full results retain all probes.",
  "",
  "The comparison with the prior Welch analysis uses exact probe IDs. Welch BH-adjusted p-values were calculated across its 38 target probes, while limma BH-adjusted p-values were calculated across all 54,673 probes; significance differences should be interpreted with these different correction families in mind. Significance is adjusted p < 0.05.",
  "",
  "This is a local limma analysis of the downloaded GEO matrix using the planned sample assignments. It does not claim exact reproduction of the GEO2R website.",
  "",
  "## Limma versus Welch: probe-level comparison",
  "",
  paste0(
    "- Directional agreement with Welch: repA ",
    sum(welch_comparison$direction_agrees_repA), "/38; repB ",
    sum(welch_comparison$direction_agrees_repB), "/38."
  ),
  paste0(
    "- Adjusted-significance comparison (both / Welch only / limma only / neither): repA ",
    paste(sig_a, collapse = " / "), "; repB ",
    paste(sig_b, collapse = " / "), "."
  ),
  "",
  "| Probe ID | Gene | Welch tumour - healthy | limma repA logFC | limma repB logFC | Welch BH p | repA BH p | repB BH p | Direction agrees repA | Direction agrees repB | Welch significant | repA significant | repB significant |",
  "|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|---|",
  table_rows,
  ""
)
writeLines(report_lines, output_files[4L], useBytes = TRUE)

cat("Validation passed: 54,673 probes, 178 source GSMs, and 141 assignments per sheet.\n")
cat("Group sizes per run: Tumour = 130; Healthy = 11.\n")
cat("Changed paired-patient selections:", length(changed_patients), "; unchanged singleton selections:", length(same_patients), ".\n")
cat("Significant probes across full matrix (adjusted p < 0.05): repA =",
    sum(result_a$adj.P.Val < 0.05), "; repB =", sum(result_b$adj.P.Val < 0.05), "\n")
cat("Significant target probes (adjusted p < 0.05): repA =",
    sum(significant(target_a$adj.P.Val)), "; repB =", sum(significant(target_b$adj.P.Val)),
    "; Welch =", sum(significant(welch$bh_adjusted_p_value)), "\n")
cat("Welch directional agreement (of 38): repA =",
    sum(welch_comparison$direction_agrees_repA), "; repB =",
    sum(welch_comparison$direction_agrees_repB), "\n")
cat("Created:\n", paste(output_files, collapse = "\n"), "\n")
