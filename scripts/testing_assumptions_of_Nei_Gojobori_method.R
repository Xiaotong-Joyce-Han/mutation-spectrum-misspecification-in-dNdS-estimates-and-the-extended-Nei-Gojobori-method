# input: data/MA_from_Long_et_al_and_Liu_and_Zhang_6type_counts.xlsm
# outfolder: results/ (can be modified in script)
# output: 
# - Nei_Gojobori_six_mutation_types.csv
# - Nei_Gojobori_two_transitions_equal_rates.csv
# - Nei_Gojobori_four_transversions_equal_rates.csv




# Global parameters
alpha = 0.05


if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Package 'readxl' is required to read the .xlsm input file. ",
       "Install it with install.packages('readxl').")
}


command_line <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_line, value = TRUE)
if (length(script_argument) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_argument), mustWork = FALSE)
  project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
} else {
  project_dir <- normalizePath(getwd(), mustWork = FALSE)
}


input_file <- file.path(project_dir, "data/MA_from_Long_et_al_and_Liu_and_Zhang_6type_counts.xlsm")
results_dir <- file.path(project_dir, "results")


if (!file.exists(input_file)) {
  stop("Input file was not found: ", input_file)
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


# The first row is a title and the second row contains the actual column names.
data <- as.data.frame(readxl::read_excel(
  input_file, sheet = 1, skip = 1, .name_repair = "minimal"
))
names(data) <- trimws(names(data))


required_columns <- c(
  "AT to CG", "AT to GC", "AT to TA",
  "CG to AT", "CG to TA", "CG to GC", "# AT sites", "#CG sites"
)
if (!all(required_columns %in% names(data))) {
  stop("The input table is missing required columns: ",
       paste(setdiff(required_columns, names(data)), collapse = ", "))
}


# Remove empty rows and the final note
has_species <- !is.na(data$Species) & trimws(as.character(data$Species)) != ""
has_numbers <- complete.cases(data[required_columns])
data <- data[has_species & has_numbers, , drop = FALSE]


numeric_columns <- required_columns
data[numeric_columns] <- lapply(data[numeric_columns], as.numeric)
if (anyNA(as.matrix(data[numeric_columns]))) {
  stop("Mutation-count or site-count columns contain missing/non-numeric values.")
}


# Run one goodness-of-fit test and return one row for the output table  
# 
# Note: usual chi-square approximation may be unreliable when
# at least one expected count is below five.
run_chisq_test <- function(species, observed, probabilities, labels) {
  
  test <- suppressWarnings(stats::chisq.test(x = observed, p = probabilities))
  
  result <- data.frame(
    Species = species,
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value,
    conclusion_at_alpha_0.05 = ifelse(test$p.value < alpha, "reject", "do_not_reject"),
    expected_count_below_5 = any(test$expected < 5),
    check.names = FALSE
  )
  
  for (index in seq_along(labels)) {
    result[[paste0("observed_", labels[index])]] <- observed[index]
    result[[paste0("expected_", labels[index])]] <- test$expected[index]
  }
  
  result
}


number_rows <- nrow(data)

# 1. Test all six mutation types together

six_type_results <- NULL

for (row_index in seq_len(number_rows)) {
  
  observed <- unlist(data[row_index, c(
    "AT to CG", "AT to GC", "AT to TA",
    "CG to AT", "CG to TA", "CG to GC"
  )])
  
  site_counts <- unlist(data[row_index, c("# AT sites", "#CG sites")])
  probabilities <- c(rep(site_counts[1], 3), rep(site_counts[2], 3))
  probabilities <- probabilities / sum(probabilities)

  row_result <- run_chisq_test(
    data$Species[row_index], observed, probabilities,
    c("AT_to_CG", "AT_to_GC", "AT_to_TA", "CG_to_AT", "CG_to_TA", "CG_to_GC")
  )
  
  six_type_results <- rbind(six_type_results, row_result)
}


# The corrections are applied across the 19 row-wise six-type tests.
six_type_results$p_value_bonferroni <- stats::p.adjust(six_type_results$p_value, "bonferroni")
six_type_results$conclusion_bonferroni_at_alpha_0.05 <- ifelse(six_type_results$p_value_bonferroni < alpha, "reject", "do_not_reject")

six_type_results$p_value_BH <- stats::p.adjust(six_type_results$p_value, "BH")
six_type_results$conclusion_BH_at_alpha_0.05 <- ifelse(six_type_results$p_value_BH < alpha, "reject", "do_not_reject")

utils::write.csv(
  six_type_results,
  file.path(results_dir, "Nei_Gojobori_six_mutation_types.csv"),
  row.names = FALSE
)

message("Wrote five row-wise chi-square result tables to: ", results_dir)


# 2. Test equal per-site rates for the two transitions: A->G (AT to GC) and
# C->T (CG to TA). 
transition_equal_rate_results <- NULL

for (row_index in seq_len(number_rows)) {
  observed <- c(
    data[row_index, "AT to GC"],
    data[row_index, "CG to TA"]
  )
  
  site_counts <- unlist(data[row_index, c("# AT sites", "#CG sites")])
  probabilities <- site_counts / sum(site_counts)

  row_result <- run_chisq_test(
    data$Species[row_index], observed, probabilities,
    c("A_to_G_transition", "C_to_T_transition")
  )
  transition_equal_rate_results <- rbind(transition_equal_rate_results, row_result)
}

transition_equal_rate_results$p_value_bonferroni <- stats::p.adjust(transition_equal_rate_results$p_value, "bonferroni")
transition_equal_rate_results$conclusion_bonferroni_at_alpha_0.05 <- ifelse(transition_equal_rate_results$p_value_bonferroni < alpha, "reject", "do_not_reject")

transition_equal_rate_results$p_value_BH <- stats::p.adjust(transition_equal_rate_results$p_value, "BH")
transition_equal_rate_results$conclusion_BH_at_alpha_0.05 <- ifelse(transition_equal_rate_results$p_value_BH < alpha, "reject", "do_not_reject")

utils::write.csv(
  transition_equal_rate_results,
  file.path(results_dir, "Nei_Gojobori_two_transitions_equal_rates.csv"),
  row.names = FALSE
)


# 3. Test equal per-site rates for the four transversions. 
transversion_equal_rate_results <- NULL

for (row_index in seq_len(number_rows)) {
  
  observed <- unlist(data[row_index, c(
    "AT to CG", "AT to TA", "CG to AT", "CG to GC"
  )])
  
  site_counts <- unlist(data[row_index, c("# AT sites", "#CG sites")])
  probabilities <- c(rep(site_counts[1], 2), rep(site_counts[2], 2))
  probabilities <- probabilities / sum(probabilities)

  row_result <- run_chisq_test(
    data$Species[row_index], observed, probabilities,
    c("AT_to_CG_transversion", "AT_to_TA_transversion",
      "CG_to_AT_transversion", "CG_to_GC_transversion")
  )
  
  transversion_equal_rate_results <- rbind(transversion_equal_rate_results, row_result)
}

transversion_equal_rate_results$p_value_bonferroni <- stats::p.adjust(transversion_equal_rate_results$p_value, "bonferroni")
transversion_equal_rate_results$conclusion_bonferroni_at_alpha_0.05 <- ifelse(transversion_equal_rate_results$p_value_bonferroni < alpha, "reject", "do_not_reject")

transversion_equal_rate_results$p_value_BH <- stats::p.adjust(transversion_equal_rate_results$p_value, "BH")
transversion_equal_rate_results$conclusion_BH_at_alpha_0.05 <- ifelse(transversion_equal_rate_results$p_value_BH < alpha, "reject", "do_not_reject")

utils::write.csv(
  transversion_equal_rate_results,
  file.path(results_dir, "Nei_Gojobori_four_transversions_equal_rates.csv"),
  row.names = FALSE
)

message("Wrote two additional row-wise chi-square result tables to: ", results_dir)

