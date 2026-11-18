#Overall balance:
#Stuart-Maxwell test
#implementation: DescTools::StuartMaxwellTest
#null hypothesis: marginal hemogenity
#https://search.r-project.org/CRAN/refmans/DescTools/html/StuartMaxwellTest.html

#Detailed balance
#Bowker's test
#implement: stats::mcnemar.test
#null hypothesis: symmetry of the contingency table
#https://metricgate.com/docs/bowker-symmetry-test/
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/mcnemar.test

#input: data/MA_from_Long_et_al_and_Liu_and_Zhang_12type_counts.csv
#output dir: test_equilibirum_assumption_of_PAML (can be specified)
#output:
# - Goldman_Yang_AT_CG_equilibrium.csv
# - Goldman_Yang_Stuart_Maxwell_joint_equilibrium.csv
# - Goldman_Yang_detailed_balance_Bowker_joint_test.csv

library(DescTools)

alpha <- 0.05
bases <- c("A", "T", "C", "G")

# Resolve paths from the script location.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
project_dir <- if (length(script_arg) == 1) {
  dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
} else {
  getwd()
}
output_dir <- file.path(project_dir, "test_equilibirum_assumption_of_PAML")
input_file <- file.path("data/MA_from_Long_et_al_and_Liu_and_Zhang_12type_counts.csv")
data <- read.csv(input_file, check.names = FALSE)
data <- data[data$mutation_type != "total", , drop = FALSE]

counts <- as.matrix(data[-1])
rownames(counts) <- data$mutation_type

# Overall balance: Stuart-Maxwell testing
joint_results <- NULL
# i：species/condition index
for (i in seq_len(ncol(counts))) {
  mutations <- matrix(0, 4, 4, dimnames = list(bases, bases))

  for (from in bases) {
    for (to in bases[bases != from]) {
      mutations[from, to] <- counts[paste(from, "to", to, sep = "_"), i]
    }
  }
  
  print(mutations)
  
  test <- DescTools::StuartMaxwellTest(mutations)
  
  joint_results <- rbind(joint_results, data.frame(
    Species = colnames(counts)[i],
    statistic = as.numeric(test$statistic),
    df = as.numeric(test$parameter),
    p_value = test$p.value
  ))
}

joint_results$conclusion_at_alpha_0.05 <-
  ifelse(is.na(joint_results$p_value), "not_testable",
         ifelse(joint_results$p_value < alpha, "reject", "do_not_reject"))

for (method in c("bonferroni", "BH")) {
  adjusted <- p.adjust(joint_results$p_value, method = method, n = ncol(counts))
  joint_results[[paste0("p_value_", method)]] <- adjusted
  joint_results[[paste0("conclusion_", method, "_at_alpha_0.05")]] <-
    ifelse(is.na(adjusted), "not_testable",
           ifelse(adjusted < alpha, "reject", "do_not_reject"))
}
write.csv(joint_results, file.path(output_dir, "Goldman_Yang_Stuart_Maxwell_joint_equilibrium.csv"), row.names = FALSE)


# Detailed balance: six unique reciprocal pairs.
pairs <- combn(bases, 2)

# Bowker symmetry: mcnemar.test generalizes to the 4 x 4 table.
detailed_joint <- NULL

# different species/conditions
for (i in seq_len(ncol(counts))) {
  
  mutations <- matrix(0, 4, 4, dimnames = list(bases, bases))
  for (from in bases) {
    for (to in bases[bases != from]) {
      mutations[from, to] <- counts[paste(from, "to", to, sep = "_"), i]
    }
  }
  
  print(mutations)
  
  # extract the 6 total numbers as a int list
  pair_totals <- (mutations + t(mutations))[upper.tri(mutations)]
  
  low_counts <- sum(pair_totals / 2 < 5)
  df <- sum(pair_totals > 0)
  
  if (df == 6) {
    test <- stats::mcnemar.test(mutations, correct = FALSE)
    statistic <- as.numeric(test$statistic)
  } else {
    # Full-table mcnemar.test returns NaN for zero-total pairs.
    statistic <- 0
    for (j in seq_len(ncol(pairs))) {
      # 2*2 contingency block
      block <- mutations[pairs[, j], pairs[, j]]
      if (sum(block) > 0) {
        # If the block has 0 sum, that means the pair has 0 total.
        # And in that case, a NaN will be returned for the chi-square statistic.
        # It's excluded.
        statistic <- statistic + as.numeric(stats::mcnemar.test(block, correct = FALSE)$statistic)
      }
    }
  }
  if (df == 0) statistic <- NA_real_
  detailed_joint <- rbind(detailed_joint, data.frame(
    Species = colnames(counts)[i],
    statistic = statistic,
    df = df,
    p_value = pchisq(statistic, df, lower.tail = FALSE),
    n_pairs_expected_count_below_5 = low_counts
  ))
}

detailed_joint$conclusion_at_alpha_0.05 <-
  ifelse(is.na(detailed_joint$p_value), "not_testable",
         ifelse(detailed_joint$p_value < alpha, "reject", "do_not_reject"))

for (method in c("bonferroni", "BH")) {
  adjusted <- p.adjust(detailed_joint$p_value, method = method, n = ncol(counts))
  detailed_joint[[paste0("p_value_", method)]] <- adjusted
  detailed_joint[[paste0("conclusion_", method, "_at_alpha_0.05")]] <-
    ifelse(is.na(adjusted), "not_testable",
           ifelse(adjusted < alpha, "reject", "do_not_reject"))
}

write.csv(detailed_joint, file.path(output_dir, "Goldman_Yang_detailed_balance_Bowker_joint_test.csv"), row.names = FALSE)
