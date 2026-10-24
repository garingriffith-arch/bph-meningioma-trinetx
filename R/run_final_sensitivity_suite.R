# ============================================================
# run_final_sensitivity_suite_v2.R
# Safe master runner for the final sensitivity pipeline.
# ============================================================

R_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"

scripts <- c(
  "07_build_core_final_sensitivities.R",
  "08_build_permissive_and_latency_cohorts.R",
  "09_build_secondary_endpoints_and_bmi.R",
  "10_run_final_sensitivity_models.R"
)

for (s in scripts) {
  cat("\n\n============================================================\n")
  cat("RUNNING:", s, "\n")
  cat("============================================================\n\n")

  # Each sourced file now defines + immediately executes a stage function.
  # The DuckDB connection therefore remains valid for the full stage and is
  # safely disconnected by on.exit() only when that stage function returns.
  source(file.path(R_DIR, s), echo = FALSE, local = FALSE)
}

cat("\nAll final-sensitivity stages completed.\n")
