# Print the source dependency order. Never executes analysis scripts.
files <- list.files("R", pattern = "^[0-9]{2}_.*\\.R$", full.names = FALSE)
files <- files[order(as.integer(sub("^([0-9]{2}).*", "\\1", files)))]
cat("Historical analysis/validation order (inspection only):\n")
cat(paste(sprintf("  %s", files), collapse = "\n"), "\n")
cat("Stages 00-17 can overwrite DuckDB tables or analytical outputs.\n")
cat("Stage 13 is superseded by Stage 17. No scripts were run.\n")
