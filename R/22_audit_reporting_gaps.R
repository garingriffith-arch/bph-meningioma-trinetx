# Stage 22: aggregate-only inventory of remaining reporting gaps.
# No inferential models are fit and the DuckDB is opened read-only.
options(stringsAsFactors = FALSE, scipen = 999)
for (pkg in c("DBI", "duckdb", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg)
}
suppressPackageStartupMessages({ library(DBI); library(duckdb); library(data.table) })
run_stage_22 <- function() {
work <- Sys.getenv("BPH_WORK_DIR", "C:/TriNetX/bph_m_work")
db <- file.path(work, "duckdb", "bph_m.duckdb")
out <- file.path(work, "reporting_gaps_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(db)) stop("Missing DuckDB: ", db)
con <- dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
tables <- dbListTables(con)
for (tab in c("analysis_cohort_primary_boundary_fixed", "analysis_cohort_tte_index_itt_zero")) {
  if (!(tab %in% tables)) next
  fields <- dbListFields(con, tab)
  if (!("patient_regional_location" %in% fields) || !("exposure" %in% fields)) next
  q <- paste0("SELECT exposure AS arm, ",
      "COALESCE(CAST(patient_regional_location AS VARCHAR), '[NULL]') AS location_code, ",
      "COUNT(*) AS n FROM ", tab,
      " GROUP BY exposure, COALESCE(CAST(patient_regional_location AS VARCHAR), '[NULL]')")
  x <- as.data.table(dbGetQuery(con, q))
  x[n > 0 & n < 11, n := NA_integer_]
  fwrite(x, file.path(out, paste0(tab, "_location_code_counts.csv")), na = "SUPPRESSED")
}
# Locate saved diagnostics only; do not substitute QC file presence for a PH test.
scan_dirs <- file.path(work, c("analysis", "analysis_final",
  "analysis_methodologic_validity", "analysis_stage16_high_specificity", "qc", "qc_methodologic_validity"))
candidates <- unlist(lapply(scan_dirs[file.exists(scan_dirs)], function(d)
  list.files(d, pattern = "(ph_|cox.zph|proportional|diagnostic)",
             full.names = TRUE, recursive = TRUE, ignore.case = TRUE)), use.names = FALSE)
if (length(candidates)) {
  inv <- data.table(file = basename(candidates), directory = dirname(candidates),
                    bytes = as.numeric(file.info(candidates)$size))
} else inv <- data.table(file = character(), directory = character(), bytes = numeric())
fwrite(inv, file.path(out, "available_model_diagnostics_inventory.csv"))
writeLines(c(
  "Reporting gap inventory; database read-only; aggregate exports only.",
  "Location values represent source codes. DO NOT interpret 'Ex-US' geography without the data dictionary.",
  "A listed PH-diagnostics file does NOT imply that the assumption passed.",
  "Review model-specific values and preserve decisions about transformations.",
  "Table S9 rare outcomes remain separate from the boundary-aligned primary denominator until refitted.",
  "Never share person-level data or small cells. Check all CSVs before external distribution."
), file.path(out, "READ_ME_FIRST.txt"))
cat("Stage 22 reporting inventory: ", out, "\n", sep = "")

}
run_stage_22()
