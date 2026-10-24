# Stage 21: TARGET reporting for the day-zero initiation analysis.
# Read-only with respect to the source DuckDB. Exports aggregate data only.
# Reuses the cohort and definitions created by Stage 14; does not fit models.

options(stringsAsFactors = FALSE, scipen = 999)
required <- c("DBI", "duckdb", "data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(data.table)
})
run_stage_21 <- function() {
work_dir <- Sys.getenv("BPH_WORK_DIR", "C:/TriNetX/bph_m_work")
db_path <- file.path(work_dir, "duckdb", "bph_m.duckdb")
out_dir <- file.path(work_dir, "reporting_TARGET_day0")
if (!file.exists(db_path)) stop("DuckDB not found: ", db_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
con <- dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = TRUE)
on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
needed <- c("permissive_eligibility", "method_dx_summary",
            "analysis_cohort_tte_index_itt_zero", "candidate_observation")
absent <- setdiff(needed, dbListTables(con))
if (length(absent)) stop("Missing Stage 14 tables: ", paste(absent, collapse = ", "))

save_aggregate <- function(x, name) {
  stopifnot(!("patient_id" %in% names(x)))
  fwrite(as.data.table(x), file.path(out_dir, name), na = "")
}
count0 <- function(sql) as.numeric(dbGetQuery(con, sql)[[1L]][1L])
cutoff <- as.Date(dbGetQuery(con,
  "SELECT MAX(last_encounter_date) AS cutoff FROM candidate_observation")$cutoff[1])
if (is.na(cutoff)) stop("No cutoff in candidate_observation")
end_date <- format(cutoff - 1L, "%Y-%m-%d")

# Every step is cumulative, with the exact Stage 14 inclusion criteria and order.
criteria <- c(
  "p.exposure IN ('5ARI','URO_ALPHA')",
  "NOT (p.first_day_5ari=1 AND p.first_day_uro_alpha=1)",
  "p.sex='M'",
  "p.age_at_index BETWEEN 40 AND 89",
  "p.index_date >= DATE '2004-01-01'",
  paste0("p.index_date <= DATE '", end_date, "'"),
  "p.first_encounter_date <= p.index_date - INTERVAL 365 DAY",
  "p.n_encounters_prior_365 >= 1",
  "COALESCE(m.bph_preindex_365,0)=1",
  "COALESCE(m.prior_prostate_cancer,0)=0",
  "m.broad_preindex_or_same_day IS NULL",
  "COALESCE(p.known_finasteride_1mg,0)=0",
  "p.last_encounter_date > p.index_date"
)
labels <- c(
  "Study-class initiators", "Exclude simultaneous cross-class initiations",
  "Male sex", "Age 40–89 y", "Index date 2004 or later",
  "Index by observed data cutoff", "At least 365 d pre-index observation",
  "At least one pre-index encounter", "BPH in prior 365 d or index day",
  "No prior prostate cancer", "No meningeal tumor before or on index",
  "No recorded finasteride 1-mg initiation", "Observed after day 0"
)
flow <- rbindlist(lapply(seq_along(criteria), function(i) {
  w <- paste(criteria[seq_len(i)], collapse = " AND ")
  sql <- paste0("SELECT p.exposure AS arm, COUNT(DISTINCT p.patient_id) AS n ",
    "FROM permissive_eligibility p INNER JOIN method_dx_summary m USING (patient_id) ",
    "WHERE ", w, " GROUP BY p.exposure")
  q <- as.data.table(dbGetQuery(con, sql))
  q[, `:=`(step = i, criterion = labels[i])]
  q[, .(step, criterion, arm, n)]
}), fill = TRUE)
actual <- as.data.table(dbGetQuery(con,
  "SELECT exposure AS arm, COUNT(*) AS n, COUNT(DISTINCT patient_id) AS unique_patients FROM analysis_cohort_tte_index_itt_zero GROUP BY exposure"))
last <- flow[step == length(criteria), .(arm, expected_n = n)]
comparison <- merge(actual, last, by = "arm", all = TRUE)
comparison[, passed := n == expected_n & n == unique_patients]
save_aggregate(flow, "TARGET_day0_cumulative_flow.csv")
save_aggregate(comparison, "TARGET_day0_flow_validation.csv")
if (anyNA(comparison$passed) || !all(comparison$passed)) {
  stop("Stage 21 day-zero flow does not reconcile. Review TARGET_day0_flow_validation.csv")
}

# Baseline summaries. Pull only analytic variables into transient memory;
# do not write patient identifiers or individual-level records to disk.
fields <- dbListFields(con, "analysis_cohort_tte_index_itt_zero")
core <- c("exposure", "meningioma_event_strict", "followup_days_strict")
if (!all(core %in% fields)) stop("Day-zero table missing essential columns")
continuous <- intersect(c("age_at_index", "index_year", "n_encounters_prior_365",
                          "n_encounters_prior_1095", "n_ambulatory_prior_365",
                          "n_ed_prior_365", "n_inpatient_prior_365"), fields)
binary <- intersect(c("hypertension", "diabetes", "obesity", "ckd", "smoking",
                      "cad", "heart_failure", "cerebrovascular_disease",
                      "hyperlipidemia", "osteoporosis", "alopecia",
                      "urinary_retention", "uti", "hematuria", "bladder_stone"), fields)
categorical <- intersect(c("race", "ethnicity", "marital_status",
                           "patient_regional_location"), fields)
keep <- unique(c(core, continuous, binary, categorical))
q <- paste0("SELECT ", paste(keep, collapse = ", "),
            " FROM analysis_cohort_tte_index_itt_zero")
d <- as.data.table(dbGetQuery(con, q))
if ("index_year" %in% names(d)) d[, index_year := suppressWarnings(as.numeric(index_year))]

rows <- list()
add <- function(variable, level, group, n, detail, missing_n = NA_real_) {
  rows[[length(rows) + 1L]] <<- data.table(
    variable = variable, level = level, arm = group, n = n,
    display = detail, missing_n = missing_n)
}
for (arm in c("5ARI", "URO_ALPHA")) {
  x <- d[exposure == arm]
  N <- nrow(x); if (!N) stop("Missing treatment arm: ", arm)
  add("Participants", "", arm, N, format(N, big.mark = ","))
  y <- suppressWarnings(as.numeric(x$followup_days_strict)) / 365.25
  events <- sum(as.numeric(x$meningioma_event_strict) == 1, na.rm = TRUE)
  years <- sum(y[is.finite(y) & y >= 0])
  add("Outcome events", "", arm, events, as.character(events))
  add("Person-years", "", arm, years, sprintf("%.1f", years))
  add("Crude rate per 100,000 person-years", "", arm, NA_real_,
      if (years > 0) sprintf("%.2f", events / years * 100000) else "NA")
  add("Follow-up, years", "Median [IQR]", arm, NA_real_,
      sprintf("%.2f [%.2f–%.2f]", median(y, na.rm = TRUE),
              quantile(y, .25, na.rm = TRUE), quantile(y, .75, na.rm = TRUE)))
  for (v in continuous) {
    vals <- suppressWarnings(as.numeric(x[[v]])); ok <- is.finite(vals)
    add(v, "Mean (SD)", arm, sum(ok), if (any(ok))
      sprintf("%.2f (%.2f)", mean(vals[ok]), sd(vals[ok])) else "NA", sum(!ok))
  }
  for (v in binary) {
    vals <- suppressWarnings(as.numeric(x[[v]]))
    good <- !is.na(vals); n1 <- sum(vals == 1, na.rm = TRUE)
    add(v, "Documented", arm, n1,
        if (n1 > 0 && n1 < 11) "Suppressed (<11)" else
          sprintf("%d (%.1f%%)", n1, if (sum(good)) 100*n1 / sum(good) else NA_real_), sum(!good))
  }
  for (v in categorical) {
    val <- as.character(x[[v]]); val[is.na(val) | !nzchar(trimws(val))] <- "[missing]"
    ncat <- sort(table(val), decreasing = TRUE)
    for (lab in names(ncat)) {
      n <- unname(as.integer(ncat[[lab]]))
      add(v, lab, arm, n,
          if (n > 0 && n < 11) "Suppressed (<11)" else sprintf("%d (%.1f%%)",n,100*n/N),
          sum(val == "[missing]"))
    }
  }
}
baseline <- rbindlist(rows, fill = TRUE)
# Conservative cell suppression: raw counts less than 11 are not exported.
baseline[n > 0 & n < 11, n := NA_real_]
baseline[missing_n > 0 & missing_n < 11, missing_n := NA_real_]
save_aggregate(baseline, "TARGET_day0_baseline_by_arm.csv")

context <- data.table(
  quantity = c("Observed data cutoff", "Day-zero analysis end of index inclusion",
    "Day-zero analytical patients", "Published Stage 20 expectation",
    "Flow matches final cohort", "Database access"),
  value = c(as.character(cutoff), end_date, as.character(nrow(d)), "660012 (verify against original output)",
    as.character(all(comparison$passed)), "Read-only")
)
save_aggregate(context, "TARGET_day0_run_context.csv")
writeLines(c(
  "TARGET day-zero reporting: SOURCE QC",
  paste("Date:", Sys.Date()),
  paste("Observed cutoff:", cutoff),
  paste("Analytical cohort:", nrow(d)),
  paste("Cumulative flow matches source:", all(comparison$passed)),
  "Outputs are aggregate only; small cells suppressed in baseline display.",
  "Cohort flow reflects the order of Stage 14 inclusion criteria.",
  "These are observed cohort counts, not a randomized target trial.",
  "Do not publish until verified against institutional disclosure requirements."
), file.path(out_dir, "TARGET_day0_QC_README.txt"))
cat("Stage 21 finished: ", out_dir, "\n", sep = "")

}
run_stage_21()
