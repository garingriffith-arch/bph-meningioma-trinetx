# ============================================================
# 13_make_manuscript_tables_figures.R
#
# Finasteride / Dutasteride and Incident Meningioma
# Manuscript-ready tables + figures
# PATCH v4: final manuscript review; removed uncorrected confirmed-code sensitivity and separated Figure 3 estimate labels
#
# Produces:
#   Table 1. Baseline characteristics
#   Table 2. Primary and robustness analyses
#   Table 3. Exploratory persistence / drug-specific analyses
#   Table S1. Covariate balance
#   Table S2. Latency analyses
#   Table S3. High-specificity / clinically consequential outcomes
#
#   Figure 1. Cohort flow
#   Figure 2. Primary / robustness forest plot
#   Figure 3. Exploratory persistence / drug-specific forest plot
#   Figure S1. Love plot
#
# Exports:
#   - editable CSVs
#   - a Word document containing all tables
#   - PNG (600 dpi), TIFF (600 dpi/LZW), and vector PDF figures
#
# Word formatting:
#   Table titles: Times New Roman 12
#   Table body/notes: Arial 10 / Arial 8
#
# IMPORTANT:
# This script DOES NOT refit the inferential models in Stages 06/10/12
# except that it re-estimates the primary overlap weights to construct
# Table 1 descriptive summaries. Inferential estimates are read directly
# from the finalized stage outputs.
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

# ------------------------------------------------------------
# 0. PATHS
# ------------------------------------------------------------
R_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"
WORK_DIR <- "C:/TriNetX/bph_m_work"
DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")

ANALYSIS_DIR <- file.path(WORK_DIR, "analysis")
FINAL_DIR <- file.path(WORK_DIR, "analysis_final")
SIGNAL_DIR <- file.path(WORK_DIR, "analysis_signal_interrogation")
QC_DIR <- file.path(WORK_DIR, "qc")

OUT_DIR <- file.path(WORK_DIR, "manuscript_ready")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR <- file.path(OUT_DIR, "figures")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------
required_pkgs <- c(
  "DBI",
  "duckdb",
  "data.table",
  "survival",
  "splines",
  "ggplot2",
  "scales",
  "flextable",
  "officer"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(data.table)
  library(survival)
  library(splines)
  library(ggplot2)
  library(scales)
  library(flextable)
  library(officer)
})

HELPER_FILE <- file.path(R_DIR, "final_analysis_helpers.R")
if (!file.exists(HELPER_FILE)) {
  stop(
    "Missing final_analysis_helpers.R at:\n",
    HELPER_FILE,
    "\nCopy it from the final sensitivity-analysis package into R_DIR."
  )
}
source(HELPER_FILE)

# ------------------------------------------------------------
# 2. FILE HELPERS
# ------------------------------------------------------------
first_existing <- function(paths, required = TRUE) {
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0L) return(hit[1])
  if (required) {
    stop(
      "Could not find any of:\n",
      paste0("  - ", paths, collapse = "\n")
    )
  }
  NA_character_
}

read_csv_required <- function(paths) {
  fread(first_existing(paths))
}

p06_file <- first_existing(c(
  file.path(ANALYSIS_DIR, "06_model_results.csv"),
  file.path(WORK_DIR, "06_model_results.csv")
))
inc06_file <- first_existing(c(
  file.path(ANALYSIS_DIR, "06_primary_incidence.csv"),
  file.path(WORK_DIR, "06_primary_incidence.csv")
))
bal06_file <- first_existing(c(
  file.path(ANALYSIS_DIR, "06_balance_all.csv"),
  file.path(WORK_DIR, "06_balance_all.csv")
))

final_model_file <- first_existing(c(
  file.path(FINAL_DIR, "10_final_model_results.csv")
))
final_inc_file <- first_existing(c(
  file.path(FINAL_DIR, "10_final_incidence_rates.csv")
))
latency_file <- first_existing(c(
  file.path(FINAL_DIR, "10_latency_interval_results.csv")
))
competing_file <- first_existing(c(
  file.path(FINAL_DIR, "10_competing_risk_result.csv")
))

signal_model_file <- first_existing(c(
  file.path(SIGNAL_DIR, "12_signal_model_results.csv")
))
signal_inc_file <- first_existing(c(
  file.path(SIGNAL_DIR, "12_signal_incidence_rates.csv")
))

flow_file <- first_existing(c(
  file.path(QC_DIR, "02_cohort_flow.csv"),
  file.path(WORK_DIR, "02_cohort_flow.csv")
))
revised_base_file <- first_existing(c(
  file.path(QC_DIR, "05_revised_base_counts.csv"),
  file.path(WORK_DIR, "05_revised_base_counts.csv")
))
revised_primary_file <- first_existing(c(
  file.path(QC_DIR, "05_revised_primary_counts.csv"),
  file.path(WORK_DIR, "05_revised_primary_counts.csv")
))

m06 <- fread(p06_file)
inc06 <- fread(inc06_file)
bal06 <- fread(bal06_file)
m10 <- fread(final_model_file)
inc10 <- fread(final_inc_file)
latency <- fread(latency_file)
competing <- fread(competing_file)
m12 <- fread(signal_model_file)
inc12 <- fread(signal_inc_file)

# ------------------------------------------------------------
# 3. MANUSCRIPT FORMATTING HELPERS
# ------------------------------------------------------------
fmt_n <- function(x) {
  ifelse(is.na(x), "\u2014", format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE))
}

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "\u2014", sprintf(paste0("%.", digits, "f"), x))
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    "\u2014",
    ifelse(
      p < 0.001,
      "<0.001",
      ifelse(p > 0.9, ">0.9", sprintf("%.3f", p))
    )
  )
}

fmt_hr_ci <- function(hr, lo, hi) {
  ifelse(
    is.na(hr) | is.na(lo) | is.na(hi),
    "\u2014",
    sprintf("%.2f (%.2f\u2013%.2f)", hr, lo, hi)
  )
}

fmt_rate <- function(x) {
  ifelse(is.na(x), "\u2014", sprintf("%.1f", x))
}

fmt_mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return("\u2014")
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

fmt_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return("\u2014")
  q <- quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE)
  sprintf("%.1f [%.1f\u2013%.1f]", q[2], q[1], q[3])
}

fmt_n_pct <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return("\u2014")
  n <- sum(x == 1, na.rm = TRUE)
  p <- 100 * mean(x == 1, na.rm = TRUE)
  sprintf("%s (%.1f)", format(n, big.mark = ","), p)
}

safe_smd <- function(x, z, w) {
  out <- try(smd_numeric(as.numeric(x), z = z, w = w), silent = TRUE)
  if (inherits(out, "try-error") || !is.finite(out)) return(NA_real_)
  abs(as.numeric(out))
}

pretty_covariate <- function(x) {
  x <- gsub("^race_f", "Race: ", x)
  x <- gsub("^ethnicity_f", "Ethnicity: ", x)
  x <- gsub("^marital_f", "Marital status: ", x)
  x <- gsub("^region_f", "Region: ", x)
  x <- gsub("^source_f", "Source: ", x)

  lookup <- c(
    age_at_index = "Age",
    index_year = "Index year",
    hypertension = "Hypertension",
    diabetes = "Diabetes",
    obesity = "Obesity",
    ckd = "Chronic kidney disease",
    smoking = "Smoking",
    cad = "Coronary artery disease",
    heart_failure = "Heart failure",
    cerebrovascular_disease = "Cerebrovascular disease",
    hyperlipidemia = "Hyperlipidemia",
    osteoporosis = "Osteoporosis",
    alopecia = "Alopecia",
    urinary_retention = "Urinary retention",
    uti = "Urinary tract infection",
    hematuria = "Hematuria",
    bladder_stone = "Bladder stone",
    log_enc_365 = "Encounters in prior year",
    log_enc_1095 = "Encounters in prior 3 years",
    log_amb_365 = "Ambulatory encounters in prior year",
    log_ed_365 = "ED encounters in prior year",
    log_ip_365 = "Inpatient encounters in prior year"
  )

  idx <- match(x, names(lookup))
  x[!is.na(idx)] <- unname(lookup[idx[!is.na(idx)]])
  x
}

# ------------------------------------------------------------
# 4. WORD TABLE STYLE
# ------------------------------------------------------------
make_ft <- function(dt, widths = NULL) {

  ft <- flextable(dt)
  ft <- theme_booktabs(ft)
  ft <- font(ft, fontname = "Arial", part = "all")
  ft <- fontsize(ft, size = 10, part = "body")
  ft <- fontsize(ft, size = 10, part = "header")
  ft <- bold(ft, bold = TRUE, part = "header")
  ft <- align(ft, align = "left", part = "all")
  ft <- valign(ft, valign = "center", part = "all")
  ft <- padding(
    ft,
    padding.top = 2,
    padding.bottom = 2,
    padding.left = 3,
    padding.right = 3,
    part = "all"
  )
  ft <- set_table_properties(
    ft,
    layout = "autofit",
    width = 1
  )

  if (!is.null(widths)) {
    for (nm in intersect(names(widths), colnames(dt))) {
      ft <- width(ft, j = nm, width = widths[[nm]])
    }
  }

  ft
}

title_fpar <- function(number, title) {
  fpar(
    ftext(
      paste0("Table ", number, ". "),
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 12,
        italic = TRUE
      )
    ),
    ftext(
      title,
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 12
      )
    )
  )
}

note_fpar <- function(note) {
  fpar(
    ftext(
      note,
      prop = fp_text(
        font.family = "Arial",
        font.size = 8
      )
    )
  )
}

add_table_block <- function(doc, number, title, ft, note = NULL) {
  doc <- body_add_fpar(doc, title_fpar(number, title))
  doc <- body_add_flextable(doc, ft)
  if (!is.null(note)) {
    doc <- body_add_fpar(doc, note_fpar(note))
  }
  doc <- body_add_par(doc, "")
  doc
}

# ------------------------------------------------------------
# 5. TABLE 1: BASELINE CHARACTERISTICS
# ------------------------------------------------------------
cat("\nBuilding primary overlap weights for Table 1...\n")

con <- dbConnect(
  duckdb::duckdb(),
  dbdir = DB_PATH,
  read_only = TRUE
)
on.exit(
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE),
  add = TRUE
)

primary <- as.data.table(
  dbGetQuery(
    con,
    "SELECT * FROM analysis_cohort_primary_revised;"
  )
)

primary[, trt := as.integer(exposure == "5ARI")]

wt <- estimate_overlap_weights(
  primary,
  label = "Manuscript Table 1 overlap weights"
)
d <- wt$data

z <- d$trt
wu <- rep(1, nrow(d))
ww <- d$ow

table1_rows <- list()

add_cont_row <- function(label, var, summary = c("mean_sd", "median_iqr")) {
  summary <- match.arg(summary)

  x <- as.numeric(d[[var]])

  v1 <- if (summary == "mean_sd") {
    fmt_mean_sd(x[z == 1])
  } else {
    fmt_median_iqr(x[z == 1])
  }

  v0 <- if (summary == "mean_sd") {
    fmt_mean_sd(x[z == 0])
  } else {
    fmt_median_iqr(x[z == 0])
  }

  table1_rows[[length(table1_rows) + 1L]] <<- data.table(
    Characteristic = label,
    `5-ARI` = v1,
    `Uroselective alpha-blocker` = v0,
    `SMD before` = sprintf("%.3f", safe_smd(x, z, wu)),
    `SMD after` = sprintf("%.3f", safe_smd(x, z, ww))
  )
}

add_binary_row <- function(label, var) {
  x <- as.integer(d[[var]])

  table1_rows[[length(table1_rows) + 1L]] <<- data.table(
    Characteristic = label,
    `5-ARI` = fmt_n_pct(x[z == 1]),
    `Uroselective alpha-blocker` = fmt_n_pct(x[z == 0]),
    `SMD before` = sprintf("%.3f", safe_smd(x, z, wu)),
    `SMD after` = sprintf("%.3f", safe_smd(x, z, ww))
  )
}

add_factor_block <- function(label, var, preferred_levels = NULL) {

  x <- as.character(d[[var]])
  x[is.na(x) | trimws(x) == ""] <- "Unknown"

  levs <- if (is.null(preferred_levels)) {
    sort(unique(x))
  } else {
    unique(c(
      preferred_levels[preferred_levels %in% unique(x)],
      setdiff(sort(unique(x)), preferred_levels)
    ))
  }

  table1_rows[[length(table1_rows) + 1L]] <<- data.table(
    Characteristic = label,
    `5-ARI` = "",
    `Uroselective alpha-blocker` = "",
    `SMD before` = "",
    `SMD after` = ""
  )

  for (lv in levs) {
    b <- as.integer(x == lv)

    table1_rows[[length(table1_rows) + 1L]] <<- data.table(
      Characteristic = paste0("   ", lv),
      `5-ARI` = fmt_n_pct(b[z == 1]),
      `Uroselective alpha-blocker` = fmt_n_pct(b[z == 0]),
      `SMD before` = sprintf("%.3f", safe_smd(b, z, wu)),
      `SMD after` = sprintf("%.3f", safe_smd(b, z, ww))
    )
  }
}

add_cont_row("Age, mean (SD), y", "age_at_index", "mean_sd")
add_cont_row("Index year, median [IQR]", "index_year", "median_iqr")

add_factor_block(
  "Race, No. (%)",
  "race",
  c(
    "White",
    "Black or African American",
    "Asian",
    "American Indian or Alaska Native",
    "Native Hawaiian or Other Pacific Islander",
    "Other Race",
    "Unknown"
  )
)

add_factor_block(
  "Ethnicity, No. (%)",
  "ethnicity",
  c(
    "Not Hispanic or Latino",
    "Hispanic or Latino",
    "Unknown"
  )
)

add_binary_row("Hypertension, No. (%)", "hypertension")
add_binary_row("Diabetes, No. (%)", "diabetes")
add_binary_row("Obesity, No. (%)", "obesity")
add_binary_row("Chronic kidney disease, No. (%)", "ckd")
add_binary_row("Smoking, No. (%)", "smoking")
add_binary_row("Coronary artery disease, No. (%)", "cad")
add_binary_row("Heart failure, No. (%)", "heart_failure")
add_binary_row("Cerebrovascular disease, No. (%)", "cerebrovascular_disease")
add_binary_row("Hyperlipidemia, No. (%)", "hyperlipidemia")
add_binary_row("Osteoporosis, No. (%)", "osteoporosis")
add_binary_row("Alopecia, No. (%)", "alopecia")
add_binary_row("Urinary retention, No. (%)", "urinary_retention")
add_binary_row("Urinary tract infection, No. (%)", "uti")
add_binary_row("Hematuria, No. (%)", "hematuria")
add_binary_row("Bladder stone, No. (%)", "bladder_stone")

add_cont_row(
  "Encounters in prior year, median [IQR]",
  "n_encounters_prior_365",
  "median_iqr"
)
add_cont_row(
  "Ambulatory encounters in prior year, median [IQR]",
  "n_ambulatory_prior_365",
  "median_iqr"
)
add_cont_row(
  "ED encounters in prior year, median [IQR]",
  "n_ed_prior_365",
  "median_iqr"
)
add_cont_row(
  "Inpatient encounters in prior year, median [IQR]",
  "n_inpatient_prior_365",
  "median_iqr"
)

table1 <- rbindlist(table1_rows, fill = TRUE)

setnames(
  table1,
  c("5-ARI", "Uroselective alpha-blocker"),
  c(
    paste0("5-ARI\n(n=", fmt_n(sum(z == 1)), ")"),
    paste0(
      "Uroselective alpha-blocker\n(n=",
      fmt_n(sum(z == 0)),
      ")"
    )
  )
)

fwrite(table1, file.path(TABLE_DIR, "Table_1_Baseline_characteristics.csv"))

# ------------------------------------------------------------
# 6. NORMALIZE MODEL OUTPUTS
# ------------------------------------------------------------
normalize_models <- function(dt, source) {
  x <- copy(dt)

  if ("n_5ari" %in% names(x)) {
    setnames(
      x,
      c(
        "n_5ari",
        "n_uro_alpha",
        "events_5ari",
        "events_uro_alpha",
        "ESS_5ari",
        "ESS_uro_alpha"
      ),
      c(
        "n_treated",
        "n_control",
        "events_treated",
        "events_control",
        "ESS_treated",
        "ESS_control"
      ),
      skip_absent = TRUE
    )
  }

  x[, source := source]
  x
}

nm06 <- normalize_models(m06, "Stage 06")
nm10 <- normalize_models(m10, "Stage 10")
nm12 <- normalize_models(m12, "Stage 12")

inc06n <- copy(inc06)
inc06n[
  ,
  `:=`(
    analysis = "Primary: 1-y lag, strict meningioma, all 5ARI vs uro-alpha",
    trt = as.integer(exposure == "5ARI")
  )
]

get_rate <- function(inc_dt, analysis_name, trt_value) {
  x <- inc_dt[
    analysis == analysis_name &
      trt == trt_value
  ]
  if (nrow(x) == 0L) return(NA_real_)
  x$weighted_rate_per_100k_py[1]
}

extract_model <- function(
    dt,
    analysis_name,
    display,
    incidence_dt = NULL,
    treated_label = "5-ARI",
    control_label = "Uro-alpha"
) {

  x <- dt[analysis == analysis_name]

  if (nrow(x) != 1L) {
    stop(
      "Expected one result for:\n",
      analysis_name,
      "\nFound: ",
      nrow(x)
    )
  }

  rt <- if (!is.null(incidence_dt)) {
    get_rate(incidence_dt, analysis_name, 1)
  } else {
    NA_real_
  }

  rc <- if (!is.null(incidence_dt)) {
    get_rate(incidence_dt, analysis_name, 0)
  } else {
    NA_real_
  }

  data.table(
    Analysis = display,
    N = x$n,
    Exposed = treated_label,
    Comparator = control_label,
    Events_exposed = x$events_treated,
    Events_comparator = x$events_control,
    Rate_exposed = rt,
    Rate_comparator = rc,
    HR = x$weighted_HR,
    CI_low = x$weighted_CI_low,
    CI_high = x$weighted_CI_high,
    p = x$weighted_p
  )
}

# ------------------------------------------------------------
# 7. TABLE 2: PRIMARY + ROBUSTNESS ANALYSES
# ------------------------------------------------------------
t2_raw <- rbindlist(
  list(
    extract_model(
      nm06,
      "Primary: 1-y lag, strict meningioma, all 5ARI vs uro-alpha",
      "Primary: 1-y landmark, strict meningioma",
      inc06n
    ),
    extract_model(
      nm06,
      "Sensitivity: 5-y administrative horizon",
      "5-y administrative follow-up",
      NULL
    ),
    extract_model(
      nm06,
      "Sensitivity: confirmed meningioma (>=2 codes >=30 d apart)",
      "Confirmed meningioma (\u22652 codes \u226530 d apart)",
      NULL
    ),
    extract_model(
      nm06,
      "Sensitivity: 2-y lag",
      "2-y landmark",
      NULL
    ),
    extract_model(
      nm10,
      "High-specificity intracranial meningioma: D32.0 + ICD-9 225.2",
      "High-specificity intracranial meningioma",
      inc10
    ),
    extract_model(
      nm10,
      "ICD-10-era D32.0-only meningioma",
      "ICD-10-era D32.0 only",
      inc10
    ),
    extract_model(
      nm10,
      "Measured-BMI sensitivity (median imputation + missingness indicator)",
      "Measured-BMI adjusted sensitivity",
      inc10
    )
  ),
  fill = TRUE
)

# Replace incorrectly inherited primary rates for rows whose dedicated
# incidence was not stored in Stage 06 with em dash rather than imply equality.
t2_raw[
  Analysis %in% c(
    "5-y administrative follow-up",
    "Confirmed meningioma (\u22652 codes \u226530 d apart)",
    "2-y landmark"
  ),
  `:=`(
    Rate_exposed = NA_real_,
    Rate_comparator = NA_real_
  )
]

# Add Fine-Gray competing-risk sensitivity.
if (nrow(competing) >= 1L) {
  t2_raw <- rbind(
    t2_raw,
    data.table(
      Analysis = "Fine-Gray model with death as competing event",
      N = competing$n[1],
      Exposed = "5-ARI",
      Comparator = "Uro-alpha",
      Events_exposed = NA_real_,
      Events_comparator = NA_real_,
      Rate_exposed = NA_real_,
      Rate_comparator = NA_real_,
      HR = competing$subdistribution_HR[1],
      CI_low = competing$CI_low[1],
      CI_high = competing$CI_high[1],
      p = competing$p[1]
    ),
    fill = TRUE
  )
}

table2 <- t2_raw[
  ,
  .(
    Analysis,
    `N` = fmt_n(N),
    `Events, exposed/comparator` = ifelse(
      is.na(Events_exposed),
      "\u2014",
      paste0(fmt_n(Events_exposed), "/", fmt_n(Events_comparator))
    ),
    `Weighted incidence/100,000 PY, exposed` = fmt_rate(Rate_exposed),
    `Weighted incidence/100,000 PY, comparator` = fmt_rate(Rate_comparator),
    `Adjusted HR (95% CI)` = fmt_hr_ci(HR, CI_low, CI_high),
    `P value` = fmt_p(p)
  )
]

fwrite(table2, file.path(TABLE_DIR, "Table_2_Primary_robustness_results.csv"))

# ------------------------------------------------------------
# 8. TABLE 3: EXPLORATORY PERSISTENCE / DRUG-SPECIFIC
# ------------------------------------------------------------
t3_specs <- list(
  c(
    "Sustained 1-y + landmark utilization: all 5ARI vs uro-alpha",
    "Sustained 1-y: all 5-ARI vs uro-alpha",
    "5-ARI",
    "Uro-alpha"
  ),
  c(
    "Sustained 1-y + landmark utilization: finasteride vs uro-alpha",
    "Sustained 1-y: finasteride vs uro-alpha",
    "Finasteride",
    "Uro-alpha"
  ),
  c(
    "Sustained 1-y + landmark utilization: dutasteride vs uro-alpha",
    "Sustained 1-y: dutasteride vs uro-alpha",
    "Dutasteride",
    "Uro-alpha"
  ),
  c(
    "Strict sustained 1-y + landmark utilization: all 5ARI vs uro-alpha",
    "Strict sustained 1-y: all 5-ARI vs uro-alpha",
    "5-ARI",
    "Uro-alpha"
  ),
  c(
    "Strict sustained 1-y + landmark utilization: finasteride vs uro-alpha",
    "Strict sustained 1-y: finasteride vs uro-alpha",
    "Finasteride",
    "Uro-alpha"
  ),
  c(
    "Strict sustained 1-y + landmark utilization: dutasteride vs uro-alpha",
    "Strict sustained 1-y: dutasteride vs uro-alpha",
    "Dutasteride",
    "Uro-alpha"
  ),
  c(
    "Persistent 3-y + landmark utilization: all 5ARI vs uro-alpha",
    "Persistent 3-y: all 5-ARI vs uro-alpha",
    "5-ARI",
    "Uro-alpha"
  ),
  c(
    "Persistent 3-y + landmark utilization: finasteride vs uro-alpha",
    "Persistent 3-y: finasteride vs uro-alpha",
    "Finasteride",
    "Uro-alpha"
  ),
  c(
    "Persistent 3-y + landmark utilization: dutasteride vs uro-alpha",
    "Persistent 3-y: dutasteride vs uro-alpha",
    "Dutasteride",
    "Uro-alpha"
  ),
  c(
    "Head-to-head homogeneous: US, 2012+, exact BPH doses, finasteride vs dutasteride",
    "US 2012+, exact dose: finasteride vs dutasteride",
    "Finasteride",
    "Dutasteride"
  ),
  c(
    "Head-to-head homogeneous: US, ICD-10 era, exact BPH doses",
    "US ICD-10 era, exact dose: finasteride vs dutasteride",
    "Finasteride",
    "Dutasteride"
  )
)

t3_raw <- rbindlist(
  lapply(
    t3_specs,
    function(s) {
      extract_model(
        nm12,
        s[1],
        s[2],
        inc12,
        treated_label = s[3],
        control_label = s[4]
      )
    }
  ),
  fill = TRUE
)

table3 <- t3_raw[
  ,
  .(
    Analysis,
    `N` = fmt_n(N),
    `Events, exposed/comparator` = paste0(
      fmt_n(Events_exposed),
      "/",
      fmt_n(Events_comparator)
    ),
    `Weighted incidence/100,000 PY, exposed` = fmt_rate(Rate_exposed),
    `Weighted incidence/100,000 PY, comparator` = fmt_rate(Rate_comparator),
    `Adjusted HR (95% CI)` = fmt_hr_ci(HR, CI_low, CI_high),
    `P value` = fmt_p(p)
  )
]

fwrite(table3, file.path(TABLE_DIR, "Table_3_Exploratory_persistence_results.csv"))

# ------------------------------------------------------------
# 9. SUPPLEMENTARY TABLE S1: COVARIATE BALANCE
# ------------------------------------------------------------
primary_balance <- bal06[
  analysis ==
    "Primary: 1-y lag, strict meningioma, all 5ARI vs uro-alpha"
]

tableS1 <- primary_balance[
  ,
  .(
    Covariate = pretty_covariate(covariate),
    `SMD before weighting` = sprintf("%.3f", abs_smd_unweighted),
    `SMD after overlap weighting` = sprintf("%.3f", abs_smd_weighted)
  )
][order(-as.numeric(`SMD before weighting`))]

fwrite(tableS1, file.path(TABLE_DIR, "Table_S1_Covariate_balance.csv"))

# ------------------------------------------------------------
# 10. SUPPLEMENTARY TABLE S2: LATENCY
# ------------------------------------------------------------
latency[, interval_display := fcase(
  interval == "0-<1 y", "0–<1 y",
  interval == "1-<2 y", "1–<2 y",
  interval == "2-<5 y", "2–<5 y",
  interval == ">=5 y", "≥5 y",
  default = interval
)]

tableS2 <- latency[
  ,
  .(
    `Time since initiation` = interval_display,
    `N at interval entry` = fmt_n(n),
    `Events` = fmt_n(events),
    `Events, 5-ARI/uro-alpha` = paste0(
      fmt_n(events_5ari),
      "/",
      fmt_n(events_uro_alpha)
    ),
    `Adjusted HR (95% CI)` = fmt_hr_ci(HR, CI_low, CI_high),
    `P value` = fmt_p(p)
  )
]

fwrite(tableS2, file.path(TABLE_DIR, "Table_S2_Latency_results.csv"))

# ------------------------------------------------------------
# 11. SUPPLEMENTARY TABLE S3: RARE / HARD ENDPOINTS
# ------------------------------------------------------------
rare_labels <- c(
  "ICD-O tumor-table meningioma",
  "Cancer-registry-source ICD-O meningioma",
  "Clinically treated meningioma: resection or cranial SRS within 180 d",
  "Oncology-treatment-table ICD-O meningioma"
)

rare_display <- c(
  "ICD-O tumor-table meningioma",
  "Cancer Registry ICD-O meningioma",
  "Meningioma followed by resection/SRS within 180 d",
  "Oncology-treatment-table ICD-O meningioma"
)

rare_raw <- rbindlist(
  Map(
    function(a, dlabel) {
      extract_model(
        nm10,
        a,
        dlabel,
        inc10
      )
    },
    rare_labels,
    rare_display
  ),
  fill = TRUE
)

tableS3 <- rare_raw[
  ,
  .(
    Analysis,
    `N` = fmt_n(N),
    `Events` = fmt_n(Events_exposed + Events_comparator),
    `Adjusted HR (95% CI)` = fmt_hr_ci(HR, CI_low, CI_high),
    `P value` = fmt_p(p)
  )
]

fwrite(tableS3, file.path(TABLE_DIR, "Table_S3_Rare_hard_endpoints.csv"))

# ------------------------------------------------------------
# 12. WORD DOCUMENT WITH ALL TABLES
# ------------------------------------------------------------
ft1 <- make_ft(table1)
ft2 <- make_ft(table2)
ft3 <- make_ft(table3)
ftS1 <- make_ft(tableS1)
ftS2 <- make_ft(tableS2)
ftS3 <- make_ft(tableS3)

# Bold the factor-block labels in Table 1.
factor_headers <- which(
  table1$Characteristic %in% c(
    "Race, No. (%)",
    "Ethnicity, No. (%)"
  )
)
if (length(factor_headers) > 0L) {
  ft1 <- bold(ft1, i = factor_headers, bold = TRUE, part = "body")
}

doc <- read_docx()

doc <- add_table_block(
  doc,
  "1",
  "Baseline Characteristics of the Primary Active-Comparator Cohort",
  ft1,
  paste0(
    "Values are No. (%) unless otherwise indicated. ",
    "SMD indicates standardized mean difference. ",
    "Overlap weighting was estimated using the finalized primary propensity-score model. ",
    "The propensity score additionally included marital status and geographic region, ",
    "which are reported in Supplementary Table 1."
  )
)

doc <- add_table_block(
  doc,
  "2",
  "Primary and Robustness Analyses of Incident Meningioma",
  ft2,
  paste0(
    "HRs are overlap-weighted Cox proportional hazards estimates except for the ",
    "Fine-Gray analysis, which reports a subdistribution HR with death as a competing event. ",
    "PY indicates person-years."
  )
)

doc <- add_table_block(
  doc,
  "3",
  "Exploratory Sustained-Exposure and Drug-Specific Analyses",
  ft3,
  paste0(
    "These analyses were exploratory. Sustained/persistent landmark models incorporated ",
    "healthcare utilization measured during the exposure-assessment period into the ",
    "propensity-score model. Outcome follow-up began after the corresponding landmark."
  )
)

doc <- add_table_block(
  doc,
  "S1",
  "Covariate Balance Before and After Overlap Weighting",
  ftS1,
  paste0(
    "Absolute standardized mean differences are shown. ",
    "Overlap weighting provides exact mean balance for covariates included in a ",
    "logistic propensity-score model up to numerical precision."
  )
)

doc <- add_table_block(
  doc,
  "S2",
  "Interval-Specific Analyses by Time Since Treatment Initiation",
  ftS2,
  "HR indicates overlap-weighted hazard ratio."
)

doc <- add_table_block(
  doc,
  "S3",
  "High-Specificity and Clinically Consequential Meningioma Endpoints",
  ftS3,
  paste0(
    "Rare-event analyses are presented as supportive sensitivity analyses and ",
    "were not used to define the primary outcome."
  )
)

print(
  doc,
  target = file.path(TABLE_DIR, "Manuscript_Tables.docx")
)


# Figure font compatibility:
# Use the generic "sans" family for graphics. On Windows this avoids the
# PostScript/PDF error caused by requesting "Arial" from a device that does
# not have Arial registered in its PostScript font database. Word tables
# remain explicitly formatted in Arial.

# ------------------------------------------------------------
# 13. FIGURE EXPORT HELPER
# ------------------------------------------------------------
save_figure <- function(
    p,
    stem,
    width,
    height
) {

  ggsave(
    file.path(FIG_DIR, paste0(stem, ".png")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )

  ggsave(
    file.path(FIG_DIR, paste0(stem, ".tiff")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )

  ggsave(
    file.path(FIG_DIR, paste0(stem, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "white"
  )
}

theme_manuscript <- function(base_size = 10) {
  theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 1,
        margin = margin(b = 8)
      ),
      axis.text = element_text(size = base_size),
      axis.title = element_text(size = base_size),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(8, 18, 8, 8)
    )
}

# ------------------------------------------------------------
# 14. FIGURE 1: COHORT FLOW
# ------------------------------------------------------------
flow <- fread(flow_file)
revbase <- fread(revised_base_file)
revprimary <- fread(revised_primary_file)

get_flow_n <- function(step_number) {
  flow[step == step_number, total][1]
}

n_initial <- get_flow_n(1)
n_calendar <- get_flow_n(4)
n_history <- get_flow_n(5)
n_bph <- get_flow_n(6)
n_after_core_exclusions <- get_flow_n(10)
n_post_landmark <- sum(revbase$patients)
n_final <- sum(revprimary$patients)
n_final_5ari <- revprimary[exposure == "5ARI", patients]
n_final_uro <- revprimary[exposure == "URO_ALPHA", patients]

flow_nodes <- data.table(
  y = 7:1,
  label = c(
    paste0(
      "Initial unambiguous 5-ARI or uroselective\n",
      "alpha-blocker initiators\nN = ",
      fmt_n(n_initial)
    ),
    paste0(
      "Age 40\u201389 y and calendar eligible\nN = ",
      fmt_n(n_calendar)
    ),
    paste0(
      "\u2265365 d observable history and encounter\n",
      "in the prior year\nN = ",
      fmt_n(n_history)
    ),
    paste0(
      "BPH diagnosis near treatment initiation\nN = ",
      fmt_n(n_bph)
    ),
    paste0(
      "After exclusions for prior prostate cancer,\n",
      "prior meningeal tumor,\n",
      "early opposite-class use, and known\n",
      "finasteride 1-mg initiation\nN = ",
      fmt_n(n_after_core_exclusions)
    ),
    paste0(
      "Observed beyond the 1-y landmark\nN = ",
      fmt_n(n_post_landmark)
    ),
    paste0(
      "Final primary 1-y landmark cohort\nN = ",
      fmt_n(n_final),
      "\n5-ARI: ",
      fmt_n(n_final_5ari),
      "   |   Uro-alpha: ",
      fmt_n(n_final_uro)
    )
  )
)

flow_exclusions <- data.table(
  y = c(6.5, 5.5, 4.5, 3.5, 2.5, 1.5),
  label = c(
    paste0("Excluded: ", fmt_n(n_initial - n_calendar)),
    paste0("Excluded: ", fmt_n(n_calendar - n_history)),
    paste0("Excluded: ", fmt_n(n_history - n_bph)),
    paste0("Excluded: ", fmt_n(n_bph - n_after_core_exclusions)),
    paste0("Excluded: ", fmt_n(n_after_core_exclusions - n_post_landmark)),
    paste0(
      "Lag-period meningeal tumor / no eligible risk time: ",
      fmt_n(n_post_landmark - n_final)
    )
  )
)

p_flow <- ggplot() +
  geom_rect(
    data = flow_nodes,
    aes(
      xmin = 0.05,
      xmax = 0.95,
      ymin = y - 0.34,
      ymax = y + 0.34
    ),
    fill = "white",
    color = "black",
    linewidth = 0.45
  ) +
  geom_text(
    data = flow_nodes,
    aes(x = 0.5, y = y, label = label),
    family = "sans",
    size = 3.05,
    lineheight = 0.95
  ) +
  geom_segment(
    data = data.table(
      y1 = 6.69:1.69,
      y2 = 6.31:1.31
    ),
    aes(
      x = 0.5,
      xend = 0.5,
      y = y1,
      yend = y2
    ),
    linewidth = 0.4,
    arrow = grid::arrow(
      length = grid::unit(0.08, "inches"),
      type = "closed"
    )
  ) +
  geom_text(
    data = flow_exclusions,
    aes(x = 0.99, y = y, label = label),
    hjust = 0,
    family = "sans",
    size = 2.6
  ) +
  coord_cartesian(
    xlim = c(0, 1.62),
    ylim = c(0.55, 7.45),
    clip = "off"
  ) +
  theme_void(base_family = "sans") +
  theme(
    plot.margin = margin(10, 150, 10, 10)
  )

save_figure(
  p_flow,
  "Figure_1_Cohort_flow",
  width = 7.2,
  height = 8.8
)

# ------------------------------------------------------------
# 15. FIGURE 2: PRIMARY / ROBUSTNESS FOREST
# ------------------------------------------------------------
f2 <- copy(t2_raw)
f2[, label := Analysis]
f2[, estimate_text := fmt_hr_ci(HR, CI_low, CI_high)]

f2[, label := factor(label, levels = rev(label))]

p_forest_primary <- ggplot(
  f2,
  aes(
    x = HR,
    y = label,
    xmin = CI_low,
    xmax = CI_high
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_errorbarh(
    height = 0,
    linewidth = 0.55
  ) +
  geom_point(
    shape = 18,
    size = 2.4
  ) +
  geom_text(
    aes(
      x = 3.35,
      label = estimate_text
    ),
    hjust = 1,
    size = 3.0,
    family = "sans"
  ) +
  scale_x_log10(
    limits = c(0.55, 3.6),
    breaks = c(0.6, 0.75, 1, 1.5, 2, 3)
  ) +
  labs(
    x = "Hazard ratio",
    y = NULL
  ) +
  theme_manuscript(9.5) +
  theme(
    axis.text.y = element_text(hjust = 1),
    panel.grid.major.y = element_blank()
  )

save_figure(
  p_forest_primary,
  "Figure_2_Primary_robustness_forest",
  width = 7.2,
  height = 4.8
)

# ------------------------------------------------------------
# 16. FIGURE 3: EXPLORATORY SIGNAL FOREST
# ------------------------------------------------------------
f3 <- copy(t3_raw)

f3[
  grepl("^Sustained 1-y:", Analysis),
  group := "1-y sustained exposure"
]
f3[
  grepl("^Strict sustained 1-y:", Analysis),
  group := "Strict 1-y sustained exposure"
]
f3[
  grepl("^Persistent 3-y:", Analysis),
  group := "3-y persistent exposure"
]
f3[
  grepl("^US ", Analysis),
  group := "Finasteride vs dutasteride"
]

f3[, short_label := Analysis]
f3[, short_label := gsub(
  "^Sustained 1-y: ",
  "",
  short_label
)]
f3[, short_label := gsub(
  "^Strict sustained 1-y: ",
  "",
  short_label
)]
f3[, short_label := gsub(
  "^Persistent 3-y: ",
  "",
  short_label
)]
f3[, estimate_text := fmt_hr_ci(HR, CI_low, CI_high)]

group_order <- c(
  "1-y sustained exposure",
  "Strict 1-y sustained exposure",
  "3-y persistent exposure",
  "Finasteride vs dutasteride"
)

f3[, group := factor(group, levels = group_order)]

# Within each facet preserve manuscript order.
f3[, row_id := .I]
f3[, short_label := factor(short_label, levels = rev(unique(short_label)))]

p_forest_signal <- ggplot(
  f3,
  aes(
    x = HR,
    y = short_label,
    xmin = CI_low,
    xmax = CI_high
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_errorbarh(
    height = 0,
    linewidth = 0.5
  ) +
  geom_point(
    shape = 18,
    size = 2.3
  ) +
  geom_text(
    aes(
      x = 18.0,
      label = estimate_text
    ),
    hjust = 1,
    size = 2.8,
    family = "sans"
  ) +
  facet_grid(
    group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_x_log10(
    limits = c(0.35, 20.0),
    breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 5)
  ) +
  labs(
    x = "Hazard ratio",
    y = NULL
  ) +
  theme_manuscript(9.0) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(size = 8.5),
    panel.spacing.y = grid::unit(0.45, "lines")
  )

save_figure(
  p_forest_signal,
  "Figure_3_Exploratory_persistence_forest",
  width = 8.5,
  height = 7.8
)

# ------------------------------------------------------------
# 17. SUPPLEMENTARY FIGURE S1: LOVE PLOT
# ------------------------------------------------------------
love <- primary_balance[
  order(-abs_smd_unweighted)
][1:min(.N, 25)]

love[, Covariate := pretty_covariate(covariate)]

love_long <- rbind(
  love[
    ,
    .(
      Covariate,
      SMD = abs_smd_unweighted,
      Timing = "Before weighting"
    )
  ],
  love[
    ,
    .(
      Covariate,
      SMD = abs_smd_weighted,
      Timing = "After overlap weighting"
    )
  ]
)

love_long[, Timing := factor(
  Timing,
  levels = c("Before weighting", "After overlap weighting")
)]

love_order <- love[
  order(abs_smd_unweighted),
  Covariate
]
love_long[
  ,
  Covariate := factor(Covariate, levels = love_order)
]

p_love <- ggplot(
  love_long,
  aes(
    x = SMD,
    y = Covariate,
    shape = Timing
  )
) +
  geom_vline(
    xintercept = 0.10,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(size = 2.3) +
  scale_shape_manual(values = c("Before weighting" = 1, "After overlap weighting" = 16)) +
  scale_x_continuous(
    limits = c(0, max(0.45, max(love_long$SMD, na.rm = TRUE) * 1.05)),
    breaks = seq(0, 0.5, 0.1)
  ) +
  labs(
    x = "Absolute standardized mean difference",
    y = NULL
  ) +
  theme_manuscript(9.0) +
  theme(
    legend.position = "bottom"
  )

save_figure(
  p_love,
  "Figure_S1_Love_plot",
  width = 7.2,
  height = 6.8
)

# ------------------------------------------------------------
# 18. CAPTIONS / NOTES FILE
# ------------------------------------------------------------
captions <- c(
  "FIGURE LEGENDS",
  "",
  "Figure 1. Cohort Construction for the Primary Active-Comparator Analysis.",
  paste0(
    "The primary cohort comprised men with benign prostatic hyperplasia initiating ",
    "a 5alpha-reductase inhibitor (5-ARI) or uroselective alpha-blocker. ",
    "Outcome follow-up began at the 1-year landmark after treatment initiation."
  ),
  "",
  "Figure 2. Primary and Robustness Analyses of Incident Meningioma.",
  paste0(
    "Points indicate adjusted hazard ratios and horizontal lines indicate 95% confidence intervals. ",
    "The vertical dashed line indicates a hazard ratio of 1.0. ",
    "Estimates are overlap-weighted Cox proportional hazards models except the Fine-Gray ",
    "competing-risk sensitivity."
  ),
  "",
  "Figure 3. Exploratory Sustained-Exposure and Drug-Specific Analyses.",
  paste0(
    "Exploratory landmark analyses evaluated repeated or persistent medication documentation ",
    "and incorporated healthcare utilization during the exposure-assessment period into the ",
    "propensity-score model. The finasteride-dutasteride analyses were restricted to ",
    "contemporary US patients receiving exact benign-prostatic-hyperplasia doses."
  ),
  "",
  "Supplementary Figure 1. Covariate Balance Before and After Overlap Weighting.",
  paste0(
    "Shown are the 25 covariate terms with the largest absolute standardized mean differences ",
    "before weighting. The dashed line indicates an absolute standardized mean difference of 0.10."
  )
)

writeLines(
  captions,
  file.path(OUT_DIR, "MANUSCRIPT_FIGURE_LEGENDS.txt")
)

# ------------------------------------------------------------
# 19. OUTPUT INDEX
# ------------------------------------------------------------
index <- c(
  "MANUSCRIPT OUTPUTS",
  "",
  "TABLES",
  "Table_1_Baseline_characteristics.csv",
  "Table_2_Primary_robustness_results.csv",
  "Table_3_Exploratory_persistence_results.csv",
  "Table_S1_Covariate_balance.csv",
  "Table_S2_Latency_results.csv",
  "Table_S3_Rare_hard_endpoints.csv",
  "Manuscript_Tables.docx",
  "",
  "FIGURES (each exported as PNG, TIFF, and PDF)",
  "Figure_1_Cohort_flow",
  "Figure_2_Primary_robustness_forest",
  "Figure_3_Exploratory_persistence_forest",
  "Figure_S1_Love_plot",
  "",
  "CAPTIONS",
  "MANUSCRIPT_FIGURE_LEGENDS.txt"
)

writeLines(index, file.path(OUT_DIR, "OUTPUT_INDEX.txt"))

cat("\n============================================================\n")
cat("MANUSCRIPT OUTPUTS COMPLETE\n")
cat("============================================================\n")
cat("Output directory:\n", OUT_DIR, "\n\n", sep = "")
cat(paste(index, collapse = "\n"))
cat("\n")
