# ============================================================
# 04_primary_analysis_v2.R
# Finasteride / Dutasteride and Incident Meningioma Study
#
# PATCHED: dynamically drops zero-variance / one-level covariates before
# each propensity-score model, preventing factor-contrast failures.
#
# Primary causal analysis:
#   5ARI initiation vs uroselective alpha-blocker initiation
#
# Estimand:
#   Average Treatment Effect in the Overlap population (ATO)
#
# Method:
#   logistic propensity score -> overlap weights ->
#   weighted Cox proportional hazards model with robust SE
#
# Also runs:
#   - unweighted crude Cox
#   - 5-year administratively censored weighted Cox
#   - >=2-code confirmed meningioma sensitivity
#   - finasteride vs alpha blocker
#   - dutasteride vs alpha blocker
#   - finasteride 5 mg vs alpha blocker
#   - 2-year lag sensitivity
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)

required_pkgs <- c("DBI", "duckdb", "data.table", "survival", "splines")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)

library(DBI)
library(duckdb)
library(data.table)
library(survival)
library(splines)

WORK_DIR <- "C:/TriNetX/bph_m_work"
DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")
OUT_DIR <- file.path(WORK_DIR, "analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_THREADS <- min(24L, max(4L, parallel::detectCores(logical = TRUE) - 6L))
MEMORY_LIMIT <- "48GB"

con <- dbConnect(duckdb::duckdb(), dbdir = DB_PATH, read_only = FALSE)
on.exit({
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
}, add = TRUE)

dbExecute(con, sprintf("PRAGMA threads=%d;", N_THREADS))
dbExecute(con, sprintf("PRAGMA memory_limit='%s';", MEMORY_LIMIT))

needed <- c("analysis_cohort_primary_corrected", "analysis_cohort_lag2")
missing <- setdiff(needed, dbListTables(con))
if (length(missing) > 0) {
  stop("Run 03_refresh_meningioma_outcomes.R first. Missing: ",
       paste(missing, collapse=", "))
}

# ------------------------------------------------------------
# 1. LOAD COMPACT ANALYSIS COHORT
# ------------------------------------------------------------
cols <- c(
  "patient_id", "exposure", "index_drug",
  "index_finasteride_strength", "index_dutasteride_strength",
  "index_date", "age_at_index",
  "race", "ethnicity", "marital_status",
  "patient_regional_location", "patient_source_id",
  "hypertension", "diabetes", "obesity", "ckd", "smoking",
  "cad", "heart_failure", "cerebrovascular_disease",
  "hyperlipidemia", "osteoporosis", "alopecia",
  "urinary_retention", "uti", "hematuria", "bladder_stone",
  "n_encounters_prior_365", "n_encounters_prior_1095",
  "n_ambulatory_prior_365", "n_ed_prior_365", "n_inpatient_prior_365",
  "meningioma_event_strict", "meningioma_event_confirmed",
  "followup_days_strict"
)

select_sql <- paste(cols, collapse=", ")

dat <- as.data.table(dbGetQuery(
  con,
  paste0("SELECT ", select_sql, " FROM analysis_cohort_primary_corrected;")
))

lag2 <- as.data.table(dbGetQuery(
  con,
  paste0("SELECT ", select_sql, " FROM analysis_cohort_lag2;")
))

dbDisconnect(con, shutdown = TRUE)

# ------------------------------------------------------------
# 2. CLEAN VARIABLES
# ------------------------------------------------------------
clean_cat <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unknown"
  factor(x)
}

prep_data <- function(d) {
  d <- copy(d)

  d[, trt := as.integer(exposure == "5ARI")]
  d[, index_year := as.integer(format(as.Date(index_date), "%Y"))]
  d[, time_y := as.numeric(followup_days_strict) / 365.25]
  d[, event := as.integer(meningioma_event_strict)]
  d[, event_confirmed := as.integer(meningioma_event_confirmed)]

  d[, race_f := clean_cat(race)]
  d[, ethnicity_f := clean_cat(ethnicity)]
  d[, marital_f := clean_cat(marital_status)]
  d[, region_f := clean_cat(patient_regional_location)]
  d[, source_f := clean_cat(patient_source_id)]

  # Stabilize highly skewed utilization variables.
  d[, log_enc_365 := log1p(as.numeric(n_encounters_prior_365))]
  d[, log_enc_1095 := log1p(as.numeric(n_encounters_prior_1095))]
  d[, log_amb_365 := log1p(as.numeric(n_ambulatory_prior_365))]
  d[, log_ed_365 := log1p(as.numeric(n_ed_prior_365))]
  d[, log_ip_365 := log1p(as.numeric(n_inpatient_prior_365))]

  d
}

dat <- prep_data(dat)
lag2 <- prep_data(lag2)

# ------------------------------------------------------------
# 3. DYNAMIC PROPENSITY-SCORE SPECIFICATION
# ------------------------------------------------------------
# The prior version used a fixed formula containing all categorical variables.
# If any factor has only one observed value in an analytic cohort (for example,
# patient_source_id may be uniformly EHR), glm/model.matrix fails with:
#
#   contrasts can be applied only to factors with 2 or more levels
#
# We therefore evaluate covariate variation separately within EVERY analytic
# cohort/subgroup and drop only variables that contain no information there.

factor_candidates <- c(
  "race_f",
  "ethnicity_f",
  "marital_f",
  "region_f",
  "source_f"
)

binary_numeric_candidates <- c(
  "hypertension",
  "diabetes",
  "obesity",
  "ckd",
  "smoking",
  "cad",
  "heart_failure",
  "cerebrovascular_disease",
  "hyperlipidemia",
  "osteoporosis",
  "alopecia",
  "urinary_retention",
  "uti",
  "hematuria",
  "bladder_stone"
)

utilization_candidates <- c(
  "log_enc_365",
  "log_enc_1095",
  "log_amb_365",
  "log_ed_365",
  "log_ip_365"
)

n_observed_unique <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- x[!is.na(x)]
  if (is.numeric(x)) x <- x[is.finite(x)]
  data.table::uniqueN(x)
}

build_ps_spec <- function(d, label) {

  d <- copy(d)

  # Drop unused factor levels after each drug-specific subset.
  for (v in factor_candidates) {
    if (v %in% names(d) && is.factor(d[[v]])) {
      d[[v]] <- droplevels(d[[v]])
    }
  }

  all_candidates <- c(
    "age_at_index",
    "index_year",
    factor_candidates,
    binary_numeric_candidates,
    utilization_candidates
  )

  diagnostics <- rbindlist(lapply(all_candidates, function(v) {
    nuniq <- if (v %in% names(d)) n_observed_unique(d[[v]]) else 0L

    data.table(
      analysis = label,
      variable = v,
      n_unique = as.integer(nuniq),
      action = ifelse(nuniq >= 2L, "KEEP", "DROP_ZERO_VARIANCE")
    )
  }))

  keep_factor <- factor_candidates[
    vapply(
      factor_candidates,
      function(v) v %in% names(d) && n_observed_unique(d[[v]]) >= 2L,
      logical(1)
    )
  ]

  keep_binary <- binary_numeric_candidates[
    vapply(
      binary_numeric_candidates,
      function(v) v %in% names(d) && n_observed_unique(d[[v]]) >= 2L,
      logical(1)
    )
  ]

  keep_util <- utilization_candidates[
    vapply(
      utilization_candidates,
      function(v) v %in% names(d) && n_observed_unique(d[[v]]) >= 2L,
      logical(1)
    )
  ]

  ps_terms <- character(0)
  balance_terms <- character(0)

  # Age: spline when sufficiently variable, otherwise linear if possible.
  n_age <- n_observed_unique(d$age_at_index)
  if (n_age >= 6L) {
    ps_terms <- c(ps_terms, "splines::ns(age_at_index, df=5)")
    balance_terms <- c(balance_terms, "age_at_index")
  } else if (n_age >= 2L) {
    ps_terms <- c(ps_terms, "age_at_index")
    balance_terms <- c(balance_terms, "age_at_index")
  }

  # Calendar year: same logic.
  n_year <- n_observed_unique(d$index_year)
  if (n_year >= 6L) {
    ps_terms <- c(ps_terms, "splines::ns(index_year, df=5)")
    balance_terms <- c(balance_terms, "index_year")
  } else if (n_year >= 2L) {
    ps_terms <- c(ps_terms, "index_year")
    balance_terms <- c(balance_terms, "index_year")
  }

  ps_terms <- c(
    ps_terms,
    keep_factor,
    keep_binary,
    keep_util
  )

  balance_terms <- c(
    balance_terms,
    keep_factor,
    keep_binary,
    keep_util
  )

  if (length(ps_terms) == 0L) {
    stop(label, ": no varying propensity-score covariates remain.")
  }

  ps_formula <- as.formula(
    paste("trt ~", paste(ps_terms, collapse=" + "))
  )

  balance_formula <- as.formula(
    paste("~", paste(balance_terms, collapse=" + "), "- 1")
  )

  list(
    data = d,
    ps_formula = ps_formula,
    balance_formula = balance_formula,
    diagnostics = diagnostics,
    kept = ps_terms,
    dropped = diagnostics[action != "KEEP", variable]
  )
}

# ------------------------------------------------------------
# 4. WEIGHTED BALANCE HELPERS
# ------------------------------------------------------------
wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

wvar <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (sum(ok) < 2) return(NA_real_)
  m <- wmean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])
}

smd_numeric <- function(x, z, w=NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  m1 <- wmean(x[z==1], w[z==1])
  m0 <- wmean(x[z==0], w[z==0])
  v1 <- wvar(x[z==1], w[z==1])
  v0 <- wvar(x[z==0], w[z==0])
  denom <- sqrt((v1 + v0) / 2)
  if (!is.finite(denom) || denom == 0) return(0)
  (m1 - m0) / denom
}

make_balance_matrix <- function(d, balance_formula) {
  mm <- model.matrix(balance_formula, data=d)

  # Defensive removal of any columns that became constant after expansion.
  if (ncol(mm) > 0) {
    keep <- apply(mm, 2, function(x) {
      x <- x[is.finite(x)]
      length(unique(x)) >= 2L
    })
    mm <- mm[, keep, drop=FALSE]
  }

  mm
}

balance_table <- function(d, balance_formula) {
  mm <- make_balance_matrix(d, balance_formula)
  z <- d$trt

  if (ncol(mm) == 0L) {
    return(data.table(
      covariate=character(),
      smd_unweighted=numeric(),
      smd_weighted=numeric(),
      abs_smd_unweighted=numeric(),
      abs_smd_weighted=numeric()
    ))
  }

  unweighted <- apply(mm, 2, smd_numeric, z=z, w=rep(1, nrow(d)))
  weighted <- apply(mm, 2, smd_numeric, z=z, w=d$ow)

  data.table(
    covariate=colnames(mm),
    smd_unweighted=as.numeric(unweighted),
    smd_weighted=as.numeric(weighted),
    abs_smd_unweighted=abs(as.numeric(unweighted)),
    abs_smd_weighted=abs(as.numeric(weighted))
  )[order(-abs_smd_weighted)]
}

# ------------------------------------------------------------
# 5. MODEL RESULT HELPERS
# ------------------------------------------------------------
extract_hr <- function(fit, term="trt") {
  b <- coef(fit)[term]
  se <- sqrt(diag(vcov(fit)))[term]
  data.table(
    HR=exp(b),
    CI_low=exp(b - 1.96*se),
    CI_high=exp(b + 1.96*se),
    p=2*pnorm(-abs(b/se))
  )
}

ess <- function(w) {
  (sum(w)^2) / sum(w^2)
}

incidence_table <- function(d) {
  rbindlist(lapply(c(0L,1L), function(z) {
    x <- d[trt==z]
    data.table(
      exposure=ifelse(z==1, "5ARI", "URO_ALPHA"),
      n=nrow(x),
      events=sum(x$event),
      person_years=sum(x$time_y),
      crude_rate_per_100k_py=
        100000 * sum(x$event) / sum(x$time_y),
      weighted_events=sum(x$ow*x$event),
      weighted_person_years=sum(x$ow*x$time_y),
      weighted_rate_per_100k_py=
        100000 * sum(x$ow*x$event) / sum(x$ow*x$time_y),
      weight_sum=sum(x$ow),
      ESS=ess(x$ow)
    )
  }))
}

# ------------------------------------------------------------
# 6. CORE ANALYSIS FUNCTION
# ------------------------------------------------------------
run_ato <- function(d, label, confirmed=FALSE, admin_5y=FALSE) {

  d <- copy(d)

  if (length(unique(d$trt)) != 2) {
    stop(label, ": both treatment groups are required.")
  }

  # Build a cohort-specific PS formula, automatically dropping
  # one-level factors and other zero-variance variables.
  spec <- build_ps_spec(d, label)
  d <- spec$data

  if (length(spec$dropped) > 0) {
    cat("\n", label, "\n", sep="")
    cat("Dropped zero-variance propensity covariates: ",
        paste(spec$dropped, collapse=", "), "\n", sep="")
  }

  # Check for missingness in the actual model frame. We do not silently
  # change the estimand by allowing glm() to remove rows without reporting it.
  mf <- model.frame(spec$ps_formula, data=d, na.action=na.pass)
  cc <- complete.cases(mf)

  if (any(!cc)) {
    warning(
      label, ": removing ", sum(!cc),
      " rows with missing PS-model covariates."
    )
    d <- d[cc]
  }

  # Propensity score model
  ps_fit <- glm(
    spec$ps_formula,
    data=d,
    family=binomial(),
    control=glm.control(maxit=75)
  )

  if (!isTRUE(ps_fit$converged)) {
    warning(label, ": propensity-score logistic regression did not converge.")
  }

  d[, ps := predict(ps_fit, newdata=d, type="response")]
  d[, ps := pmin(pmax(ps, 1e-6), 1-1e-6)]

  # ATO overlap weights:
  # treated -> 1 - PS
  # control -> PS
  d[, ow := ifelse(trt==1, 1-ps, ps)]

  # Outcome selection
  if (confirmed) {
    d[, analysis_event := event_confirmed]
  } else {
    d[, analysis_event := event]
  }

  d[, analysis_time := time_y]

  if (admin_5y) {
    d[, analysis_event := as.integer(analysis_event==1 & analysis_time <= 5)]
    d[, analysis_time := pmin(analysis_time, 5)]
  }

  # Crude model
  crude <- coxph(
    Surv(analysis_time, analysis_event) ~ trt,
    data=d,
    ties="efron"
  )

  # Primary weighted marginal Cox model
  weighted <- coxph(
    Surv(analysis_time, analysis_event) ~ trt,
    data=d,
    weights=ow,
    robust=TRUE,
    ties="efron"
  )

  crude_hr <- extract_hr(crude)
  weighted_hr <- extract_hr(weighted)

  result <- data.table(
    analysis=label,
    n=nrow(d),
    n_5ari=sum(d$trt==1),
    n_uro_alpha=sum(d$trt==0),
    events=sum(d$analysis_event),
    events_5ari=sum(d$analysis_event[d$trt==1]),
    events_uro_alpha=sum(d$analysis_event[d$trt==0]),
    crude_HR=crude_hr$HR,
    crude_CI_low=crude_hr$CI_low,
    crude_CI_high=crude_hr$CI_high,
    crude_p=crude_hr$p,
    weighted_HR=weighted_hr$HR,
    weighted_CI_low=weighted_hr$CI_low,
    weighted_CI_high=weighted_hr$CI_high,
    weighted_p=weighted_hr$p,
    ESS_5ari=ess(d$ow[d$trt==1]),
    ESS_uro_alpha=ess(d$ow[d$trt==0]),
    ps_min=min(d$ps),
    ps_p01=quantile(d$ps, .01),
    ps_median=median(d$ps),
    ps_p99=quantile(d$ps, .99),
    ps_max=max(d$ps)
  )

  bal <- balance_table(d, spec$balance_formula)
  bal[, analysis := label]

  # PH diagnostic for treatment effect
  zph <- try(cox.zph(weighted, transform="km"), silent=TRUE)
  zph_out <- NULL
  if (!inherits(zph, "try-error")) {
    ztab <- as.data.frame(zph$table)
    ztab$term <- rownames(ztab)
    zph_out <- as.data.table(ztab)
    zph_out[, analysis := label]
  }

  list(
    result=result,
    balance=bal,
    data=d,
    ps_fit=ps_fit,
    crude_fit=crude,
    weighted_fit=weighted,
    zph=zph_out,
    predictor_diagnostics=spec$diagnostics,
    ps_formula=deparse(spec$ps_formula)
  )
}

# ------------------------------------------------------------
# 7. PRIMARY ANALYSIS
# ------------------------------------------------------------
primary <- run_ato(
  dat,
  label="Primary: 1-y lag, strict meningioma, all 5ARI vs uro-alpha"
)

# ------------------------------------------------------------
# 8. PRESPECIFIED SENSITIVITIES
# ------------------------------------------------------------

# 5-year administrative horizon
primary_5y <- run_ato(
  dat,
  label="Sensitivity: 5-y administrative horizon",
  admin_5y=TRUE
)

# More specific outcome: >=2 strict meningioma diagnosis dates >=30 days apart
confirmed <- run_ato(
  dat,
  label="Sensitivity: confirmed meningioma (>=2 codes >=30 d apart)",
  confirmed=TRUE
)

# 2-year lag
lag2_res <- run_ato(
  lag2,
  label="Sensitivity: 2-y lag"
)

# Finasteride vs uro-alpha
fin <- dat[
  exposure=="URO_ALPHA" |
  (exposure=="5ARI" & index_drug=="finasteride")
]
fin[, trt := as.integer(exposure=="5ARI")]

fin_res <- run_ato(
  fin,
  label="Drug-specific: finasteride vs uro-alpha"
)

# Dutasteride vs uro-alpha
dut <- dat[
  exposure=="URO_ALPHA" |
  (exposure=="5ARI" & index_drug=="dutasteride")
]
dut[, trt := as.integer(exposure=="5ARI")]

dut_res <- run_ato(
  dut,
  label="Drug-specific: dutasteride vs uro-alpha"
)

# Exact BPH-dose finasteride sensitivity
fin5 <- dat[
  exposure=="URO_ALPHA" |
  (
    exposure=="5ARI" &
    index_drug=="finasteride" &
    !is.na(index_finasteride_strength) &
    toupper(trimws(index_finasteride_strength))=="5 MG"
  )
]
fin5[, trt := as.integer(exposure=="5ARI")]

fin5_res <- run_ato(
  fin5,
  label="Dose-specific: finasteride 5 mg vs uro-alpha"
)

# ------------------------------------------------------------
# 9. COMBINE / SAVE RESULTS
# ------------------------------------------------------------
results <- rbindlist(list(
  primary$result,
  primary_5y$result,
  confirmed$result,
  lag2_res$result,
  fin_res$result,
  dut_res$result,
  fin5_res$result
), fill=TRUE)

fwrite(results, file.path(OUT_DIR, "04_model_results.csv"))

balances <- rbindlist(list(
  primary$balance,
  primary_5y$balance,
  confirmed$balance,
  lag2_res$balance,
  fin_res$balance,
  dut_res$balance,
  fin5_res$balance
), fill=TRUE)

fwrite(balances, file.path(OUT_DIR, "04_balance_all.csv"))


predictor_diagnostics <- rbindlist(list(
  primary$predictor_diagnostics,
  primary_5y$predictor_diagnostics,
  confirmed$predictor_diagnostics,
  lag2_res$predictor_diagnostics,
  fin_res$predictor_diagnostics,
  dut_res$predictor_diagnostics,
  fin5_res$predictor_diagnostics
), fill=TRUE)

fwrite(
  predictor_diagnostics,
  file.path(OUT_DIR, "04_predictor_diagnostics.csv")
)

# Primary balance summary
balance_summary <- balances[, .(
  max_abs_smd_unweighted=max(abs_smd_unweighted, na.rm=TRUE),
  max_abs_smd_weighted=max(abs_smd_weighted, na.rm=TRUE),
  n_weighted_smd_gt_0_10=sum(abs_smd_weighted > 0.10, na.rm=TRUE),
  n_weighted_smd_gt_0_05=sum(abs_smd_weighted > 0.05, na.rm=TRUE)
), by=analysis]

fwrite(balance_summary, file.path(OUT_DIR, "04_balance_summary.csv"))

# Propensity score distributions in primary analysis
ps_qc <- primary$data[, .(
  n=.N,
  ps_min=min(ps),
  ps_p01=quantile(ps,.01),
  ps_p05=quantile(ps,.05),
  ps_p25=quantile(ps,.25),
  ps_median=median(ps),
  ps_p75=quantile(ps,.75),
  ps_p95=quantile(ps,.95),
  ps_p99=quantile(ps,.99),
  ps_max=max(ps),
  weight_min=min(ow),
  weight_median=median(ow),
  weight_max=max(ow),
  ESS=ess(ow)
), by=exposure]

fwrite(ps_qc, file.path(OUT_DIR, "04_primary_ps_qc.csv"))

# Incidence table for primary weighted population
inc <- incidence_table(primary$data)
fwrite(inc, file.path(OUT_DIR, "04_primary_incidence.csv"))

# PH diagnostics
zph_all <- rbindlist(Filter(Negate(is.null), list(
  primary$zph,
  primary_5y$zph,
  confirmed$zph,
  lag2_res$zph,
  fin_res$zph,
  dut_res$zph,
  fin5_res$zph
)), fill=TRUE)

if (nrow(zph_all) > 0) {
  fwrite(zph_all, file.path(OUT_DIR, "04_ph_diagnostics.csv"))
}

# ------------------------------------------------------------
# 10. WEIGHTED KM RISK AT FIXED TIME POINTS
# ------------------------------------------------------------
km <- survfit(
  Surv(time_y, event) ~ exposure,
  data=primary$data,
  weights=ow
)

times_requested <- c(3,5,10)
times_use <- times_requested[times_requested <= max(primary$data$time_y)]

if (length(times_use) > 0) {
  sm <- summary(km, times=times_use, extend=FALSE)

  km_tab <- data.table(
    stratum=as.character(sm$strata),
    time=sm$time,
    survival=sm$surv,
    cumulative_risk=1-sm$surv,
    std_error=sm$std.err,
    lower_survival=sm$lower,
    upper_survival=sm$upper
  )

  fwrite(km_tab, file.path(OUT_DIR, "04_primary_weighted_KM_risk.csv"))
}

# ------------------------------------------------------------
# 11. PRINT KEY OUTPUT
# ------------------------------------------------------------
cat("\n================ DROPPED PS COVARIATES =============\n")
print(predictor_diagnostics[action != "KEEP"])

cat("\n================ MODEL RESULTS =====================\n")
print(results)

cat("\n================ BALANCE SUMMARY ===================\n")
print(balance_summary)

cat("\n================ PRIMARY PS QC =====================\n")
print(ps_qc)

cat("\n================ PRIMARY INCIDENCE =================\n")
print(inc)

cat("\nAnalysis complete.\n")
cat("Outputs:", OUT_DIR, "\n")
