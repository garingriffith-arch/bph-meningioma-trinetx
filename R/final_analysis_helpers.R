# ============================================================
# final_analysis_helpers.R
# Shared propensity-score / overlap-weighted survival helpers
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

clean_cat <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unknown"
  factor(x)
}

n_observed_unique <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  x <- x[!is.na(x)]
  if (is.numeric(x)) x <- x[is.finite(x)]
  data.table::uniqueN(x)
}

prepare_common <- function(d) {
  d <- copy(as.data.table(d))

  if (!"trt" %in% names(d)) {
    stop("prepare_common(): data must contain binary trt (0/1).")
  }

  if ("index_date" %in% names(d)) {
    d[, index_date := as.Date(index_date)]
    d[, index_year := as.integer(format(index_date, "%Y"))]
  }

  # Standard categorical fields.
  cat_map <- c(
    race = "race_f",
    ethnicity = "ethnicity_f",
    marital_status = "marital_f",
    patient_regional_location = "region_f",
    patient_source_id = "source_f"
  )

  for (src in names(cat_map)) {
    dst <- cat_map[[src]]
    if (src %in% names(d)) {
      d[, (dst) := clean_cat(get(src))]
    }
  }

  # Stabilize skewed utilization variables when present.
  util_map <- c(
    n_encounters_prior_365 = "log_enc_365",
    n_encounters_prior_1095 = "log_enc_1095",
    n_ambulatory_prior_365 = "log_amb_365",
    n_ed_prior_365 = "log_ed_365",
    n_inpatient_prior_365 = "log_ip_365"
  )

  for (src in names(util_map)) {
    dst <- util_map[[src]]
    if (src %in% names(d)) {
      d[, (dst) := log1p(as.numeric(get(src)))]
    }
  }

  d
}

FACTOR_CANDIDATES <- c(
  "race_f",
  "ethnicity_f",
  "marital_f",
  "region_f",
  "source_f"
)

BINARY_CANDIDATES <- c(
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
  "bladder_stone",
  # permissive-model / optional covariates:
  "prior_prostate_cancer",
  "bph_near_index",
  "known_finasteride_1mg",
  "bmi_missing"
)

UTILIZATION_CANDIDATES <- c(
  "log_enc_365",
  "log_enc_1095",
  "log_amb_365",
  "log_ed_365",
  "log_ip_365"
)

build_ps_spec <- function(d, label, extra_numeric = character()) {

  d <- copy(d)

  # Drop unused levels after subgrouping.
  for (v in FACTOR_CANDIDATES) {
    if (v %in% names(d) && is.factor(d[[v]])) {
      d[[v]] <- droplevels(d[[v]])
    }
  }

  candidate_names <- c(
    "age_at_index",
    "index_year",
    FACTOR_CANDIDATES,
    BINARY_CANDIDATES,
    UTILIZATION_CANDIDATES,
    extra_numeric
  )

  diagnostics <- rbindlist(lapply(candidate_names, function(v) {
    nuniq <- if (v %in% names(d)) n_observed_unique(d[[v]]) else 0L
    data.table(
      analysis = label,
      variable = v,
      n_unique = as.integer(nuniq),
      action = ifelse(nuniq >= 2L, "KEEP", "DROP_ZERO_VARIANCE_OR_ABSENT")
    )
  }))

  keep_factor <- FACTOR_CANDIDATES[
    vapply(
      FACTOR_CANDIDATES,
      function(v) v %in% names(d) && n_observed_unique(d[[v]]) >= 2L,
      logical(1)
    )
  ]

  keep_binary <- BINARY_CANDIDATES[
    vapply(
      BINARY_CANDIDATES,
      function(v) v %in% names(d) && n_observed_unique(d[[v]]) >= 2L,
      logical(1)
    )
  ]

  keep_util <- UTILIZATION_CANDIDATES[
    vapply(
      UTILIZATION_CANDIDATES,
      function(v) v %in% names(d) && n_observed_unique(d[[v]]) >= 2L,
      logical(1)
    )
  ]

  keep_extra <- extra_numeric[
    vapply(
      extra_numeric,
      function(v) v %in% names(d) && n_observed_unique(d[[v]]) >= 2L,
      logical(1)
    )
  ]

  ps_terms <- character(0)
  bal_terms <- character(0)

  if ("age_at_index" %in% names(d)) {
    n_age <- n_observed_unique(d$age_at_index)
    if (n_age >= 6L) {
      ps_terms <- c(ps_terms, "splines::ns(age_at_index, df=5)")
      bal_terms <- c(bal_terms, "age_at_index")
    } else if (n_age >= 2L) {
      ps_terms <- c(ps_terms, "age_at_index")
      bal_terms <- c(bal_terms, "age_at_index")
    }
  }

  if ("index_year" %in% names(d)) {
    n_year <- n_observed_unique(d$index_year)
    if (n_year >= 6L) {
      ps_terms <- c(ps_terms, "splines::ns(index_year, df=5)")
      bal_terms <- c(bal_terms, "index_year")
    } else if (n_year >= 2L) {
      ps_terms <- c(ps_terms, "index_year")
      bal_terms <- c(bal_terms, "index_year")
    }
  }

  # Optional continuous variables use a flexible 4-df spline when possible.
  for (v in keep_extra) {
    if (n_observed_unique(d[[v]]) >= 5L) {
      ps_terms <- c(ps_terms, sprintf("splines::ns(%s, df=4)", v))
    } else {
      ps_terms <- c(ps_terms, v)
    }
    bal_terms <- c(bal_terms, v)
  }

  ps_terms <- c(ps_terms, keep_factor, keep_binary, keep_util)
  bal_terms <- c(bal_terms, keep_factor, keep_binary, keep_util)

  ps_terms <- unique(ps_terms)
  bal_terms <- unique(bal_terms)

  if (length(ps_terms) == 0L) {
    stop(label, ": no varying propensity covariates remain.")
  }

  list(
    data = d,
    ps_formula = as.formula(
      paste("trt ~", paste(ps_terms, collapse = " + "))
    ),
    balance_formula = as.formula(
      paste("~", paste(bal_terms, collapse = " + "), "- 1")
    ),
    diagnostics = diagnostics,
    dropped = diagnostics[action != "KEEP", variable]
  )
}

wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

wvar <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (sum(ok) < 2L) return(NA_real_)
  m <- wmean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])
}

smd_numeric <- function(x, z, w = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  m1 <- wmean(x[z == 1], w[z == 1])
  m0 <- wmean(x[z == 0], w[z == 0])
  v1 <- wvar(x[z == 1], w[z == 1])
  v0 <- wvar(x[z == 0], w[z == 0])
  denom <- sqrt((v1 + v0) / 2)
  if (!is.finite(denom) || denom == 0) return(0)
  (m1 - m0) / denom
}

make_balance_matrix <- function(d, balance_formula) {
  mm <- model.matrix(balance_formula, data = d)

  if (ncol(mm) > 0L) {
    keep <- apply(mm, 2, function(x) {
      x <- x[is.finite(x)]
      length(unique(x)) >= 2L
    })
    mm <- mm[, keep, drop = FALSE]
  }

  mm
}

balance_table <- function(d, balance_formula) {
  mm <- make_balance_matrix(d, balance_formula)
  z <- d$trt

  if (ncol(mm) == 0L) {
    return(data.table(
      covariate = character(),
      smd_unweighted = numeric(),
      smd_weighted = numeric(),
      abs_smd_unweighted = numeric(),
      abs_smd_weighted = numeric()
    ))
  }

  unweighted <- apply(
    mm, 2, smd_numeric, z = z, w = rep(1, nrow(d))
  )
  weighted <- apply(
    mm, 2, smd_numeric, z = z, w = d$ow
  )

  data.table(
    covariate = colnames(mm),
    smd_unweighted = as.numeric(unweighted),
    smd_weighted = as.numeric(weighted),
    abs_smd_unweighted = abs(as.numeric(unweighted)),
    abs_smd_weighted = abs(as.numeric(weighted))
  )[order(-abs_smd_weighted)]
}

ess <- function(w) {
  if (length(w) == 0L || sum(w^2) == 0) return(NA_real_)
  (sum(w)^2) / sum(w^2)
}

extract_hr <- function(fit, term = "trt") {
  b <- coef(fit)[term]
  se <- sqrt(diag(vcov(fit)))[term]
  data.table(
    HR = exp(b),
    CI_low = exp(b - 1.96 * se),
    CI_high = exp(b + 1.96 * se),
    p = 2 * pnorm(-abs(b / se))
  )
}

estimate_overlap_weights <- function(
    d,
    label,
    extra_numeric = character()
) {
  d <- prepare_common(d)

  if (!all(sort(unique(d$trt)) == c(0, 1))) {
    stop(label, ": both treatment groups are required.")
  }

  spec <- build_ps_spec(d, label, extra_numeric = extra_numeric)
  d <- spec$data

  mf <- model.frame(spec$ps_formula, data = d, na.action = na.pass)
  cc <- complete.cases(mf)

  if (any(!cc)) {
    warning(
      label, ": removing ", sum(!cc),
      " rows with missing PS-model covariates."
    )
    d <- d[cc]
  }

  ps_fit <- glm(
    spec$ps_formula,
    data = d,
    family = binomial(),
    control = glm.control(maxit = 75)
  )

  if (!isTRUE(ps_fit$converged)) {
    warning(label, ": propensity-score model did not converge.")
  }

  d[, ps := predict(ps_fit, newdata = d, type = "response")]
  d[, ps := pmin(pmax(ps, 1e-6), 1 - 1e-6)]
  d[, ow := ifelse(trt == 1, 1 - ps, ps)]

  bal <- balance_table(d, spec$balance_formula)
  bal[, analysis := label]

  list(
    data = d,
    ps_fit = ps_fit,
    balance = bal,
    diagnostics = spec$diagnostics,
    ps_formula = paste(deparse(spec$ps_formula), collapse = " ")
  )
}

run_ato_generic <- function(
    d,
    label,
    event_col,
    time_days_col,
    extra_numeric = character(),
    admin_horizon_years = NULL
) {
  d <- copy(as.data.table(d))

  if (!event_col %in% names(d)) {
    stop(label, ": missing event column ", event_col)
  }
  if (!time_days_col %in% names(d)) {
    stop(label, ": missing time column ", time_days_col)
  }

  d[, analysis_event := as.integer(get(event_col))]
  d[, analysis_time_y := as.numeric(get(time_days_col)) / 365.25]

  d <- d[
    is.finite(analysis_time_y) &
    analysis_time_y > 0 &
    !is.na(analysis_event)
  ]

  if (!is.null(admin_horizon_years)) {
    h <- as.numeric(admin_horizon_years)
    d[, analysis_event := as.integer(
      analysis_event == 1 & analysis_time_y <= h
    )]
    d[, analysis_time_y := pmin(analysis_time_y, h)]
  }

  wt <- estimate_overlap_weights(
    d,
    label = label,
    extra_numeric = extra_numeric
  )
  d <- wt$data

  crude <- coxph(
    Surv(analysis_time_y, analysis_event) ~ trt,
    data = d,
    ties = "efron"
  )

  weighted <- coxph(
    Surv(analysis_time_y, analysis_event) ~ trt,
    data = d,
    weights = ow,
    robust = TRUE,
    ties = "efron"
  )

  crude_hr <- extract_hr(crude)
  weighted_hr <- extract_hr(weighted)

  result <- data.table(
    analysis = label,
    n = nrow(d),
    n_treated = sum(d$trt == 1),
    n_control = sum(d$trt == 0),
    events = sum(d$analysis_event),
    events_treated = sum(d$analysis_event[d$trt == 1]),
    events_control = sum(d$analysis_event[d$trt == 0]),
    crude_HR = crude_hr$HR,
    crude_CI_low = crude_hr$CI_low,
    crude_CI_high = crude_hr$CI_high,
    crude_p = crude_hr$p,
    weighted_HR = weighted_hr$HR,
    weighted_CI_low = weighted_hr$CI_low,
    weighted_CI_high = weighted_hr$CI_high,
    weighted_p = weighted_hr$p,
    ESS_treated = ess(d$ow[d$trt == 1]),
    ESS_control = ess(d$ow[d$trt == 0]),
    ps_min = min(d$ps),
    ps_p01 = as.numeric(quantile(d$ps, 0.01)),
    ps_median = median(d$ps),
    ps_p99 = as.numeric(quantile(d$ps, 0.99)),
    ps_max = max(d$ps)
  )

  incidence <- rbindlist(lapply(c(0L, 1L), function(z) {
    x <- d[trt == z]
    data.table(
      analysis = label,
      trt = z,
      n = nrow(x),
      events = sum(x$analysis_event),
      person_years = sum(x$analysis_time_y),
      crude_rate_per_100k_py =
        100000 * sum(x$analysis_event) / sum(x$analysis_time_y),
      weighted_events = sum(x$ow * x$analysis_event),
      weighted_person_years = sum(x$ow * x$analysis_time_y),
      weighted_rate_per_100k_py =
        100000 * sum(x$ow * x$analysis_event) /
        sum(x$ow * x$analysis_time_y),
      weight_sum = sum(x$ow),
      ESS = ess(x$ow)
    )
  }))

  zph_out <- NULL
  zph <- try(cox.zph(weighted, transform = "km"), silent = TRUE)
  if (!inherits(zph, "try-error")) {
    ztab <- as.data.frame(zph$table)
    ztab$term <- rownames(ztab)
    zph_out <- as.data.table(ztab)
    zph_out[, analysis := label]
  }

  list(
    result = result,
    incidence = incidence,
    balance = wt$balance,
    diagnostics = wt$diagnostics,
    data = d,
    ps_fit = wt$ps_fit,
    weighted_fit = weighted,
    crude_fit = crude,
    zph = zph_out
  )
}

balance_summary <- function(balance_dt) {
  balance_dt[, .(
    max_abs_smd_unweighted = max(abs_smd_unweighted, na.rm = TRUE),
    max_abs_smd_weighted = max(abs_smd_weighted, na.rm = TRUE),
    n_weighted_smd_gt_0_10 = sum(abs_smd_weighted > 0.10, na.rm = TRUE),
    n_weighted_smd_gt_0_05 = sum(abs_smd_weighted > 0.05, na.rm = TRUE)
  ), by = analysis]
}
