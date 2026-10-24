# =============================================================================
# Stage 20. RECORD-PE + TARGET reporting audit and descriptive supplements
# Study: 5-alpha-reductase inhibitors and incident meningioma in BPH
# =============================================================================
# PURPOSE
#   Generate only reporting/descriptive analyses and editable checklists.
#   NO inferential models are refit; NO source tables are written or replaced.
#   The canonical final primary cohort is Stage 16's boundary-fixed cohort.
#   The target-trial-aligned day-zero analysis is a SENSITIVITY analysis.
#   A checklist worksheet is NOT proof of reporting-guideline compliance.
#
# INPUT: local DuckDB (read only), Stage 15 and 16 aggregate outputs,
#        existing R scripts (if available).
# OUTPUT: C:/TriNetX/bph_m_work/reporting_RECORD_PE_TARGET/
# REQUIRED packages: DBI, duckdb, data.table (already used in the project)
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)
R_DIR <- "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R"
WORK_DIR <- "C:/TriNetX/bph_m_work"
DB_PATH <- file.path(WORK_DIR, "duckdb", "bph_m.duckdb")
OUT <- file.path(WORK_DIR, "reporting_RECORD_PE_TARGET")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

required <- c("DBI", "duckdb", "data.table")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install missing packages before Stage 20: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(data.table)
})

log_lines <- character()
add_log <- function(level, item, explanation) {
  log_lines <<- c(log_lines, paste(level, item, explanation, sep = " | "))
  cat(sprintf("[%s] %s: %s\n", level, item, explanation))
}
write_csv <- function(x, filename) {
  data.table::fwrite(as.data.table(x), file.path(OUT, filename), na = "")
}
write_lines <- function(x, filename) {
  writeLines(enc2utf8(x), file.path(OUT, filename), useBytes = TRUE)
}
existing <- function(...) {
  p <- c(...)
  h <- p[file.exists(p)]
  if (length(h)) h[1] else NA_character_
}
maybe_csv <- function(path) {
  if (!is.na(path) && file.exists(path)) data.table::fread(path) else NULL
}

# These reporting labels are condensed summaries, not reproductions of the
# copyrighted TARGET checklist. Consult the official checklists when approving.
# RECORD-PE is an extension: the parent STROBE + RECORD items ALSO apply.
make_entries <- function(guideline, ids, topics, where, actions) {
  stopifnot(length(ids) == length(topics), length(ids) == length(where),
            length(ids) == length(actions))
  data.table(guideline = guideline, item = ids, topic = topics,
             proposed_location = where, evidence_or_action = actions,
             manuscript_page = "", supplement_page = "",
             author_assessment = "NOT_YET_VERIFIED", reviewer_notes = "")
}

strobe <- make_entries(
  "STROBE (base for RECORD-PE)", as.character(1:22),
  c("Title/abstract", "Background", "Objectives", "Study design", "Setting and dates",
    "Eligibility/selection", "Variables", "Measurement", "Bias", "Study size",
    "Quantitative variables", "Statistical methods/missingness/censoring", "Participant flow",
    "Baseline/missingness/follow-up", "Outcome events", "Adjusted and absolute results",
    "Other analyses", "Key findings", "Limitations", "Interpretation",
    "Generalisability", "Funding"),
  c("Title/Abstract", "Introduction", "Introduction", "Methods/Design", "Methods/Data source",
    "Methods/Population + Fig 1", "Methods/Outcomes + Supplementary Methods",
    "Methods/Covariates", "Methods/Robustness", "Methods/Study size", "Methods/Model",
    "Methods/Statistical analysis", "Fig 1", "Table 1 + Table S8 + Missingness",
    "Table 2 + Table S8", "Table 2/Results", "Tables 2-3/Supplement",
    "Discussion", "Discussion", "Discussion", "Discussion", "Declarations"),
  c("Mention design and TriNetX", "Set rationale", "Differentiate primary and exploratory hypotheses",
    "State primary landmark and day-zero sensitivity distinction", "Give study start/end and extract vintage",
    "Confirm entry order, once-only enrollment, washout", "Publish codes and phenotype algorithms",
    "Clarify variable provenance and comparison", "Explain bias mitigations and unresolved sources",
    "All eligible data, not sample size chosen by power", "Explain spline and covariate coding",
    "Specify model, complete missingness method, censoring, competing mortality",
    "Check all counts against final primary cohort", "Add group follow-up, missingness and PY",
    "Report meningiomas by treatment", "Report crude and weighted estimates, do not equate rates with risk",
    "Report all exploratory and sensitivity analyses", "Summarize null primary finding",
    "Discuss exposure misclassification, selection, residual confounding",
    "Avoid causal statements unsupported by design", "Overlap-population and health-network generalisability",
    "Declare no external study funding and funder role"))

record <- make_entries(
  "RECORD extension", c("1.1","1.2","1.3","6.1","6.2","6.3","7.1","12.1",
                        "12.2","12.3","13.1","19.1","22.1"),
  c("Routinely collected data type","Geography/period","Database linkage",
    "Population selection algorithms","Validation evidence","Linkage flow",
    "All codes/algorithms","Database access extent","Data cleaning",
    "Linkage methods","Filtering flow","EHR-specific limitations",
    "Protocol/data/code availability"),
  c("Abstract","Abstract/Methods","Methods","Supplementary Methods",
    "Supplementary Methods","Methods/Fig 1","Codebook supplement","Methods",
    "Supplementary Methods","Methods","Figure 1","Discussion","Declarations"),
  c("Say deidentified TriNetX EHR", "Verify US-only meaning and data dates",
    "State whether different databases were linked; N/A only if confirmed",
    "Source Stages 00-06 recovered; corroborate cohort data dictionary and chronology",
    "Distinguish internal corroboration from external validation",
    "N/A if no linkage; confirm", "Disclose med/diagnosis/outcome code lists",
    "Describe access to source extract, permission and restrictions",
    "Document joins, dates, duplicates, episode logic", "Specify if any linkage occurred",
    "Report each exclusion step and denominator", "Discuss encounter capture and coding",
    "Give code/protocol availability and controlled-data limitations"))

pe <- make_entries(
  "RECORD-PE extension",
  c("4.a","4.b","6.1.a","7.1.a","7.1.b","7.1.c","7.1.d","7.1.e",
    "7.1.f","7.1.g","8.a","12.1.a","12.1.b","19.1.a","20.a"),
  c("Multiple study designs", "Study design timeline", "Eligibility order and re-entry",
    "Medication code identification", "Medication recording source", "Exposure windows",
    "Outcome attribution to exposure", "Duration/dose approach", "Comparator choice",
    "Cross-class switching", "Prescribing environment", "Assumptions/diagnostics",
    "Justify distinct analyses", "Drug capture completeness", "Healthy-adherer/indication bias"),
  c("Methods/Design", "Figure S2", "Methods/Population", "Supplementary Methods",
    "Methods/Data source", "Methods/Exposure", "Methods/Outcome",
    "Methods/Exposure", "Methods/Population", "Methods/Exposure",
    "Methods/Data source", "Supplementary Methods", "Methods/Robustness",
    "Discussion", "Discussion"),
  c("Separate 1-y primary and day-0 sensitivity", "Schematic produced; verify timing against source code",
    "Stage 02/05 source recovered; verify once-only eligibility and captured history",
    "Verify original medication dictionary and strength exclusions",
    "Differentiate orders from fills/dispensations", "Exact 1-y and 3-y recorded-medication rules recovered in Stage 07",
    "Specify initiator-based ITT-like vs persistent landmark associations",
    "Explain why fills/days supply are unavailable", "Justify uroselective alpha-blocker",
    "Report dual initiation, switching and combination management",
    "Describe how participating providers record medications",
    "Report PS balance, overlap, IPCW, PH assessments and limits",
    "Explain distinct cohorts and their causal questions",
    "State EHR medication documentation may be incomplete",
    "Discuss selection into persistent subgroup and positive negative control"))

target <- make_entries(
  "TARGET (day-0 sensitivity only)",
  c("1a","1b","1c","2","3","4","5", paste0("6",letters[1:8]),
    paste0("7",letters[1:8]),"8","9","10","11","12","13","14","15","16",
    "17","18","19","20","21"),
  c("Abstract design", "Abstract data source", "Abstract summary", "Background",
    "Causal question", "Reason for emulation", "Data source", "Eligibility specification",
    "Treatment strategies specification", "Assignment specification", "Start/end follow-up",
    "Outcome specification", "Causal contrasts", "Identification assumptions",
    "Analysis plan", "Observed eligibility criteria", "Observed strategies", "Observed assignment",
    "Observed follow-up", "Observed outcomes", "Observed contrasts", "Observational assumptions",
    "Observed analysis", "Flow", "Baseline by group", "Follow-up by group",
    "Missingness by group", "Events by group", "Effect estimates", "Sensitivity analyses",
    "Interpretation", "Differences and assumptions", "Ethics determination",
    "Protocol registration", "Materials sharing", "Funding", "Disclosures"),
  c(rep("Abstract",3),rep("Introduction",3),"Methods/Data source",
    rep("Table S1 TARGET specification",16),"Figure 1", "Table 1",
    "Table S8", "Missingness supplement", "Table 2", "Table 2",
    "Table 2 and supplements", "Discussion", "Discussion",
    "Declarations", "Declarations", "Declarations", "Declarations", "Declarations"),
  c("Label day-zero component as target-trial-aligned sensitivity, not entire paper",
    "Say observational EHR network", "Give concise assumptions + estimate",
    "Clinical uncertainty", "Population, treatment, outcome, follow-up",
    "Why observational comparisons informative", "Network geography, vintage, record purpose",
    "Specify eligibility for hypothetical trial", "Define planned initiator strategies",
    "Hypothetical random assignment at day 0", "Set time 0 and endpoints",
    "Meningioma phenotype", "Define initiation-strategy effect",
    "State conditional exchangeability, positivity, consistency and censoring",
    "Specify causal effect estimator and missingness strategy",
    "Operationalize BPH, age, diagnosis exclusion", "Medication classification from raw records",
    "Day-0 classification without future treatment info", "Follow-up starts at day 0; verify censoring",
    "ICD codes and ascertainment", "ITT-like observational contrast; overlap population",
    "State confounding and record-availability assumptions", "Overlap-weighted Cox, IPCW sensitivities",
    "Day-zero selection and events", "Observed covariates by day-zero arm",
    "Generate day-zero descriptive table if cohort available", "Generate day-zero missingness",
    "Provide events per arm", "Adjusted HR and CI; absolute measures when defensible",
    "List variations and purpose", "Cautious interpretation", "No randomization; measurement limits",
    "RFD: non-human-subjects determination", "State if no registry/protocol registration",
    "Document controlled data and code sharing", "No funding", "Declarations"))

checklist <- rbindlist(list(strobe, record, pe, target), fill = TRUE)
write_csv(checklist, "RECORD_PE_TARGET_item_crosswalk_AUTHOR_REVIEW.csv")
stopifnot(nrow(pe) == 15L, nrow(record) == 13L, nrow(strobe) == 22L,
          length(unique(target$item)) == nrow(target))
add_log("CREATED", "Guideline crosswalk", paste(nrow(checklist), "checklist subitems; ALL require human sign-off and page numbers"))

# Source-code provenance; never infer exposure definitions from filenames alone.
need_stages <- 0:19
code_inventory <- data.table(stage = need_stages,
  expected_pattern = sprintf("%02d_*.R", need_stages),
  present = vapply(need_stages, function(i) {
    length(list.files(R_DIR, pattern = sprintf("^%02d_.*\\.R$", i))) > 0L
  }, logical(1)))
write_csv(code_inventory, "R_stage_inventory.csv")
for (s in c(0:6, 11, 12)) {
  if (!code_inventory[stage == s, present])
    add_log("ACTION", paste0("Stage ", sprintf("%02d", s)), "Source script absent; cannot certify all original cohort/exposure definitions")
}

# Manually verifiable rather than fabricated methods declarations.
# Keep author edits across reruns.
author_path <- file.path(OUT, "AUTHOR_CONFIRMATIONS_EDIT_ME.csv")
if (!file.exists(author_path)) {
  decisions <- data.table(
    topic = c("RFD institution", "RFD determination date", "RFD documentation/reference",
              "Protocol registered before analysis?", "Protocol registration URL/ID if any",
              "Were analyses prespecified, post hoc, or mixed?",
              "Data-access conditions and public-code permissions",
              "Medication record source and dictionary", "New-user washout and prior cross-class rules",
              "Stage 07 one-year sustained-record rule", "Stage 07 three-year persistent-record rule",
              "Geography variable patient_regional_location meaning",
              "External validation of meningioma codes?", "Can protocol be shared?",
              "Can analysis code be shared?"),
    response = c("Institutional RFD completed; non-human-subjects research", rep("", 14)),
    evidence_location = "", verified_by = "", verified_on = ""
  )
  write_csv(decisions, "AUTHOR_CONFIRMATIONS_EDIT_ME.csv")
}

rfd_language <- c(
  "PROPOSED ETHICS DECLARATION (review against the RFD record before use)",
  "This study underwent an institutional Research Feasibility Determination (RFD)",
  "and was determined not to constitute human subjects research. Accordingly,",
  "institutional review board review was not required.",
  "Add the institution, determination date, and reference number if assigned.",
  "DO NOT describe this as an IRB exemption or IRB approval unless that is",
  "what the institutional documentation specifically states.")
write_lines(rfd_language, "RFD_ethics_wording_TO_VERIFY.txt")

# Official documentation and reporting URLs; review the source checklist.
write_lines(c(
  "RECORD-PE (2018): https://pmc.ncbi.nlm.nih.gov/articles/PMC6234471/",
  "STROBE and RECORD items are inherited; this worksheet includes all three layers.",
  "TARGET (2025): https://www.bmj.com/content/390/bmj-2025-087179",
  "TARGET is mapped only to the target-trial-aligned DAY-0 sensitivity.",
  "The hypothetical trial and observed emulation require distinct descriptions.",
  "Do not mark any checklist item complete without verifying manuscript page references."
), "GUIDELINE_SOURCES.txt")

# DRAFT target-trial specification: contrast hypothetical protocol vs observed.
spec <- data.table(
  component = c("Eligibility", "Treatment strategies", "Assignment/time zero",
    "Follow-up", "Outcome", "Causal contrast", "Identification assumptions",
    "Analysis", "Treatment adherence", "Censoring"),
  hypothetical_trial_DRAFT = c(
    "Men aged 40-89 with BPH, no relevant prior meningeal tumor or prostate cancer; verify washout",
    "Initiate a 5-ARI vs initiate uroselective alpha-blocker, with assignment rules defined before follow-up",
    "Random allocation at treatment initiation (day 0)",
    "From randomization to outcome, death, loss or administrative end; exact definitions to confirm",
    "Incident meningioma with prespecified verification approach",
    "Effect of initial treatment strategy, independent of later adherence; define target population",
    "Randomization plus accurate event ascertainment, consistency and censoring assumptions",
    "Predefined survival analysis with effect estimates and absolute risk if estimable",
    "Not required for assignment-effect contrast; separate per-protocol question if pursued",
    "Specify competing events, nonadherence and censoring scheme"),
  observed_emulation_DRAFT = c(
    "BPH within preceding 365 days, no day 1-30 exposure used for assignment; Stage 14",
    "Index class based on observed finasteride/dutasteride versus tamsulosin/alfuzosin/silodosin",
    "Nonrandom observed treatment classification strictly at day 0, excluding same-day cross-class starts",
    "Zero-lag day-0 cohort; 1-y-lag sensitivity is separate; verify end-of-follow-up rules",
    "First ICD-10 D32.0/D32.9 or ICD-9 225.2 after risk-set entry",
    "ITT-like initiator association within the measured-covariate overlap population; not randomized ITT",
    "Conditional exchangeability, positivity, consistency, model specification, adequate EHR capture",
    "Propensity-score overlap-weighted Cox; confirm missingness handling and PH diagnostics",
    "Later switching not used to redefine initial day-0 assignment",
    "Censor at last observed EHR, with IPCW sensitivity; competing-death analysis"),
  verification_status = "AUTHOR_VALIDATION_REQUIRED"
)
write_csv(spec, "TARGET_Table_S1_expansion_WORKSHEET.csv")

# Original medication and persistence algorithms recovered in Stages 01/02/07.
codebook <- data.table(
  concept = c("BPH", "Prostate cancer", "Strict benign meningioma",
              "High-specificity intracranial meningioma", "Benign cranial-nerve neoplasm",
              "Benign spinal-cord neoplasm", "5-ARI drug class", "Comparator drugs",
              "Drug strength and dose", "One-year persistence", "Three-year persistence"),
  code_or_rule = c("ICD-10 N40*; ICD-9 600*", "ICD-10 C61; ICD-9 185",
    "ICD-10 D32.0, D32.9; ICD-9 225.2", "ICD-10 D32.0; ICD-9 225.2",
    "ICD-10 D33.3; ICD-9 225.1", "ICD-10 D33.4; ICD-9 225.3",
    "Finasteride RxNorm 25025; dutasteride RxNorm 228790; verify extract semantics",
    "Tamsulosin 77492; alfuzosin 17300; silodosin 720825 (RxNorm)",
    "Finasteride 1 mg excluded per Methods; verify code and attribution",
    "Stage 07: see explicit index-relative same-class date windows and switching exclusions", "Stage 07: see explicit index-relative same-class date windows and switching exclusions"),
  source = c(rep("Stage 14 / Methods", 6), rep("Stage 01/02 source recovered; verify TriNetX extract documentation", 3),
             "Stage 07 source recovered; Stage 11 adds utilization", "Stage 07 source recovered; Stage 11 adds utilization"),
  verified_against_source = c(rep("Stage 14 code reviewed; verify final raw source", 6),
                               rep("Recovered source; verify export dictionary", 5))
)
write_csv(codebook, "CODEBOOK_verify_before_supplement.csv")

# ---------------------------------------------------------------------------
# Aggregate-only read-only database reporting
# ---------------------------------------------------------------------------
if (!file.exists(DB_PATH)) {
  add_log("ACTION", "DuckDB", "Not available. Checklist drafts generated; no numerical follow-up/missingness output")
} else {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = DB_PATH, read_only = TRUE)
  tables <- DBI::dbListTables(con)
  primary_table <- "analysis_cohort_primary_boundary_fixed"
  if (!(primary_table %in% tables)) {
    add_log("FAIL", "Primary cohort", paste(primary_table, "missing; run completed Stage 16 first"))
  } else {
    fields <- DBI::dbListFields(con, primary_table)
    expected_fields <- c("patient_id", "exposure", "followup_days_strict_boundary",
                         "meningioma_event_strict_boundary")
    if (!all(expected_fields %in% fields)) {
      add_log("FAIL", "Primary fields", paste(setdiff(expected_fields, fields), collapse=", "))
    } else {
      # The primary cohort has 501k rows; calculations remain inside DuckDB.
      followup <- as.data.table(DBI::dbGetQuery(con, paste0(
        "SELECT exposure, COUNT(*) AS n, COUNT(DISTINCT patient_id) AS n_distinct, ",
        "SUM(CASE WHEN meningioma_event_strict_boundary=1 THEN 1 ELSE 0 END) AS events, ",
        "SUM(followup_days_strict_boundary)/365.25 AS person_years, ",
        "QUANTILE_CONT(followup_days_strict_boundary,0.25)/365.25 AS followup_q1_years, ",
        "MEDIAN(followup_days_strict_boundary)/365.25 AS followup_median_years, ",
        "QUANTILE_CONT(followup_days_strict_boundary,0.75)/365.25 AS followup_q3_years, ",
        "SUM(CASE WHEN meningioma_event_strict_boundary=0 THEN 1 ELSE 0 END) AS ",
        "ended_without_recorded_event, ",
        "SUM(CASE WHEN followup_days_strict_boundary IS NULL THEN 1 ELSE 0 END) AS missing_followup ",
        "FROM ", primary_table, " GROUP BY exposure ORDER BY exposure;")))
      followup[, crude_rate_per_100k_py := 1e5 * events / person_years]
      followup[, cohort := "Primary 1-y boundary-fixed"]
      setcolorder(followup, c("cohort", setdiff(names(followup), "cohort")))
      write_csv(followup, "Table_S8_primary_followup_BY_ARM.csv")
      # Numerical prose to check and edit into Results; not an adjusted risk.
      if (nrow(followup)==2L && sum(followup$n)==501445L &&
          all(c("5ARI") %in% followup$exposure)) {
        arm_a <- followup[exposure == "5ARI"][1]
        arm_b <- followup[exposure != "5ARI"][1]
        prose <- sprintf(
          paste0("In the boundary-corrected primary cohort, the 5-ARI group ",
          "contributed %.0f person-years of follow-up (median %.2f, ",
          "interquartile range %.2f-%.2f years) with %d recorded meningioma ",
          "events (crude rate %.1f per 100,000 person-years). The uroselective ",
          "alpha-blocker comparator contributed %.0f person-years ",
          "(median %.2f, interquartile range %.2f-%.2f years) ",
          "with %d events (crude rate %.1f per 100,000 person-years)."),
          arm_a$person_years, arm_a$followup_median_years,
          arm_a$followup_q1_years, arm_a$followup_q3_years,
          arm_a$events, arm_a$crude_rate_per_100k_py,
          arm_b$person_years, arm_b$followup_median_years,
          arm_b$followup_q1_years, arm_b$followup_q3_years,
          arm_b$events, arm_b$crude_rate_per_100k_py)
        write_lines(c("DRAFT descriptive results; VERIFY BEFORE INSERTING:", prose),
                    "RESULTS_followup_paragraph_FROM_DUCKDB.txt")
      }
      if (nrow(followup) != 2L || any(followup$n != followup$n_distinct) ||
          sum(followup$n) != 501445L || sum(followup$events) != 1456L) {
        add_log("FAIL", "Primary follow-up", "Treatment labels, one-entry-per-patient, or final counts do not reconcile")
      } else {
        add_log("PASS", "Primary follow-up", "N=501445, events=1456, two groups and one row/patient")
      }
      add_log("NOTE", "Outcome-free observation endings",
        "NOT equivalent to adjudicated deaths or confirmed dropout; last encounter cannot resolve cause")

      # Audit missing values of raw baseline fields. A coded 0 is NOT missing.
      # For categorical fields, report NULL/empty and Unknown/Unspecified separately.
      covars <- c("age_at_index","race","ethnicity","marital_status",
        "patient_regional_location","patient_source_id","hypertension", "diabetes",
        "obesity","ckd","smoking","cad","heart_failure","cerebrovascular_disease",
        "hyperlipidemia","osteoporosis","alopecia","urinary_retention","uti",
        "hematuria","bladder_stone","n_encounters_prior_365",
        "n_encounters_prior_1095","n_ambulatory_prior_365","n_ed_prior_365",
        "n_inpatient_prior_365")
      available <- intersect(covars, fields)
      missing_covars <- setdiff(covars, available)
      if (length(missing_covars)) add_log("ACTION", "Missingness inventory",
          paste("Columns absent:", paste(missing_covars, collapse = ", ")))
      categorical <- c("race","ethnicity","marital_status","patient_regional_location",
                       "patient_source_id")
      # Each group/field aggregation is SQL, not patient-level export.
      miss <- rbindlist(lapply(available, function(v) {
        var_sql <- paste0('"', v, '"')
        sql <- paste0("SELECT exposure, COUNT(*) AS n, ",
          "SUM(CASE WHEN ", var_sql, " IS NULL OR TRIM(CAST(",var_sql,
          " AS VARCHAR))='' THEN 1 ELSE 0 END) AS null_or_blank, ",
          "SUM(CASE WHEN UPPER(TRIM(CAST(",var_sql,
          " AS VARCHAR))) IN ('UNKNOWN','UNSPECIFIED','NOT KNOWN','NOT REPORTED') ",
          "THEN 1 ELSE 0 END) AS explicitly_unknown ",
          "FROM ", primary_table, " GROUP BY exposure;")
        dd <- as.data.table(DBI::dbGetQuery(con, sql))
        dd[, `:=`(variable = v,
                  interpretation = if (v %in% categorical)
                    "Null/blank and explicit unknown separately; inspect category codes"
                  else "NULL/blank only; binary 0 is not treated as missing")]
        dd
      }), use.names = TRUE, fill = TRUE)
      miss[, `:=`(null_or_blank_pct = 100*null_or_blank/n,
                  unknown_pct = 100*explicitly_unknown/n)]
      # Suppress small nonzero count cells in exported summaries, and withhold
      # associated proportions. Institutional disclosure review is still needed.
      miss[, missing_report := ifelse(null_or_blank %in% 1:10, "<11", as.character(null_or_blank))]
      miss[, unknown_report := ifelse(explicitly_unknown %in% 1:10, "<11", as.character(explicitly_unknown))]
      miss[null_or_blank %in% 1:10, null_or_blank_pct := NA_real_]
      miss[explicitly_unknown %in% 1:10, unknown_pct := NA_real_]
      write_csv(miss[, .(variable, exposure, n, missing_report,
                          null_or_blank_pct, unknown_report, unknown_pct, interpretation)],
                "MISSINGNESS_BY_ARM_REVIEW_DISCLOSURE.csv")
      add_log("CREATED", "Missingness", paste(length(available), "available baseline variables by arm"))

      # Geographic coding QC: do not assume a label implies country of care.
      if ("patient_regional_location" %in% fields) {
        regions <- as.data.table(DBI::dbGetQuery(con, paste0(
          "SELECT CAST(patient_regional_location AS VARCHAR) AS region_code, ",
          "COUNT(*) AS n FROM ",primary_table," GROUP BY 1 ORDER BY n DESC;")))
        # Hide small counts; do not disclose suppressed cells via totals.
        regions[, n_release_view := ifelse(n < 11L, "<11 / suppressed", as.character(n))]
        write_csv(regions[, .(region_code, n_release_view)],
                  "REGION_CODING_REVIEW_no_small_cells.csv")
        add_log("ACTION", "Geography", "Verify meaning of any Ex-US/nonstandard category in TriNetX data dictionary")
      }
    }
  }

  # Describe the separate day-zero target-trial-aligned cohort, if available.
  day0_table <- "analysis_cohort_tte_index_itt_zero"
  if (day0_table %in% tables) {
    f0 <- DBI::dbListFields(con, day0_table)
    if (all(c("exposure","meningioma_event_strict","followup_days_strict") %in% f0)) {
      d0 <- as.data.table(DBI::dbGetQuery(con, paste0(
        "SELECT exposure, COUNT(*) AS n, SUM(meningioma_event_strict) AS events, ",
        "SUM(followup_days_strict)/365.25 AS person_years, ",
        "MEDIAN(followup_days_strict)/365.25 AS median_followup_years, ",
        "QUANTILE_CONT(followup_days_strict,0.25)/365.25 AS followup_q1_years, ",
        "QUANTILE_CONT(followup_days_strict,0.75)/365.25 AS followup_q3_years ",
        "FROM ",day0_table," GROUP BY exposure ORDER BY exposure;")))
      d0[, crude_rate_per_100k_py := 1e5*events/person_years]
      write_csv(d0,"TARGET_day0_followup_BY_ARM.csv")
      # Day-zero missingness cannot be assumed identical to primary cohort.
      d0_covars <- intersect(c("age_at_index", "race", "ethnicity",
                              "marital_status", "patient_regional_location",
                              "patient_source_id", "n_encounters_prior_365"), f0)
      if (length(d0_covars)) {
        d0_missing <- rbindlist(lapply(d0_covars, function(v) {
          q <- paste0('"', v, '"')
          z <- paste0("SELECT exposure, COUNT(*) AS n, ",
            "SUM(CASE WHEN ", q, " IS NULL OR TRIM(CAST(",q,
            " AS VARCHAR))='' THEN 1 ELSE 0 END) AS missing_n ",
            "FROM ",day0_table," GROUP BY exposure;")
          ans <- as.data.table(DBI::dbGetQuery(con,z))
          ans[, variable := v]
          ans
        }), use.names=TRUE, fill=TRUE)
        d0_missing[, missing_pct := 100*missing_n/n]
        d0_missing[, missing_n_report := ifelse(missing_n %in% 1:10,"<11",as.character(missing_n))]
        d0_missing[missing_n %in% 1:10, missing_pct := NA_real_]
        write_csv(d0_missing[,.(variable,exposure,n,missing_n_report,missing_pct)],
                  "TARGET_day0_missingness_BY_ARM.csv")
      }
      add_log("CREATED", "Day-zero TARGET descriptives", "Not to be pooled with the 1-y primary landmark cohort")
    } else add_log("ACTION", "Day-zero fields", "Missing required columns in Stage 14 table")
  } else add_log("ACTION", "TARGET day-zero", "Stage 14 day-zero cohort missing; cannot generate its descriptive table")
  DBI::dbDisconnect(con, shutdown = TRUE)
}

# ---------------------------------------------------------------------------
# Validate existing aggregate numerical files without refitting models.
# ---------------------------------------------------------------------------
# Stage 14 already created a manuscript target-trial specification table.
trial14 <- maybe_csv(file.path(WORK_DIR,"qc_methodologic_validity",
                               "14_target_trial_specification.csv"))
if (!is.null(trial14)) {
  write_csv(trial14,"SOURCE_14_target_trial_specification.csv")
  add_log("CREATED","Stage 14 target-trial source","Review beside TARGET S1 expansion worksheet")
} else add_log("ACTION","Stage 14 target-trial source","Missing 14_target_trial_specification.csv")
incfile <- file.path(WORK_DIR,"analysis_stage16_high_specificity","16_incidence_rates.csv")
inc <- maybe_csv(incfile)
if (!is.null(inc)) {
  inc <- inc[grepl("^Boundary-fixed primary 1-y landmark$",analysis)]
  if (nrow(inc) == 2L) {
    write_csv(inc,"PRIMARY_crude_and_weighted_incidence_FROM_STAGE16.csv")
    if (exists("followup")) {
      check_n <- sum(inc$n)==sum(followup$n) && sum(inc$events)==sum(followup$events)
      check_py <- isTRUE(all.equal(sum(inc$person_years), sum(followup$person_years), tolerance=1e-7))
      add_log(if (check_n && check_py) "PASS" else "FAIL", "Stage 16 incidence",
              "Cross-check group counts/events and total person-years against Stage 20")
    }
  }
} else add_log("ACTION","Stage 16 incidence","Missing 16_incidence_rates.csv")

fdrpath <- file.path(WORK_DIR,"analysis_stage16_high_specificity","16_expanded_exploratory_FDR.csv")
fdr <- maybe_csv(fdrpath)
if (!is.null(fdr)) {
  mandatory <- c("analysis","weighted_p","q_BH_expanded")
  if (all(mandatory %in% names(fdr)) && nrow(fdr) == 16L &&
      !anyDuplicated(fdr$analysis) && all(is.finite(fdr$weighted_p))) {
    fdr[, independently_recalculated_q := stats::p.adjust(weighted_p, method="BH")]
    delta <- max(abs(fdr$independently_recalculated_q-fdr$q_BH_expanded))
    write_csv(fdr[, .(analysis, weighted_p, q_BH_expanded,
                       independently_recalculated_q)], "FDR_AUDIT_UNROUNDED.csv")
    add_log(if (delta < 1e-10) "PASS" else "FAIL", "FDR 16-model family",
            sprintf("Maximum difference: %.3g", delta))
  } else add_log("FAIL","FDR input","Missing columns, wrong family size, duplicated name or invalid p")
} else add_log("ACTION","FDR","Stage 16 exploratory FDR file absent")

for (nm in c("15_ipcw_diagnostics.csv","15_ipcw_results.csv")) {
  p <- file.path(WORK_DIR, "analysis_methodologic_validity", nm)
  d <- maybe_csv(p)
  if (!is.null(d)) {
    write_csv(d,paste0("SOURCE_",nm))
    add_log("CREATED",nm,"Copied aggregate-only diagnostics; interpret weight truncation")
  } else add_log("ACTION",nm,"Missing; cannot report IPCW diagnostics")
}

# Audit presence of PH diagnostics only; never pretend model assumptions passed.
ph_candidates <- unlist(lapply(c(file.path(WORK_DIR,"analysis"),
  file.path(WORK_DIR,"analysis_final"),
  file.path(WORK_DIR,"analysis_methodologic_validity"),
  file.path(WORK_DIR,"analysis_stage16_high_specificity")),
  function(d) if (dir.exists(d)) list.files(d, pattern="(ph|zph).*\\.csv$",
                                            full.names=TRUE, ignore.case=TRUE)))
if (length(ph_candidates)) {
  write_lines(ph_candidates, "PH_diagnostic_file_inventory.txt")
  add_log("ACTION","Proportional hazards","Diagnostic filenames located; inspect model-specific results before claiming assumption met")
} else add_log("ACTION","Proportional hazards","No PH diagnostic CSV found. Check model outputs/cox.zph before claiming assumptions met")

# Stage 17 v8+ code should be the authoritative table builder; flag Table S7.
rare <- maybe_csv(file.path(WORK_DIR,"manuscript_ready_final","tables",
                            "Table_S7_Rare_hard_endpoints.csv"))
if (!is.null(rare) && "N" %in% names(rare)) {
  rare_n <- suppressWarnings(as.integer(gsub(",", "", rare$N)))
  if (any(is.finite(rare_n) & rare_n != 501445L))
    add_log("ACTION","Table S7","Some rare-endpoint models retain pre-boundary N=501446. Correct/refit or document transparently")
}

# Study-design schematic (not a new estimate); Stage 07 defines persistence windows.
try({
  grDevices::pdf(file.path(OUT,"Figure_S2_design_schematic_DRAFT.pdf"),
                 width=10.5,height=4.3, onefile=FALSE, useDingbats=FALSE)
  old <- par(mar=c(1.3,1.3,1.8,1.3))
  plot.new(); plot.window(xlim=c(0,10),ylim=c(0,3.7))
  text(0.1,3.43,"Study designs and time zero",adj=0,font=2,cex=1.2)
  rect(0.1,2.35,2.5,2.72,col="gray90",border=NA)
  rect(2.5,2.35,5.5,2.72,col="gray78",border=NA)
  rect(5.5,2.35,9.7,2.72,col="gray90",border=NA)
  text(1.3,2.53,"Baseline",cex=.92)
  text(4.0,2.53,"1-year eligibility/latency",cex=.8)
  text(7.6,2.53,"Primary outcome follow-up",cex=.85)
  text(.1,3.0,"Primary active-comparator 1-year landmark",adj=0,font=2,cex=.95)
  text(2.5,2.10,"Index / day 0",cex=.83); text(5.5,2.10,"Year 1 risk-set entry",cex=.83)
  rect(.1,.9,2.5,1.27,col="gray90",border=NA)
  rect(2.5,.9,9.7,1.27,col="gray82",border=NA)
  text(1.3,1.08,"Baseline",cex=.92)
  text(6.0,1.08,"Outcome follow-up begins at day 0",cex=.9)
  text(.1,1.55,"Target-trial-aligned day-zero sensitivity",adj=0,font=2,cex=.95)
  text(2.5,.68,"Index / day 0",cex=.83)
  text(.1,.18,"Schematic only. Washout and 1-/3-year persistence definitions require source-code verification.",
       adj=0,cex=.76)
  par(old); grDevices::dev.off()
}, silent=TRUE)

# Editable drafts, NOT a submission-ready supplement or a certification.
sup_methods <- c(
  "SUPPLEMENTARY METHODS: AUTHOR DRAFT / NOT VERIFIED",
  "", "Reporting frameworks",
  "The primary study was an active-comparator, new-user EHR cohort with a",
  "one-year landmark. The distinct day-zero sensitivity was constructed to",
  "align treatment classification and the start of follow-up at treatment",
  "initiation. Reporting was guided by RECORD-PE (including STROBE and RECORD)",
  "and, for the day-zero analysis only, by TARGET.",
  "", "Medication exposure definitions - INCOMPLETE UNTIL STAGE 00-06 AND 11",
  "The drug classes are finasteride/dutasteride and uroselective alpha-blockers",
  "(tamsulosin/alfuzosin/silodosin). Verify the medication vocabulary, whether",
  "records represent prescribing or dispensation, washout, index-date rules,",
  "same-class persistence windows, multi-class switching, combinations and",
  "finasteride 1-mg exclusions. Do not copy this paragraph as complete methods.",
  "", "Risk sets and follow-up",
  "The boundary-corrected primary analysis excludes tumors diagnosed on or",
  "before the 365-day landmark; events begin strictly after the landmark.",
  "For each arm, report person-years, event counts, crude incidence rates,",
  "median (IQR) follow-up and ending without recorded outcome (Table S8).",
  "An absence of a recorded event by last EHR encounter does not establish",
  "death or transfer; additional censoring reasons may not be identifiable.",
  "", "Missingness and model assumptions",
  "Use the Stage 20 missingness export to report observed missing categories",
  "separately from absent comorbidity codes. Verify original confounder coding",
  "and whether unrecorded diagnoses were assigned zero in the analysis.",
  "Describe propensity-score overlap, weighted covariate balance, proportional",
  "hazards assessment, censoring-weight tails, truncation and sensitivity",
  "analyses from original model outputs. Never claim untested assumptions pass.",
  "", "Ethics", paste(rfd_language[2:4],collapse=" "),
  "", "Limitations for interpretation",
  "A positive negative-control association warrants caution about residual",
  "confounding/differential ascertainment. Post-initiation persistence is not",
  "randomized or a verified per-protocol causal estimate. Overlap weighting",
  "targets the population with substantial treatment-assignment overlap."
)
write_lines(sup_methods, "SUPPLEMENTARY_METHODS_AUTHOR_DRAFT.txt")

insertions <- c(
  "MANUSCRIPT INSERTIONS / DO NOT PASTE BEFORE SOURCE CHECKS",
  "", "Methods - Reporting standards:",
  "Reporting followed RECORD-PE (including its STROBE and RECORD components).",
  "The separate target-trial-aligned day-zero analyses were additionally",
  "described with reference to the TARGET reporting guideline.",
  "", "Ethics:", paste(rfd_language[2:4],collapse=" "),
  "", "Results - follow-up:",
  "Insert the arm-specific person-years, median (IQR) years, events and crude",
  "rates from Table_S8_primary_followup_BY_ARM.csv after QC passes.",
  "Do not report events/total treated as a 5-year risk.",
  "", "Data availability:",
  "The underlying TriNetX data are controlled-access EHR data and are not",
  "deposited with this publication [CONFIRM TERMS]. Deidentified analytic",
  "code and supporting aggregate data may be shared to the extent permitted",
  "by institutional and TriNetX terms [AUTHOR CONFIRMATION REQUIRED].",
  "", "Protocol registration:",
  "State explicitly whether the protocol was prospectively registered; do not",
  "claim prespecification without dated contemporaneous documentation."
)
write_lines(insertions, "MANUSCRIPT_INSERTIONS_TO_VERIFY.txt")

add_log("ACTION","Checklist completion",
  "Fill item-by-item manuscript/supplement PAGE numbers after final document pagination")
add_log("ACTION","Compliance scope",
  "New descriptive statistics + worksheets do not make a study automatically RECORD-PE/TARGET compliant")
add_log("ACTION","Data disclosure",
  "Review all aggregate exports under institutional and TriNetX sharing requirements before circulation")
write_lines(c("STAGE 20 RUN REPORT", paste("Generated:",Sys.time()),
              paste("R version:",R.version.string), "", log_lines),
            "STAGE20_QC_AND_ACTION_LOG.txt")

cat("\nStage 20 completed. Review:", normalizePath(OUT, winslash="/"),"\n")
cat("Check STAGE20_QC_AND_ACTION_LOG.txt before reusing numerical exports.\n")
