# Stage 26: Aggregate-only release audit and RECORD-PE/TARGET closeout register.
# No model refits, no database writes, no fabricated IRB or checklist approvals.

run_stage_26 <- function() {
  options(stringsAsFactors = FALSE, scipen = 999)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table")
  suppressPackageStartupMessages(library(data.table))
  WORK_DIR <- Sys.getenv("BPH_WORK_DIR", "C:/TriNetX/bph_m_work")
  R_DIR <- Sys.getenv("BPH_R_DIR", "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R")
  OUT_DIR <- file.path(WORK_DIR, "reporting_final_closeout_stage26")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  csv <- function(...) file.path(WORK_DIR, ...)
  required_paths <- c(
    csv("analysis_stage16_high_specificity", "16_model_results.csv"),
    csv("analysis_stage16_high_specificity", "16_expanded_exploratory_FDR.csv"),
    csv("manuscript_ready_final", "analysis_snapshot", "FINAL_results_snapshot.csv"),
    csv("analysis_methodologic_validity", "15_ipcw_results.csv")
  )
  missing <- required_paths[!file.exists(required_paths)]
  if (length(missing)) stop("Missing finalized input:\n", paste(missing, collapse = "\n"))
  model16 <- fread(required_paths[1]); fdr <- fread(required_paths[2])
  snap <- fread(required_paths[3]); ipcw <- fread(required_paths[4])
  decisions <- list()
  add <- function(item, state, evidence, next_action="") {
    decisions[[length(decisions)+1L]] <<- data.table(
      item = item, status = state, evidence = evidence, action = next_action)
  }
  nearly <- function(a, b, tol=1e-9) is.finite(a) && is.finite(b) && abs(a-b) <= tol
  primary <- model16[analysis == "Boundary-fixed primary 1-y landmark"]
  if (nrow(primary)!=1L) stop("Primary Stage 16 result not unique")
  main <- snap[Analysis == "Primary 1-y landmark"]
  if (nrow(main)!=1L) stop("Primary snapshot result not unique")
  cohort_ok <- primary$n==501445L && primary$n_treated==51390L &&
    primary$n_control==450055L && primary$events==1456L &&
    primary$events_treated==206L && primary$events_control==1250L
  estimates_ok <- nearly(primary$weighted_HR,main$HR) &&
    nearly(primary$weighted_CI_low, main$CI_low) &&
    nearly(primary$weighted_CI_high, main$CI_high)
  if (!cohort_ok || !estimates_ok) stop("Final primary cohort/results mismatch")
  add("Final primary cohort and effect", "PASS",
      "N=501445, 206/1250 events; Stage 16 vs final snapshot exact agreement")
  day0 <- snap[grepl("day-0 ITT-like", Analysis, fixed=TRUE)]
  day0_ok <- nrow(day0)==1L && day0$N==660012L && day0$Events==1908L
  if (!day0_ok) stop("Day-zero model snapshot not reconciled")
  add("TARGET day-zero model", "PASS", "N=660012, 1908 outcome events")

  required_fdr <- c("analysis","weighted_p","q_BH_expanded")
  if (!all(required_fdr %in% names(fdr))) stop("FDR columns changed")
  if (nrow(fdr)!=16L || anyDuplicated(fdr$analysis)) stop("Expected 16 unique FDR hypotheses")
  recomputed <- p.adjust(fdr$weighted_p, method="BH")
  if (any(!is.finite(recomputed)) || max(abs(recomputed-fdr$q_BH_expanded)) > 1e-10)
    stop("Benjamini-Hochberg values do not reconcile")
  fdr[, q_recomputed := recomputed]
  fwrite(fdr, file.path(OUT_DIR, "26_FDR_family_verified_unrounded.csv"))
  add("16-model exploratory FDR family", "PASS", "BH recalculation agrees with saved values")
  if (nrow(ipcw)<1L) stop("No IPCW results")
  add("IPCW output", "PASS", "Finalized Stage 15 output present; see manuscript snapshot for reported effect")

  flowfile <- csv("reporting_TARGET_day0", "TARGET_day0_flow_validation.csv")
  if (!file.exists(flowfile)) {
    add("Day-zero participant flow", "MISSING", "Stage 21 validation not found", "Run Stage 21")
  } else {
    flow <- fread(flowfile)
    if (!"passed" %in% names(flow) || anyNA(flow$passed) || !all(as.logical(flow$passed)))
      stop("Stage 21 day-zero cumulative flow failed")
    add("Day-zero participant flow", "PASS", "All treatment arms match the Stage 14 table")
  }
  baseline <- csv("reporting_TARGET_day0", "TARGET_day0_baseline_characteristics.csv")
  # Stage 21 output filenames can vary by patch; capture actual inventory.
  files21 <- if (dir.exists(csv("reporting_TARGET_day0")))
    list.files(csv("reporting_TARGET_day0"), full.names=FALSE) else character()
  add("TARGET baseline reporting", if (any(grepl("baseline", files21, ignore.case=TRUE))) "PRESENT" else "MISSING",
      paste(files21[grepl("baseline", files21, ignore.case=TRUE)], collapse="; "),
      "Retain final day-zero table and specify baseline is unweighted")

  ph_file <- csv("reporting_boundary_PH_stage24", "24_boundary_fixed_PH_diagnostics.csv")
  if (!file.exists(ph_file)) {
    add("Boundary-aligned PH diagnostics", "MISSING", "Final Stage 16 exact-model PH not run", "Run Stage 24")
  } else {
    ph <- fread(ph_file)
    if (!all(c("analysis","p","n","events","matched_stage16") %in% names(ph)))
      stop("Stage 24 PH schema mismatch")
    if (anyNA(ph$matched_stage16) || !all(as.logical(ph$matched_stage16)) || anyNA(ph$p))
      stop("Stage 24 model reconciliation or PH results invalid")
    expected_core <- c("Boundary-fixed primary 1-y landmark",
      "Boundary-fixed target-trial day-0 assignment, 1-y landmark",
      "Boundary-fixed strict endpoint: strict sustained 1-y, finasteride vs uro-alpha",
      "Boundary-fixed strict endpoint: persistent 3-y, finasteride vs uro-alpha")
    if (!all(expected_core %in% ph$analysis)) stop("Stage 24 missing core exact-boundary PH models")
    # Tests inform interpretation; P <= 0.05 is not silently ignored.
    note <- if (any(ph$p <= 0.05)) "At least one treatment PH test P<=0.05: review time variation" else
      "No treatment-effect test P<=0.05; this does not prove proportional hazards"
    add("Boundary-aligned PH diagnostics", if (any(ph$p<=0.05)) "REVIEW" else "PASS",note)
    fwrite(ph, file.path(OUT_DIR, "26_boundary_PH_for_reporting.csv"))
  }

  rare_path <- csv("reporting_rare_stage25", "25_rare_boundary_aligned_models_unrounded.csv")
  rare_qc_path <- csv("reporting_rare_stage25", "25_archived_replay_and_boundary_reconciliation.csv")
  if (!file.exists(rare_path) || !file.exists(rare_qc_path)) {
    add("Rare endpoint denominator", "MISSING", "Table S11 still uses historical N=501446", "Run Stage 25, then replace Table S11")
  } else {
    rare <- fread(rare_path); rare_qc <- fread(rare_qc_path)
    if (nrow(rare)!=4L || any(rare$n != 501445L) ||
      nrow(rare_qc)!=4L || any(!as.logical(rare_qc$replay_pass)))
      stop("Stage 25 rare endpoint verification failed")
    add("Rare endpoint denominator", "PASS", "All four historical fits reproduced then N=501445 boundary-aligned", "Replace Supplementary Table S11")
  }

  # The local institutional data dictionary, RFD metadata and protocol history
  # are not inferable from model CSVs. These statuses are never auto-PASS.
  for (entry in list(
    c("Geographic field definition", "Obtain institutional Dataworks dictionary; do not reinterpret Ex-US"),
    c("Medication source semantics", "Confirm order/medication-list/administration versus dispensing"),
    c("RFD documentation", "Record institution, determination date, and determination identifier if issued"),
    c("Registration and prespecification history", "Report original protocol dates and exploratory timing accurately"),
    c("Authorship and declarations", "Verify author order, affiliations, corresponding author, funding and conflicts"),
    c("RECORD-PE/TARGET page mapping", "Add final manuscript and supplement page numbers after layout lock"),
    c("TriNetX permitted disclosure", "Confirm approved external-sharing terms for aggregated results and code")
  )) add(entry[1], "AUTHOR_CONFIRM", "Cannot be inferred from numerical files", entry[2])

  closed <- rbindlist(decisions, fill=TRUE)
  fwrite(closed, file.path(OUT_DIR,"26_closeout_register.csv"))
  exported <- snap[, .(Category, Analysis, N, Events, HR, CI_low, CI_high, p, q)]
  fwrite(exported, file.path(OUT_DIR,"26_final_results_locked.csv"))

  methods <- c(
    "# Manuscript text for author review (not an automatic manuscript edit)",
    "",
    "## Ethics",
    "This study underwent an institutional Research Feasibility Determination (RFD) and was determined not to constitute human subjects research; institutional review board review was therefore not required.",
    "Add institution and reference/date only from the actual determination.",
    "",
    "## Design and reporting",
    "This active-comparator, new-user observational cohort was reported with reference to RECORD-PE, incorporating STROBE and RECORD. Separately, the day-zero initiation analyses were described using the applicable TARGET items. The primary one-year landmark analysis is not described as a fully emulated day-zero randomized trial.",
    "",
    "## PH diagnostics",
    if (file.exists(ph_file)) {
      ph <- fread(ph_file)
      paste0("Treatment-effect proportional-hazards tests were evaluated on ",nrow(ph),
        " archived Stage 16 model specifications refitted on their exact boundary-corrected cohorts. ",
        "Report per-analysis results and limitations from 26_boundary_PH_for_reporting.csv; ",
        "a nonsignificant test is not confirmation of the PH assumption.")
    } else "Insert the exact-cohort Stage 24 Schoenfeld treatment-effect tests after they are run.",
    "",
    "## Rare outcomes",
    if (file.exists(rare_path))
      "Supplementary Table S11 was refitted using the 501,445-person boundary-corrected primary cohort with original rare-endpoint definitions unchanged."
    else "Rare endpoints retain historical 501,446-person denominators until Stage 25 refitting; do not imply boundary alignment.",
    "",
    "## Limitations requiring explicit discussion",
    "Primary landmark selection, incomplete medication capture, persistence-related surveillance/selection bias, unmeasured confounding, a positive cranial-nerve neoplasm falsification signal, informative EHR observation, and source-unknown geographic category definitions."
  )
  writeLines(methods, file.path(OUT_DIR, "26_manuscript_text_for_review.md"), useBytes=TRUE)

  readme <- c(
    "# Final analysis closeout: Stage 26",
    "",
    paste0("Generated: ", format(Sys.time(), tz="UTC",usetz=TRUE)),
    "Source: already finalized Stage 15/16/17/21/24/25 aggregate outputs.",
    "No clinical-data extracts, database modifications, additional hypothesis tests, or model refits.",
    "",
    "Automated PASS refers to cross-file numerical reconciliation, not institutional endorsement or validation of outcome coding.",
    "Missing Stage 24 or 25 is reported as MISSING; the script may still complete.",
    "Author confirmations are listed separately and never invented.",
    "",
    paste0("PASS: ",sum(closed$status=="PASS"),"; PRESENT: ",sum(closed$status=="PRESENT"),
           "; REVIEW: ",sum(closed$status=="REVIEW"),"; MISSING: ",sum(closed$status=="MISSING"),
           "; AUTHOR_CONFIRM: ",sum(closed$status=="AUTHOR_CONFIRM")),
    "",
    "Replace updated rare endpoint table in Word, reconcile all in-text citations, then freeze page numbering and complete checklists."
  )
  writeLines(readme, file.path(OUT_DIR,"26_SUMMARY.md"))
  message("Stage 26 completed: ", OUT_DIR)
}
run_stage_26()
