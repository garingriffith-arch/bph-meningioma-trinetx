# Final closeout orchestrator. Does not rerun Stages 00-23 or overwrite DuckDB.
# Default runs 24 (exact-boundary PH), 25 (four rare endpoints), 26 (aggregate audit).
# To review existing outputs only: set BPH_CLOSEOUT_MODE=audit.
# To run a single stage: BPH_CLOSEOUT_MODE=24, 25, or 26.
# Existing Stage 24/25 outputs are not overwritten.

R_DIR <- Sys.getenv("BPH_R_DIR", "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R")
WORK_DIR <- Sys.getenv("BPH_WORK_DIR", "C:/TriNetX/bph_m_work")
mode <- tolower(Sys.getenv("BPH_CLOSEOUT_MODE", "all"))
if (!mode %in% c("all", "audit", "24", "25", "26"))
  stop("BPH_CLOSEOUT_MODE must be all, audit, 24, 25 or 26")

run_one <- function(stage, filename, sentinel=NULL) {
  f <- file.path(R_DIR, filename)
  if (!file.exists(f)) stop("Missing stage script: ",f)
  if (!is.null(sentinel) && file.exists(sentinel)) {
    message("Stage ",stage," output already exists; reusing: ",sentinel)
    return(invisible(FALSE))
  }
  message("Running Stage ",stage," (",filename,")")
  source(f, local=new.env(parent=globalenv()))
  if (!is.null(sentinel) && !file.exists(sentinel)) stop("Stage ",stage," did not create expected output: ",sentinel)
  invisible(TRUE)
}
if (mode %in% c("all","24"))
  run_one(24,"24_verify_boundary_fixed_PH.R",
    file.path(WORK_DIR,"reporting_boundary_PH_stage24","24_boundary_fixed_PH_diagnostics.csv"))
if (mode %in% c("all","25"))
  run_one(25,"25_reconcile_rare_endpoints.R",
    file.path(WORK_DIR,"reporting_rare_stage25","25_rare_boundary_aligned_models_unrounded.csv"))
if (mode %in% c("all","audit","26"))
  run_one(26,"26_finalize_analysis_audit.R")
message("Closeout pipeline finished. Inspect reporting_final_closeout_stage26/26_closeout_register.csv")
