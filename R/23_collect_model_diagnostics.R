# 23_collect_model_diagnostics.R
# Collect existing aggregate proportional-hazards and model diagnostics.
# Does not access the patient-level database, refit models, or change source outputs.
options(stringsAsFactors = FALSE)
work <- Sys.getenv("BPH_WORK_DIR", "C:/TriNetX/bph_m_work")
out <- file.path(work, "reporting_model_diagnostics_review")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
rel <- c("analysis/06_ph_diagnostics.csv", "analysis_final/10_final_ph_diagnostics.csv",
         "analysis_methodologic_validity/15_methodologic_ph_diagnostics.csv",
         "analysis_methodologic_validity/15_ipcw_diagnostics.csv")
manifest <- lapply(rel, function(r) {
  f <- file.path(work, r)
  if (!file.exists(f)) return(data.frame(file=r, status="NOT_FOUND", bytes=NA_integer_))
  x <- utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
  # Diagnostic tables are usually aggregate, but inspect column names before exporting.
  banned <- grepl("patient_id|mrn|person_id|date_of_birth|dob|name|address", names(x), ignore.case=TRUE)
  if (any(banned)) stop("Sensitive column identified in ", r, "; refusing export")
  target <- file.path(out, basename(f))
  if (file.exists(target)) stop("Output exists; move or archive it before rerunning: ", target)
  utils::write.csv(x, target, row.names=FALSE, na="")
  data.frame(file=r, status="COPIED_AGGREGATE", bytes=file.info(f)$size)
})
utils::write.csv(do.call(rbind, manifest), file.path(out,"diagnostic_manifest.csv"), row.names=FALSE)
writeLines(c("Stage 23 aggregate diagnostic collection", "No models were refitted.",
 "Confirm PH test specification, model name, and p values before making manuscript claims.",
 "Review files for disclosure compliance before sharing."),file.path(out,"READ_ME.txt"))
cat("Diagnostic review files:",out,"\n")
