# 19_collect_final_resident_audit.R
# BPH-meningioma project: collect aggregate results for the final pre-resident audit.
# Run after Stages 15, 17, and 18. DOES NOT fit models, query DuckDB, or edit results.
# Only named aggregate CSVs and optional Stage 11/12 R source are collected.
# The raw patient data, DuckDB, and patient-level QC files are deliberately EXCLUDED.

options(stringsAsFactors = FALSE, scipen = 999)
R_DIR <- 'C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R'
WORK_DIR <- 'C:/TriNetX/bph_m_work'
DEST <- file.path(WORK_DIR, 'resident_audit_upload')
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)

wanted <- c(
  'analysis_methodologic_validity/15_ipcw_results.csv',
  'analysis_methodologic_validity/15_ipcw_diagnostics.csv',
  'analysis_methodologic_validity/15_methodologic_model_results.csv',
  'analysis_stage16_high_specificity/16_expanded_exploratory_FDR.csv',
  'analysis_stage16_high_specificity/16_model_results.csv',
  'analysis_stage16_high_specificity/16_incidence_rates.csv',
  'qc_stage16_high_specificity/16_landmark_boundary_reconciliation.csv',
  'manuscript_ready_final/analysis_snapshot/FINAL_results_snapshot.csv',
  'manuscript_ready_final/analysis_snapshot/FINAL_boundary_aligned_sensitivity_models.csv',
  'manuscript_ready_final/analysis_snapshot/FINAL_boundary_aligned_sensitivity_incidence.csv',
  'manuscript_ready_final/tables/Table_2_Primary_robustness_results.csv',
  'manuscript_ready_final/tables/Table_3_Exploratory_persistence_results.csv',
  'manuscript_ready_final/tables/Table_S6_Complete_exploratory_FDR_family.csv',
  'manuscript_ready_final/tables/Table_S7_Rare_hard_endpoints.csv',
  'manuscript_ready_final/MANUSCRIPT_FIGURE_LEGENDS.txt',
  'pre_resident_QC/QC_REPORT.txt',
  'pre_resident_QC/borderline_FDR_unrounded.csv',
  'pre_resident_QC/boundary_fixed_cohort_counts.csv'
)

log <- c('FINAL PRE-RESIDENT AGGREGATE FILE COLLECTION',
         paste('Run at', format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z')), '')
copied <- character()
for (relative in wanted) {
  source <- file.path(WORK_DIR, relative)
  target <- file.path(DEST, relative)
  if (!file.exists(source)) {
    log <- c(log, paste('[MISSING]', relative))
    next
  }
  # Guard against accidentally collecting unusually large or patient-level extracts.
  if (file.info(source)$size > 15 * 1024^2) {
    log <- c(log, paste('[SKIPPED >15MB]', relative))
    next
  }
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, target, overwrite = TRUE)) {
    log <- c(log, paste('[COPY FAILED]', relative))
  } else {
    copied <- c(copied, relative)
    log <- c(log, paste('[COPIED]', relative))
  }
}

# Source scripts are useful for auditing the postbaseline persistence definitions.
# Restrict this collector to two specifically named code stages.
if (dir.exists(R_DIR)) {
  source_scripts <- list.files(R_DIR, pattern = '^(11|12)_.*\\.[Rr]$',
                               recursive = TRUE, full.names = TRUE)
  for (src in source_scripts) {
    target <- file.path(DEST, 'R_source', basename(src))
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (file.copy(src, target, overwrite = TRUE)) {
      copied <- c(copied, file.path('R_source', basename(src)))
      log <- c(log, paste('[COPIED SOURCE]', basename(src)))
    }
  }
  if (!length(source_scripts)) {
    log <- c(log, '[MISSING SOURCE] No Stage 11 or Stage 12 R files under R_DIR')
  }
}

# Validate propagation of the Stage 15 IPCW sensitivity to Stage 17 snapshot.
ipcw_path <- file.path(WORK_DIR, 'analysis_methodologic_validity', '15_ipcw_results.csv')
snapshot_path <- file.path(WORK_DIR, 'manuscript_ready_final',
                           'analysis_snapshot', 'FINAL_results_snapshot.csv')
if (file.exists(ipcw_path) && file.exists(snapshot_path)) {
  a <- try(read.csv(ipcw_path, check.names = FALSE), silent = TRUE)
  b <- try(read.csv(snapshot_path, check.names = FALSE), silent = TRUE)
  if (!inherits(a, 'try-error') && !inherits(b, 'try-error') &&
      all(c('analysis', 'HR', 'CI_low', 'CI_high', 'p') %in% names(a)) &&
      all(c('Analysis', 'HR', 'CI_low', 'CI_high', 'p') %in% names(b))) {
    x <- a[grepl('1st/99th', a$analysis, fixed = TRUE), , drop=FALSE]
    y <- b[b$Analysis == 'EHR-continuity IPCW; 1st/99th percentile truncation',,
           drop=FALSE]
    if (nrow(x) == 1L && nrow(y) == 1L) {
      vals <- c('HR','CI_low','CI_high','p')
      delta <- abs(as.numeric(x[1, vals]) - as.numeric(y[1, vals]))
      names(delta) <- vals
      ok <- all(is.finite(delta)) && all(delta < 1e-10)
      log <- c(log, paste(if(ok) '[PASS]' else '[FAIL]',
          'Stage 15 -> Stage 17 IPCW propagation;',
          paste(sprintf('%s diff=%.12g', names(delta), delta), collapse='; ')))
    } else {
      log <- c(log, '[REVIEW] Unable to uniquely identify matching IPCW 1st/99th rows')
    }
  } else {
    log <- c(log, '[REVIEW] IPCW comparison schemas did not match expectations')
  }
} else {
  log <- c(log, '[REVIEW] Missing IPCW result and/or Stage 17 snapshot for comparison')
}

# Check that the Stage 17 snapshot retains the unrounded 3-year q-value.
fdr_path <- file.path(WORK_DIR, 'analysis_stage16_high_specificity',
                      '16_expanded_exploratory_FDR.csv')
if (file.exists(fdr_path) && file.exists(snapshot_path)) {
  fdr <- try(read.csv(fdr_path, check.names=FALSE), silent=TRUE)
  ss <- try(read.csv(snapshot_path, check.names=FALSE), silent=TRUE)
  if (!inherits(fdr, 'try-error') && !inherits(ss, 'try-error') &&
      all(c('analysis','q_BH_expanded') %in% names(fdr)) &&
      all(c('Category','q') %in% names(ss))) {
    target <- fdr[grepl('Persistent 3-y.*all 5ARI', fdr$analysis,
                        ignore.case=TRUE),,drop=FALSE]
    if (nrow(target)==1L) {
      q <- target$q_BH_expanded[1]
      found <- any(abs(ss$q - q) < 1e-12, na.rm=TRUE)
      log <- c(log, sprintf('[%s] Three-year q propagated to Stage 17 snapshot: %.12f',
                            if(found) 'PASS' else 'REVIEW', q))
    } else {
      log <- c(log, '[REVIEW] Could not uniquely find three-year all-5ARI comparison in Stage 16 FDR')
    }
  }
}

report <- file.path(DEST, 'COLLECTION_AND_PROPAGATION_REPORT.txt')
writeLines(log, report)
cat(paste(log, collapse='\n'), '\n\n')
cat('Output folder: ', DEST, '\n', sep='')

# Try to build a zip to upload, but keep the folder usable if the zip package is absent.
if (requireNamespace('zip', quietly=TRUE)) {
  output_zip <- file.path(WORK_DIR, 'BPH_final_resident_audit_upload.zip')
  # ZIP only files copied this run, never unrelated files from a prior run.
  files <- unique(c(copied, 'COLLECTION_AND_PROPAGATION_REPORT.txt'))
  files <- files[file.exists(file.path(DEST, files))]
  zip::zipr(zipfile=output_zip, files=files, root=DEST)
  cat('Upload this ZIP: ', output_zip, '\n', sep='')
} else {
  cat('\nThe optional R package "zip" is not installed.\n',
      'Right-click the resident_audit_upload folder and use Send to > Compressed (zipped) folder,\n',
      'then upload that archive. No package installation is necessary.\n', sep='')
}
