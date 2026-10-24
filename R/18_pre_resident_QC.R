# 18_pre_resident_QC.R
# Finasteride/Dutasteride and Meningioma: pre-resident consistency checks.
# Run AFTER 15_run_methodologic_validity_models.R (if re-running IPCW)
# and AFTER 17_make_final_manuscript_outputs.R (rebuild tables with q precision).
# This script does not change any analysis tables or fit inferential models.

options(stringsAsFactors=FALSE, scipen=999)
if (!requireNamespace('data.table', quietly=TRUE)) install.packages('data.table')
if (!requireNamespace('duckdb', quietly=TRUE)) install.packages('duckdb')
if (!requireNamespace('DBI', quietly=TRUE)) install.packages('DBI')
library(data.table)
library(DBI)
library(duckdb)
R_DIR <- 'C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R'
WORK_DIR <- 'C:/TriNetX/bph_m_work'
DB_PATH <- file.path(WORK_DIR,'duckdb','bph_m.duckdb')
OUT_DIR <- file.path(WORK_DIR,'pre_resident_QC')
dir.create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)
model16 <- file.path(WORK_DIR,'analysis_stage16_high_specificity','16_expanded_exploratory_FDR.csv')
final_tables <- file.path(WORK_DIR,'manuscript_ready_final','tables')
final_snapshot <- file.path(WORK_DIR,'manuscript_ready_final','analysis_snapshot','FINAL_results_snapshot.csv')
log <- character()
add <- function(status, label, detail) {
  cat(sprintf('[%s] %s: %s\n',status,label,detail))
  log <<- c(log,sprintf('[%s] %s: %s',status,label,detail))
}
need <- c('R/11 persistence builder' = length(list.files(R_DIR,pattern='^11_.*\\.R$')) > 0,
          'R/12 exploratory models' = length(list.files(R_DIR,pattern='^12_.*\\.R$')) > 0,
          'Stage 16 FDR CSV' = file.exists(model16),
          'Final manuscript snapshot' = file.exists(final_snapshot),
          'Final table S6' = file.exists(file.path(final_tables,'Table_S6_Complete_exploratory_FDR_family.csv')),
          'Final table 3' = file.exists(file.path(final_tables,'Table_3_Exploratory_persistence_results.csv')),
          'DuckDB database' = file.exists(DB_PATH))
for (i in names(need)) add(if(need[[i]]) 'FOUND' else 'MISSING',i,'')
if (file.exists(model16)) {
  fdr <- fread(model16)
  expected_cols <- c('analysis','weighted_p','q_BH_expanded')
  if (!all(expected_cols %in% names(fdr))) stop('Stage 16 FDR input schema mismatch; verify output version.')
  ok <- nrow(fdr)==16L && !anyDuplicated(fdr$analysis) &&
     all(is.finite(fdr$weighted_p)) && all(is.finite(fdr$q_BH_expanded))
  add(if(ok) 'PASS' else 'FAIL','FDR family',paste('number of rows:',nrow(fdr),'distinct:',uniqueN(fdr$analysis)))
  if (!ok) stop('Invalid Stage 16 FDR family')
  expected <- p.adjust(fdr$weighted_p,method='BH')
  difference <- max(abs(expected-fdr$q_BH_expanded))
  add(if(difference < 1e-9) 'PASS' else 'FAIL','BH q recomputation',paste('largest abs diff:',format(difference,digits=6)))
  if(difference >= 1e-9) stop('Saved Stage 16 FDR does not match p.adjust.')
  borderline <- fdr[abs(q_BH_expanded-0.05)<0.01,.(analysis,weighted_p,q_BH_expanded,
                significant_unrounded = q_BH_expanded < 0.05)]
  fwrite(borderline,file.path(OUT_DIR,'borderline_FDR_unrounded.csv'))
  if(nrow(borderline)) {
    for(i in seq_len(nrow(borderline))) {
      x <- borderline[i]
      add('REVIEW','Borderline q',sprintf('%s: q=%.12f; q<0.05=%s',x$analysis,x$q_BH_expanded,x$significant_unrounded))
    }
  }
}
if (file.exists(DB_PATH)) {
  con <- dbConnect(duckdb::duckdb(),dbdir=DB_PATH,read_only=TRUE)
  tables <- dbListTables(con)
  if ('analysis_cohort_primary_boundary_fixed' %in% tables) {
    counts <- dbGetQuery(con,"SELECT COUNT(*) n, SUM(CASE WHEN exposure='5ARI' THEN 1 ELSE 0 END) n_5ari, SUM(CASE WHEN exposure<>'5ARI' THEN 1 ELSE 0 END) n_alpha, SUM(meningioma_event_strict_boundary) n_events FROM analysis_cohort_primary_boundary_fixed")
    fwrite(as.data.table(counts),file.path(OUT_DIR,'boundary_fixed_cohort_counts.csv'))
    expected <- c(n=501445,n_5ari=51390,n_alpha=450055,n_events=1456)
    observed <- as.numeric(counts[1,names(expected)])
    ok <- all(observed==unname(expected))
    add(if(ok) 'PASS' else 'FAIL','Manuscript primary counts',paste(names(expected),observed,sep='=',collapse='; '))
  } else add('MISSING','Boundary-corrected primary table','Run Stage 16')
  dbDisconnect(con,shutdown=TRUE)
}
# Different results populations are expected by analysis design. Verify S7 still uses
# the 501446-patient population; never silently substitute 501445 without refitting.
s7 <- file.path(final_tables,'Table_S7_Rare_hard_endpoints.csv')
if (file.exists(s7)) {
  a <- fread(s7)
  if('N' %in% names(a)) add('REVIEW','S7 model denominators',paste(unique(a$N),collapse=', '))
}
# Workflow/documentation checks requiring source files or institutional action.
add('AUTHOR','Three-year persistence definition','Check Stage 11 exact medication windows, switching and treatment grace periods; reconcile manuscript')
add('AUTHOR','Prespecification','Verify dates for falsification outcomes and 16-model family against protocol/revision history')
add('AUTHOR','Ethics and data sharing','Fill IRB regulatory determination and Dataworks/contract permissions; do not assume')
add('AUTHOR','Geographic covariate','Check region Ex-US against Dataworks source-field dictionary')
add('AUTHOR','IPCW ordering','If Stage 15 was rerun with patched ordering, compare results and censoring diagnostics to archived run')
add('AUTHOR','Rare-event analyses','Stage 10 estimates have original N=501446; justify one-patient difference or refit faithfully')
writeLines(c('BPH/MENINGIOMA PRE-RESIDENT QUALITY CONTROL','',log),file.path(OUT_DIR,'QC_REPORT.txt'))
cat('\nReport at: ',file.path(OUT_DIR,'QC_REPORT.txt'),'\n')
