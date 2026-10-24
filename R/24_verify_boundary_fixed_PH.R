# 24_verify_boundary_fixed_PH.R
# Exact-cohort proportional-hazards verification for archived Stage 16 fits.
# Read-only DuckDB connection. No patient-level data or database tables are written.
# Stage 16's model specifications are preserved; only aggregate PH test results are exported.
#
# Default: four decision-critical Stage 16 models.
#   BPH_PH_SCOPE=core (default): primary, target-trial 1-y landmark,
#     1-y strict endpoint/finasteride, 3-y strict endpoint/finasteride.
#   BPH_PH_SCOPE=all: eight Stage 16 models (including high-specificity outcomes).
#
# Prerequisites: final_analysis_helpers.R, completed Stage 16 data tables and
# analysis_stage16_high_specificity/16_model_results.csv. Do NOT rerun earlier stages.

run_stage_24 <- function() {
options(stringsAsFactors = FALSE, scipen = 999)
for (pkg in c('DBI','duckdb','data.table','survival','splines')) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop('Install required package: ',pkg)
}
suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(data.table)
  library(survival); library(splines)
})
R_DIR <- Sys.getenv('BPH_R_DIR',
 'C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R')
WORK_DIR <- Sys.getenv('BPH_WORK_DIR', 'C:/TriNetX/bph_m_work')
DB_PATH <- file.path(WORK_DIR, 'duckdb', 'bph_m.duckdb')
OUT_DIR <- file.path(WORK_DIR, 'reporting_boundary_PH_stage24')
SCOPE <- tolower(Sys.getenv('BPH_PH_SCOPE', 'core'))
if (!SCOPE %in% c('core','all')) stop('BPH_PH_SCOPE must be core or all')
if (!file.exists(DB_PATH)) stop('Missing DuckDB: ',DB_PATH)
if (!file.exists(file.path(R_DIR,'final_analysis_helpers.R'))) stop('Missing helper file in ',R_DIR)
models_file <- file.path(WORK_DIR,'analysis_stage16_high_specificity','16_model_results.csv')
if (!file.exists(models_file)) stop('Missing Stage 16 reference file: ',models_file)
if (dir.exists(OUT_DIR)) stop('Output folder already exists. Archive or move: ',OUT_DIR)
dir.create(OUT_DIR,recursive=TRUE)
source(file.path(R_DIR, 'final_analysis_helpers.R'))
reference <- fread(models_file)
con <- dbConnect(duckdb::duckdb(),dbdir=DB_PATH,read_only=TRUE)
on.exit(try(dbDisconnect(con,shutdown=TRUE), silent=TRUE),add=TRUE)

base_cols <- c('patient_id','exposure','index_drug',
 'index_finasteride_strength','index_dutasteride_strength','index_date',
 'age_at_index','race','ethnicity','marital_status','patient_regional_location',
 'patient_source_id','hypertension','diabetes','obesity','ckd','smoking','cad',
 'heart_failure','cerebrovascular_disease','hyperlipidemia','osteoporosis',
 'alopecia','urinary_retention','uti','hematuria','bladder_stone',
 'n_encounters_prior_365','n_encounters_prior_1095','n_ambulatory_prior_365',
 'n_ed_prior_365','n_inpatient_prior_365')
util1 <- c('lm_encounters_total','lm_encounter_days','lm_ambulatory','lm_ed',
           'lm_inpatient','lm_y1_encounter_days')
util3 <- c(util1,'lm_y2_encounter_days','lm_y3_encounter_days')
extra1 <- paste0('log_',c('lm_encounters_total','lm_encounter_days',
                             'lm_ambulatory','lm_ed','lm_inpatient'))
extra3 <- paste0('log_',c('lm_encounters_total','lm_ambulatory','lm_ed',
                          'lm_inpatient','lm_y1_encounter_days','lm_y2_encounter_days',
                          'lm_y3_encounter_days'))
required_tables <- c('analysis_cohort_primary_boundary_fixed',
 'analysis_cohort_tte_lag1_boundary_fixed',
 'analysis_cohort_sustained_1y_strict_boundary_fixed',
 'analysis_cohort_persistent_3y_boundary_fixed')
if (!all(required_tables %in% dbListTables(con))) {
 stop('Missing Stage 16 boundary-fixed tables: ',
      paste(setdiff(required_tables,dbListTables(con)),collapse=', '))
}
read_cohort <- function(table, extras=character()) {
 existing <- dbListFields(con,table)
 requested <- unique(c(base_cols,
    'meningioma_event_strict_boundary','followup_days_strict_boundary',
    'meningioma_event_specific_boundary','followup_days_specific_boundary',
    extras))
 use <- requested[requested %in% existing]
 if (!all(c('exposure','meningioma_event_strict_boundary',
            'followup_days_strict_boundary') %in% use)) stop('Missing required fields in ',table)
 q <- paste0('SELECT ',paste(sprintf('"%s"',use),collapse=','),
             ' FROM "',table,'"')
 d <- as.data.table(dbGetQuery(con,q))
 d[,trt:=as.integer(exposure=='5ARI')]
 for (nm in extras) if (nm %in% names(d))
   d[,(paste0('log_',nm)):=log1p(as.numeric(get(nm)))]
 d
}

run_spec <- function(data, label, event, time, extra=character()) {
 message('Validating Stage 16 model: ',label)
 res <- run_ato_generic(data,label,event,time,extra_numeric=extra)
 rr <- reference[analysis==label]
 if (nrow(rr)!=1L) stop('Expected exactly one matching Stage 16 reference row: ',label)
 result <- res$result
 if (nrow(result)!=1L) stop('Unexpected result shape: ',label)
 for (nm in c('n','events')) {
   if (!identical(as.integer(result[[nm]][1]),as.integer(rr[[nm]][1])))
     stop('Count mismatch for ',label,': ',nm)
 }
 if (!is.finite(result$weighted_HR[1]) ||
     abs(result$weighted_HR[1]-rr$weighted_HR[1]) > 1e-5)
   stop('Weighted HR differs from archived Stage 16 by >1e-5 for ',label)
 if (is.null(res$zph) || nrow(res$zph)==0L) stop('PH diagnostic unavailable: ',label)
 tab <- as.data.table(res$zph)[term=='trt',.(analysis,chisq,df,p,term)]
 if (nrow(tab)!=1L || !is.finite(tab$p[1])) stop('No valid treatment PH test: ',label)
 tab[,`:=`(n=result$n[1],events=result$events[1],
          weighted_HR=result$weighted_HR[1],
          reference_weighted_HR=rr$weighted_HR[1],
          matched_stage16=TRUE)]
 tab
}
reslist<-list()
primary <- read_cohort('analysis_cohort_primary_boundary_fixed')
reslist[[length(reslist)+1]] <- run_spec(primary,
 'Boundary-fixed primary 1-y landmark',
 'meningioma_event_strict_boundary','followup_days_strict_boundary')
rm(primary);gc(verbose=FALSE)
tte <- read_cohort('analysis_cohort_tte_lag1_boundary_fixed')
reslist[[length(reslist)+1]] <- run_spec(tte,
 'Boundary-fixed target-trial day-0 assignment, 1-y landmark',
 'meningioma_event_strict_boundary','followup_days_strict_boundary')
rm(tte);gc(verbose=FALSE)
s1 <- read_cohort('analysis_cohort_sustained_1y_strict_boundary_fixed',util1)
s1fin <- s1[exposure=='URO_ALPHA' | (exposure=='5ARI' & index_drug=='finasteride')]
reslist[[length(reslist)+1]] <- run_spec(s1fin,
 'Boundary-fixed strict endpoint: strict sustained 1-y, finasteride vs uro-alpha',
 'meningioma_event_strict_boundary','followup_days_strict_boundary',extra1)
if (SCOPE=='all') {
 reslist[[length(reslist)+1]] <- run_spec(s1,
  'High-specificity: strict sustained 1-y, all 5ARI vs uro-alpha',
  'meningioma_event_specific_boundary','followup_days_specific_boundary',extra1)
 reslist[[length(reslist)+1]] <- run_spec(s1fin,
  'High-specificity: strict sustained 1-y, finasteride vs uro-alpha',
  'meningioma_event_specific_boundary','followup_days_specific_boundary',extra1)
}
rm(s1,s1fin);gc(verbose=FALSE)
p3 <- read_cohort('analysis_cohort_persistent_3y_boundary_fixed',util3)
p3fin <- p3[exposure=='URO_ALPHA' | (exposure=='5ARI' & index_drug=='finasteride')]
reslist[[length(reslist)+1]] <- run_spec(p3fin,
 'Boundary-fixed strict endpoint: persistent 3-y, finasteride vs uro-alpha',
 'meningioma_event_strict_boundary','followup_days_strict_boundary',extra3)
if (SCOPE=='all') {
 reslist[[length(reslist)+1]] <- run_spec(p3,
  'High-specificity: persistent 3-y, all 5ARI vs uro-alpha',
  'meningioma_event_specific_boundary','followup_days_specific_boundary',extra3)
 reslist[[length(reslist)+1]] <- run_spec(p3fin,
  'High-specificity: persistent 3-y, finasteride vs uro-alpha',
  'meningioma_event_specific_boundary','followup_days_specific_boundary',extra3)
}
rm(p3,p3fin);gc(verbose=FALSE)
out <- rbindlist(reslist,fill=TRUE)
if (SCOPE=='core' && nrow(out)!=4L) stop('Expected four core models')
if (SCOPE=='all' && nrow(out)!=8L) stop('Expected eight models')
fwrite(out,file.path(OUT_DIR,'24_boundary_fixed_PH_diagnostics.csv'))
writeLines(c(
 paste('Stage 24 PH model scope:',SCOPE),
 'Database access: read-only',
 'Inference: refit of existing Stage 16 models to obtain model-specific PH tests only',
 'All model n/events/weighted HR reconciled with Stage 16 saved aggregate estimates',
 'Do not interpret nonsignificant PH tests as proof of proportional hazards.'),
 file.path(OUT_DIR,'24_README.txt'))
message('Stage 24 completed: ',file.path(OUT_DIR,'24_boundary_fixed_PH_diagnostics.csv'))

}
run_stage_24()
