# Post hoc Ex-US source-geography sensitivity for the locked primary model.
# Read-only: neither the source cohort nor database tables are altered.
# No patient-level rows are exported.
#
# Run only on an authorized workstation with the completed Stage 16 DuckDB.
# BPH_R_DIR points to the local final_analysis_helpers.R directory.
# BPH_WORK_DIR points to the validated analysis workspace.

for (p in c("DBI","duckdb","data.table","survival","splines")) {
  if (!requireNamespace(p,quietly=TRUE)) stop("Missing required package: ",p)
}
suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(data.table)
  library(survival); library(splines)
})
R_DIR <- Sys.getenv("BPH_R_DIR",
 "C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R")
WORK_DIR <- Sys.getenv("BPH_WORK_DIR","C:/TriNetX/bph_m_work")
db_path <- file.path(WORK_DIR,"duckdb","bph_m.duckdb")
ref_path <- file.path(WORK_DIR,"analysis_stage16_high_specificity","16_model_results.csv")
out_dir <- file.path(WORK_DIR,"reporting_geography_sensitivity")
helper <- file.path(R_DIR,"final_analysis_helpers.R")
for (f in c(db_path,ref_path,helper)) if (!file.exists(f)) stop("Missing: ",f)
if (dir.exists(out_dir)) stop("Move or archive the existing folder: ",out_dir)
source(helper)
reference <- data.table::fread(ref_path)
label <- "Boundary-fixed primary 1-y landmark"
ref <- reference[analysis==label]
if (nrow(ref)!=1L) stop("Archived Stage 16 primary result unavailable")
con <- DBI::dbConnect(duckdb::duckdb(),dbdir=db_path,read_only=TRUE)
on.exit(try(DBI::dbDisconnect(con,shutdown=TRUE),silent=TRUE),add=TRUE)
tbl <- "analysis_cohort_primary_boundary_fixed"
if (!tbl %in% DBI::dbListTables(con)) stop("Final primary cohort table missing")
cols <- c("patient_id","exposure","index_date","age_at_index",
 "race","ethnicity","marital_status","patient_regional_location",
 "patient_source_id","hypertension","diabetes","obesity","ckd",
 "smoking","cad","heart_failure","cerebrovascular_disease",
 "hyperlipidemia","osteoporosis","alopecia","urinary_retention",
 "uti","hematuria","bladder_stone","n_encounters_prior_365",
 "n_encounters_prior_1095","n_ambulatory_prior_365","n_ed_prior_365",
 "n_inpatient_prior_365","meningioma_event_strict_boundary",
 "followup_days_strict_boundary")
if (!all(cols %in% DBI::dbListFields(con,tbl)))
  stop("Required fields missing: ",paste(setdiff(cols,DBI::dbListFields(con,tbl)),collapse=", "))
sql <- paste0("SELECT ",paste(sprintf('"%s"',cols),collapse=","),
              ' FROM "',tbl,'"')
d <- data.table::as.data.table(DBI::dbGetQuery(con,sql))
if (data.table::uniqueN(d$patient_id)!=nrow(d)) stop("Duplicate patients")
d[,trt:=as.integer(exposure=="5ARI")]
if (!all(sort(unique(d$trt))==c(0L,1L))) stop("Invalid treatment coding")
full <- run_ato_generic(d,label,
 "meningioma_event_strict_boundary","followup_days_strict_boundary")$result
if (nrow(full)!=1L ||
    as.integer(full$n)!=as.integer(ref$n) ||
    as.integer(full$events)!=as.integer(ref$events) ||
    abs(full$weighted_HR-ref$weighted_HR)>1e-5)
  stop("Archived full-cohort result did not reproduce; stop")

region <- toupper(trimws(as.character(d$patient_regional_location)))
ex_us <- !is.na(region) & region %in%
  c("EX-US","EX US","EXCEPT UNITED STATES")
if (!any(ex_us)) stop("Ex-US coding not recognized; check source codebook")
summary_by_arm <- d[,.(n=.N,events=sum(meningioma_event_strict_boundary),
  ex_us=sum(ex_us[.I]),unknown_or_blank=sum(is.na(region[.I]) |
    region[.I] %in% c("","UNKNOWN"))),by=exposure]
summary_by_arm[,ex_us_pct:=round(100*ex_us/n,2)]
restricted <- d[!ex_us]
if (!all(c("5ARI","URO_ALPHA") %in% restricted$exposure))
  stop("Treatment group absent in restricted cohort")
sensitivity <- run_ato_generic(restricted,
 "Post hoc Ex-US-excluded HCO source analysis",
 "meningioma_event_strict_boundary","followup_days_strict_boundary")$result
results <- data.table::rbindlist(list(full,sensitivity),fill=TRUE)
dir.create(out_dir,recursive=TRUE)
data.table::fwrite(results,file.path(out_dir,"geography_sensitivity_results.csv"))
data.table::fwrite(summary_by_arm,file.path(out_dir,"source_region_counts_by_arm.csv"))
writeLines(c(
 "Ex-US indicates HCO source geography and is not a patient-residence measure.",
 "The full-cohort model is first reproduced against the archived Stage 16 result.",
 "The second analysis excludes recorded Ex-US sources only; unknowns remain.",
 "Propensity scores and overlap weights were re-estimated in the restricted sample.",
 "Post hoc sensitivity; primary study results and FDR family remain unchanged."),
 file.path(out_dir,"README.txt"))
message("Read-only aggregate report written to: ",out_dir)
