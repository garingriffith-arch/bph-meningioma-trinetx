# Record package versions currently installed on an authorized workstation.
# A current snapshot cannot establish the versions used for historical model fits.

work_dir <- Sys.getenv("BPH_WORK_DIR", "C:/TriNetX/bph_m_work")
out_dir <- Sys.getenv("BPH_ENV_REPORT_DIR",
  file.path(work_dir, "reporting_environment_snapshot"))
if (dir.exists(out_dir)) stop("Archive or move existing folder: ", out_dir)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

packages <- c("DBI","duckdb","data.table","survival","splines",
              "ggplot2","scales","flextable","officer","speedglm","zip")
installed <- vapply(packages, requireNamespace, logical(1), quietly=TRUE)
versions <- vapply(packages, function(p) {
  if (!requireNamespace(p, quietly=TRUE)) return(NA_character_)
  as.character(utils::packageVersion(p))
}, character(1))
utils::write.csv(data.frame(package=packages,installed=installed,version=versions),
  file.path(out_dir,"installed_package_versions.csv"), row.names=FALSE, na="")
writeLines(c("Current workstation snapshot; NOT evidence of original model execution",
  paste("Captured:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("Platform:", R.version$platform),
  "", capture.output(utils::sessionInfo())),
  file.path(out_dir,"sessionInfo.txt"))
message("Current environment inventory: ",normalizePath(out_dir,winslash="/"))
