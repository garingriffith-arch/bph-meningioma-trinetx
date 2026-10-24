# Stage 27: Assemble document-only resident review packet from finalized files.
# No database access or statistical modeling. Does not write into GitHub.
# Default mode inventories expected documents. To copy/zip, set
# BPH_BUILD_RESIDENT_PACKAGE=1 after verifying Word files and sharing permissions.

run_stage_27 <- function() {
  if (!requireNamespace("data.table", quietly=TRUE)) stop("Install data.table")
  suppressPackageStartupMessages(library(data.table))
  WORK_DIR <- Sys.getenv("BPH_WORK_DIR", "C:/TriNetX/bph_m_work")
  DOC_DIR <- Sys.getenv("BPH_FINAL_DOC_DIR", file.path(WORK_DIR,"final_documents"))
  OUT_DIR <- file.path(WORK_DIR,"reporting_resident_release_stage27")
  dir.create(OUT_DIR, recursive=TRUE,showWarnings=FALSE)
  expected <- data.table(
    role=c("manuscript", "main_and_supplemental_tables", "supplementary_methods",
      "record_pe_checklist", "target_checklist"),
    filename=c("Manuscript_FINAL.docx", "Tables_FINAL.docx", "Supplementary_Methods_FINAL.docx",
      "RECORD_PE_Checklist_FINAL.docx", "TARGET_Checklist_FINAL.docx")
  )
  expected[, path := file.path(DOC_DIR,filename)]
  expected[, exists := file.exists(path)]
  fwrite(expected, file.path(OUT_DIR,"27_document_inventory.csv"))

  closeout_path <- file.path(WORK_DIR,"reporting_final_closeout_stage26","26_closeout_register.csv")
  if (!file.exists(closeout_path)) stop("Run Stage 26 before making a release packet")
  audit <- fread(closeout_path)
  aggregate_errors <- audit[status %in% c("MISSING","REVIEW"), .(item,status,action)]
  fwrite(aggregate_errors, file.path(OUT_DIR,"27_remaining_analytical_issues.csv"))
  if (nrow(aggregate_errors)) message("Some analytical/reporting items remain open; inspect 27_remaining_analytical_issues.csv")

  # Optional CSV for manual pages; these are reporting locations, not declarations
  # of guideline compliance. Populate final Word page refs after layout is locked.
  pages_path <- file.path(DOC_DIR,"reporting_page_locations.csv")
  if (file.exists(pages_path)) {
    pages <- fread(pages_path)
    expected_cols <- c("guideline", "item", "reported_location", "page_numbers", "verified_by")
    if (!all(expected_cols %in% names(pages))) stop("Page map schema missing fields")
    if (anyNA(pages$page_numbers) || any(!nzchar(trimws(pages$page_numbers))))
      message("Page mapping incomplete; do not call the reporting checklists final")
    fwrite(pages, file.path(OUT_DIR,"27_reporting_page_locations.csv"))
  } else {
    fwrite(data.table(guideline=c("RECORD-PE","TARGET"),item=c("Complete checklist","Applicable day-zero items"),
      reported_location="",page_numbers="",verified_by=""),
      file.path(OUT_DIR,"27_reporting_page_locations_TEMPLATE.csv"))
  }

  build <- identical(Sys.getenv("BPH_BUILD_RESIDENT_PACKAGE", "0"), "1")
  if (!build) {
    writeLines(c("Stage 27 inventory only. No release ZIP generated.",
      "Review 26_closeout_register.csv; place final Word documents in:",DOC_DIR,
      "Then set BPH_BUILD_RESIDENT_PACKAGE=1 to create a resident review ZIP.",
      "No raw exports, patient-level data or restricted RFD letter will be copied."),
      file.path(OUT_DIR,"27_README.txt"))
    message("Stage 27 inventory saved: ", OUT_DIR)
    return(invisible(NULL))
  }
  if (!all(expected$exists)) stop("Required Word files missing: ", paste(expected$filename[!expected$exists], collapse=", "))
  if (nrow(aggregate_errors)) warning("Release requested despite unresolved analytical items; see inventory")
  dest <- file.path(OUT_DIR,"resident_review_documents")
  dir.create(dest, recursive=TRUE, showWarnings=FALSE)
  for (i in seq_len(nrow(expected))) {
    ok <- file.copy(expected$path[i],file.path(dest,expected$filename[i]),overwrite=TRUE)
    if (!ok) stop("Failed copying: ",expected$path[i])
  }
  figs <- file.path(DOC_DIR,"figures")
  allowed <- c("png","tif","tiff","pdf")
  if (dir.exists(figs)) {
    image_files <- list.files(figs, full.names=TRUE,recursive=FALSE)
    image_files <- image_files[tolower(tools::file_ext(image_files)) %in% allowed]
    if (length(image_files)) {
      dir.create(file.path(dest,"figures"),showWarnings=FALSE)
      file.copy(image_files,file.path(dest,"figures"),overwrite=TRUE)
    }
  }
  writeLines(c("Draft for resident coauthor review only.",
    "Primary analysis is observational and uses a one-year landmark.",
    "Exploratory medication-persistence analyses are noncausal associations.",
    "Ethics: institutional RFD determined non-human-subjects research.",
    "This package contains curated documents and figures only; no patient-level files."),
    file.path(dest,"README_FOR_REVIEWERS.txt"))
  zip_file <- file.path(OUT_DIR,"BPH_Meningioma_Resident_Review_Documents.zip")
  oldwd <- getwd(); on.exit(setwd(oldwd),add=TRUE)
  setwd(OUT_DIR)
  files <- list.files("resident_review_documents",recursive=TRUE,full.names=TRUE)
  if (!requireNamespace("zip",quietly=TRUE)) stop("Install package 'zip' to write ZIP; collected files remain in ",dest)
  zip::zipr(zipfile=zip_file,files=files,root=OUT_DIR)
  message("Created document-only review ZIP: ",zip_file)
}
run_stage_27()
