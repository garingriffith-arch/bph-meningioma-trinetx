# Final analytical and reporting closeout

The validated analyses through Stages 00–23 remain unchanged. These stages close remaining reporting issues using source tables already in `C:/TriNetX/bph_m_work/duckdb/bph_m.duckdb`.

| Stage | Objective | Database access | New inference |
| --- | --- | --- | --- |
| 24 | Refit *previously specified* Stage 16 contrasts to obtain exact-cohort treatment PH tests | Read-only | None; diagnostic tests for prior models |
| 25 | Replay four Stage 10 rare endpoint fits, then refit on canonical N=501,445 cohort | Read-only | Denominator correction of preexisting comparisons |
| 26 | Reconcile final aggregate numerical outputs and produce author-action register | None | None |
| 27 | Inventory or package final Word manuscripts, tables and checklists | None | None |

## Run in Windows PowerShell

Back up the current R scripts and preserve the existing working database. Place `24_...R` through `27_...R` and `run_closeout.R` in the project R folder. Do **not** rerun 00–23.

```powershell
$env:BPH_R_DIR = 'C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R'
$env:BPH_WORK_DIR = 'C:/TriNetX/bph_m_work'
$env:BPH_PH_SCOPE = 'all'  # eight existing Stage 16 models; use 'core' for four key contrasts
$env:BPH_CLOSEOUT_MODE = 'all'
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' "$env:BPH_R_DIR/run_closeout.R"
```

If the R version/path changes, locate `Rscript.exe` as usual. All outputs are outside Git in `C:/TriNetX/bph_m_work/reporting_*`.

To run only the aggregate audit after Stages 24 and 25, set `BPH_CLOSEOUT_MODE=audit`; to run one stage, set `BPH_CLOSEOUT_MODE=24`, `25`, or `26`. The orchestration script skips successfully completed stage-24/25 outputs rather than silently overwriting them.

**Important:** Stage 25 demands exact replay of the original Stage 10 model (N/events/HR) and refuses to accept revised estimates if the replay differs. It also refuses to silently include zero-day follow-up events. If it stops, retain the QC logs and investigate the phenotype/time-at-risk definition rather than overriding the error.

## After Stage 26

1. Inspect `26_closeout_register.csv`, `26_boundary_PH_for_reporting.csv`, `25_rare_boundary_aligned_models_unrounded.csv`, and `25_Table_S11_Rare_endpoints_BOUNDARY_ALIGNED.csv`.
2. Replace Supplementary **Table S11** (the earlier analysis versions called this S7). Do not overwrite the existing Stage 10/17 outputs.
3. Add relevant treatment PH results and IPCW limitations to Supplementary Methods. A nonsignificant Schoenfeld test does not demonstrate PH.
4. The source category `Ex-US` is labeled as data from healthcare organizations outside the United States; do not treat it as a patient-residence indicator. Clarify why it is present in the Dataworks-USA extract if requested.
5. Confirmed: medication orders/prescriptions or administrations depending on context; RFD STUDY00030063 on April 7, 2026; no written protocol; no funding or competing interests. Still required: code disclosure authorization and journal-specific page references.
6. Update RECORD-PE and TARGET checklists with actual final page numbers **after** the manuscript and supplement pagination are frozen. To create a complete source-derived page-map worksheet, run `python scripts/create_checklist_page_map.py RECORD_PE_Checklist_FINAL.docx TARGET_Checklist_FINAL.docx reporting_page_locations.csv` in the `final_documents` folder (adjust the relative script path). Fill `reported_location`, `page_numbers`, and `verified_by` manually. This utility extracts items from the working Word checklists; it does not claim checklist compliance. Final manual page references are not generated from a CSV by guessing Word pagination.
7. Place finalized Word files under `C:/TriNetX/bph_m_work/final_documents/` with these exact names: `Manuscript_FINAL.docx`, `Tables_FINAL.docx`, `Supplementary_Methods_FINAL.docx`, `RECORD_PE_Checklist_FINAL.docx`, and `TARGET_Checklist_FINAL.docx`. Place existing images in `final_documents/figures/`.
8. Run Stage 27 once with `BPH_BUILD_RESIDENT_PACKAGE=0` to inventory documents. If they are reviewed for disclosure, set `BPH_BUILD_RESIDENT_PACKAGE=1` and rerun to package only curated Word/figure files.

Example:

```powershell
$env:BPH_BUILD_RESIDENT_PACKAGE = '1'
& 'C:\Program Files\R\R-4.6.0\bin\Rscript.exe' "$env:BPH_R_DIR/27_prepare_resident_release.R"
```

**No automatic publication, GitHub upload, institutional approval, journal submission, or coauthor email occurs.** Final distribution remains under researcher control.

## Research interpretation

- Primary study: active-comparator, new-user observational cohort with a one-year landmark.
- Day-zero TARGET analyses: separate target-trial-aligned initiation-strategy sensitivity, not a randomized trial.
- Persistence analyses: exploratory associations conditioned on postbaseline medication records; not per-protocol causal effects.
- Ethics: institutional RFD determined this was not human-subjects research. Do not call it IRB exempt or IRB approved unless documentation explicitly supports that.

## Known limitations

The closeout code was reviewed structurally and against previously supplied aggregate schemas, but cannot be executed against the restricted local TriNetX database in this workspace. The four rare endpoint refits may reveal issues requiring revised outcome timing; do not change definitions silently simply to make outputs match.
