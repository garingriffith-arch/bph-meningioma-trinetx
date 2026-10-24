# Reporting-guideline mapping

The principal observational cohort is described with reference to RECORD-PE, together with its RECORD/STROBE foundations. TARGET applies to the separate day-zero initiation-strategy sensitivity analyses; it does not describe the complete one-year-landmark study as a target-trial emulation.

This table points to the corresponding analysis source files. It is not a completed reporting checklist and does not substitute for manuscript page references.

| Reporting topic | Relevant source |
| --- | --- |
| Medication vocabulary and treatment classification | `R/01_extract_medication_timeline.R` |
| Primary cohort eligibility and landmark | `R/02_build_primary_cohort_v3.R`, `R/03_refresh_meningioma_outcomes.R`, `R/05_rebuild_without_future_followup_requirement.R` |
| Revised primary models and overlap weights | `R/06_revised_primary_analysis.R`, `R/final_analysis_helpers.R` |
| Recorded treatment persistence | `R/07_build_core_final_sensitivities.R`, `R/11_build_signal_interrogation_cohorts.R` |
| Exploratory drug-specific and persistence analyses | `R/10_run_final_sensitivity_models.R`, `R/12_run_signal_interrogation_models.R` |
| Day-zero treatment assignment and sensitivity cohorts | `R/14_build_methodologic_validity_cohorts.R` |
| Day-zero and censoring-weighted models | `R/15_run_methodologic_validity_models.R` |
| Landmark-boundary and outcome checks | `R/16_high_specificity_persistence_and_boundary_fix.R` |
| Final tabulation and plots | `R/17_make_final_manuscript_outputs.R` |
| Reporting summaries and missingness | `R/20_RECORD_PE_TARGET_reporting_QC.R` |
| Day-zero participant flow and baseline summaries | `R/21_target_dayzero_reporting.R` |
| Proportional-hazards and rare-outcome reconciliation | `R/24_verify_boundary_fixed_PH.R`, `R/25_reconcile_rare_endpoints.R` |
| Numerical reconciliation and disclosure flags | `R/26_finalize_analysis_audit.R` |
| Study declarations and controlled-access data | `reporting/study-declarations.md`, `docs/data-governance.md` |

The final manuscript and supplement, rather than this mapping, establish which analyses and results were actually reported. Source-code presence alone does not establish reporting-guideline compliance.
