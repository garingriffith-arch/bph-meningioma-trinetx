# Reproducibility and reporting audit

The public repository contains 28 numbered R analysis stages (00–27), a final figure/table builder (Stage 17), source provenance, and numerical reconciliation (Stages 24–26). Historical scripts are deliberately retained and must not be rerun against the locked database. The source release contains no patient-level extracts or local DuckDB databases. Reproducing the results requires separate authorized access to the relevant TriNetX Dataworks export. A passing repository-hygiene test does not constitute independent study reproduction.

## Data coding and missingness
The retained modeling helper converts missing or blank race, ethnicity, marital-status, region, and source indicators to an explicit `Unknown` level. Missing diagnosis-history indicators are constructed as no **recorded** diagnosis, not confirmation that a diagnosis is clinically absent. Propensity-score model frames are checked for missing values, with any excluded rows reported. Audited model fields contained no null/blank entries; explicit unknown demographic categories were retained. No multiple imputation was performed.

## Source geography
`Ex-US` indicates healthcare organizations outside the United States, not necessarily patient residence. These records were included in the locked original cohort, and geography was included among measured confounders. Excluding them would change the target population and requires re-estimation of propensity scores, not simple subtraction of outcomes. The original estimates should not be replaced by a post hoc restriction; any geographically restricted result must be labeled sensitivity analysis.

## Medication records
Clinical medication records can represent prescriptions written or medications administered depending on encounter context. They are not verified dispensing or ingestion records. Persistence criteria describe repeat documentation, not continuous supply or adherence.

## Software environments
The archived model-fitting package versions and an exact historical `renv.lock` were not recovered. Do not attribute the versions currently installed on a workstation to the original fits. Execute `scripts/record_environment.R` on an authorized analysis workstation to save *current* `sessionInfo()` and installed-package versions. This aids future reproduction but does not certify original versions.

## TARGET reporting
The day-zero initiation analysis is separate from the primary one-year landmark. Day-zero arm-level follow-up, events and crude incidence are in supplementary Table S4; cumulative flow and baseline characteristics are in S5–S6. Crude rates are not standardized adjusted absolute risks; adjusted cumulative risks and risk differences were not estimated. End of observable EHR follow-up does not establish a specific cause of censoring. The IPCW result is a sensitivity analysis, not a guarantee against informative observation.
