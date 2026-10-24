# Supplementary Methods: source-code-verified definitions

Status: scientific draft for author verification. Definitions below were recovered from the original September study scripts and were checked against later Stage 14–17 boundary corrections. Numerical results remain those previously reconciled in Stage 19/20.

## Medication exposure and comparator selection

Stage 01 identifies medication records using the exported `RxNorm` `code_system` and the ingredient vocabulary defined below. Source code (ingredient-level mapping) should not be represented as verified dispensing or patient ingestion.

| Drug | RxNorm ingredient code | Study classification |
| --- | --- | --- |
| Finasteride | 25025 | 5-ARI |
| Dutasteride | 228790 | 5-ARI |
| Tamsulosin | 77492 | Uroselective α-blocker |
| Alfuzosin | 17300 | Uroselective α-blocker |
| Silodosin | 720825 | Uroselective α-blocker |
| Doxazosin | 49276 | Other α-blocker (not study comparator) |
| Terazosin | 37798 | Other α-blocker (not study comparator) |

Stage 01 collapses multiple records to a unique patient–drug–calendar-day event for initiation timing. `route`, `brand`, `strength` and `source_id` are retained as representative metadata; they do not establish dispensing. Stage 02 excludes recorded finasteride 1-mg starts to reduce an alopecia-indication component. This is **not** proof that all finasteride 1-mg medication use was captured.

## Primary cohort (historical and revised)

Stage 02's index begins at the first observed class start after a minimum of 365 days of prior possible EHR observation with at least one recorded pre-index encounter. Eligible men were aged 40–89 years, with treatment starts from January 1, 2004. BPH could be documented from 365 days before index to 30 days after index; baseline diagnosis covariates use the previous 1,095 days. Opposite-class starts during a 30-day stabilization period influenced the primary exposure assignment/exclusions. The primary outcome enters the risk set at a 365-day landmark, with exclusion of previously recorded broad meningeal neoplasms. Stages 05–06 remove a problematic future-follow-up eligibility condition from the original cohort construction. Stage 16 corrects outcomes at **exactly** the landmark: tumors on or before entry are excluded and events strictly after entry are counted.

**Interpretation:** “new user” means first observed study-class medication start within the available exported history under the code-defined eligibility filters; it must not be equated automatically with a verified pharmacy dispensing washout. The primary cohort is not a day-zero target-trial emulation.

## Sustained 1-year exposure and 3-year persistence (Stage 07)

Stage 07 builds patient-specific same-class and opposite-class distinct medication-start dates from `study_med_events`, relative to the original treatment index date:

- **One-year sustained:** at least one further same-class date from day 90 through day 365 after index, and no opposite study-class date from day 1 through day 365; eligible patients come from the revised one-year landmark cohort.
- **Strict one-year sustained:** satisfies the above plus at least one same-class date from day 90 to before day 181 and another from day 181 through day 365.
- **Three-year persistence:** at least one same-class date in *each* of days 1–365, 366–730, and 731–1095, no opposite study-class date in days 1–1095, no recorded broad meningeal neoplasm before the three-year landmark, and observation extending beyond the landmark. Stage 16 subsequently enforces the corrected exact-landmark outcome boundary for the persistence sensitivity models.

A medication event is documentation in the EHR export, **not proof of a filled prescription, continuous supply, adherence, or exposure throughout the entire interval**. Selection into landmark-defined persistent groups may depend on healthcare utilization and other postbaseline characteristics. Persistence comparisons are exploratory associations, not effects of sustained treatment under a fully specified dynamic treatment strategy.

## Utilization and exploratory analyses (Stages 11–12)

Stage 11 augments Stage 07 sustained/persistent risk sets with pre-landmark encounters and constructs drug-specific/modern-era cohorts; Stage 12 estimates overlap-weighted exploratory models. The Stage 16 expanded multiple-testing family includes 16 distinct comparisons; the threshold-adjacent three-year class-level result has an unrounded Benjamini–Hochberg-adjusted *q* of 0.0495876444 (reported as 0.0496), whereas the one-year sustained finasteride comparison has *q* of approximately 0.0591. These labels and values were verified from the Stage 18–19 aggregate audits, not by refitting individual-level models in this repository.

## Day-zero target-trial-aligned sensitivity (Stages 14–15)

Stage 14 assigns treatment using the medication class documented **on index day** and does not use later switching or days 1–30 medication records to revise assignment. It requires BPH documented in the 365 days before or on index, excludes simultaneous day-zero 5-ARI and uroselective α-blocker starts, and excludes meningeal tumors present before or on day zero. The zero-lag analysis follows outcomes from index; a separate one-year-latency sensitivity applies a later outcome landmark. Stage 15 fits observational overlap-weighted models; these are *intention-to-treat-like* initiator contrasts, not randomized intention-to-treat effects.

The day-zero sample (660,012 patients in prior reporting) differs from the one-year primary cohort (501,445). The day-zero flow and baseline cannot be borrowed from the primary flow or Table 1. Stage 21 generates these summaries directly from the finalized Stage 14 tables with a cumulative eligibility-reconciliation check.

## Overlap weighting and model limitations

Treatment assignment was modeled from measured baseline variables with flexible age and calendar-time terms. Overlap weighting targets the treatment-overlap population, rather than the population-average treatment effect for all men with BPH. IPCW was used in sensitivity analyses to assess censoring related to EHR continuity; extreme untruncated weights and imperfect EHR observation remain limitations. The negative-control benign cranial nerve neoplasm association (HR 1.34; 95% CI 1.06–1.69) raises concern for ascertainment or residual confounding and should not be dismissed. No model fit can establish absence of unmeasured confounding.

## Outcome codes (Stage 14)

- BPH: ICD-10-CM `N40*`; ICD-9-CM `600*`.
- Prostate cancer exclusion: ICD-10-CM `C61`; ICD-9-CM `185`.
- Primary strict benign meningioma: ICD-10-CM `D32.0`, `D32.9`; ICD-9-CM `225.2`.
- High-specificity intracranial meningioma sensitivity: ICD-10-CM `D32.0`; ICD-9-CM `225.2`.
- Falsification outcomes: benign cranial nerve neoplasm (`D33.3`, `225.1`) and benign spinal cord neoplasm (`D33.4`, `225.3`).

These are algorithm definitions, not externally adjudicated phenotype validation. Use the full code list and details of source records in the final supplement.

## Ethics and data access

The study underwent an institutional Research Feasibility Determination and was determined **not to constitute human subjects research**; institutional review board review was therefore not required. RFD STUDY00030063 was determined April 7, 2026. Raw data and underlying EHR extracts remain subject to TriNetX/institutional access terms. The GitHub repository must remain private pending authorization for code release.

## Remaining sign-off before manuscript submission

1. Independently verify the source `medication.csv` record semantics (order, history, administration, or dispensing) from its TriNetX extract dictionary; code alone does not settle provenance.
2. Confirm **re-entry**, episode construction and prior-exposure wording from complete Stage 02/05 cohort logs and export documentation.
3. Review the PH diagnostics file and clarify any departures from proportional hazards; no blanket claim of satisfied assumptions.
4. Explain the Ex-US outside-US healthcare-organization category and its occurrence in the Dataworks-USA-labeled extract.
5. Refit or clearly annotate the rare-endpoint table whose historical denominator is 501,446 rather than the 501,445 boundary-fixed cohort.
6. Run Stage 21 and verify its day-zero eligibility flow and baseline table.
7. No written study protocol was used; do not claim prospective prespecification of exploratory analyses.
8. Confirm final sharing/registration/ethics statements against institutional documents.


## Confirmed source and regulatory details (October 2026)

RFD STUDY00030063 (April 7, 2026) determined non-human-subjects research; no IRB oversight. No written protocol was prepared. No funding, conflicts, or additional disclosures were reported. The source Ex-US category denotes outside-US healthcare organizations. TriNetX EHR medication records may reflect prescriptions written or administrations depending on context; dispensing and adherence are not verified. Final submission checklist page numbers remain dependent on target-journal layout.
