# Original BPH–meningioma R stages: recovered source code

**Source type: recovered original scripts, not newly reconstructed or rewritten.**

These R files were recovered byte-for-byte from the user's prior September 2026 project archives in the file library. The SHA256_SOURCE_MANIFEST.csv records which archive supplied each file and each delivered file's hash. They have **not** been executed against the local TriNetX data in this environment.

## Requested scripts (copy only if absent)

- `R/00_validate_export.R` — recovered from `bph_m_finasteride_stage1_R.zip`
- `R/01_extract_medication_timeline.R` — same Stage 1 archive
- `R/02_build_primary_cohort_v3.R` — latest recovered Stage 2 v3, fixes a parser error and retains the heavy-scan resume patch
- `R/03_refresh_meningioma_outcomes.R` — recovered Stage 3 correction
- `R/04_primary_analysis_v2.R` — later Stage 4 v2, fixes zero-variance factor failures (replaces original v1)
- `R/05_rebuild_without_future_followup_requirement.R` — revised primary cohort without a future >=365-day observed follow-up requirement
- `R/06_revised_primary_analysis.R` — revised primary overlap-weighted and sensitivity models
- `R/11_build_signal_interrogation_cohorts.R` — healthcare-utilization adjustment and exact-dose homogeneous cohorts
- `R/12_run_signal_interrogation_models.R` — exploratory signal models

**Support file included**: `SUPPORTING_STAGE_07_PERSISTENCE_DEFINITIONS/07_build_core_final_sensitivities.R` (unchanged copy of Stage 07 from the user's previously uploaded `R(3).zip`). This is the file that actually *builds* the one-year sustained and three-year persistent exposure cohorts. Stage 11 **does not** independently define those persistence thresholds; it enriches the pre-existing Stage 07 cohorts with healthcare-utilization measurements.

## Study-rule provenance now directly verifiable

- Stage 01 ingredient dictionary: finasteride and dutasteride (5-ARIs); tamsulosin, alfuzosin, silodosin (uroselective alpha blockers); doxazosin/terazosin (tracked as non-uroselective comparators in source preparation).
- Stage 02: index date >= January 1, 2004; age 40–89; >=365-day baseline observation; BPH diagnosis within [-365, +30] days of index; exclude study-class crossover within days 1–30; exclude known index-dose finasteride 1 mg. **Stage 02 originally imposed an additional follow-up observation restriction superseded by Stage 05.** Exact new-user/washout implementation must be described based on Stage 01 + 02, not invented as a fixed pharmacy-claims washout.
- Stage 05: retains >=1 day of observed risk time after the 365-day outcome landmark, rather than requiring >=365 *additional* observed days; creates `analysis_cohort_primary_revised`. Stage 16 subsequently applies exact-landmark-boundary correction.
- Stage 07 one-year sustained cohort: >=1 repeat same-class medication date on days 90–365; no opposite study class on days 1–365. Strict version additionally requires >=1 same-class date in days 90–180 and another in days 181–365.
- Stage 07 three-year persistence cohort: >=1 distinct same-class medication date in each of days 1–365, 366–730, and 731–1095; no opposite study class on days 1–1095; no broad meningeal tumor before three-year landmark; observed follow-up beyond landmark. This is repeated *EHR medication documentation*, not necessarily confirmed dispensing/adherence.
- Stage 11 measures healthcare utilization *before* the relevant one- or three-year risk-set entry to probe documentation bias; Stage 12 analyzes adjusted drug-specific/persistence cohorts.
- Stage 14 provides separate day-zero target-trial-aligned cohorts. Do not describe the primary one-year landmark cohort as itself a complete day-zero target trial emulation.
- Stage 16 refines the exact-on-landmark tumor convention and is authoritative for final cohort/event counts. The recovered upstream scripts are historical primary source code, not a claim that re-running them reproduces every final boundary-fixed result without the subsequent stages.

## Safe installation and execution guidance

Project code location:
`C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R`

1. Back up the entire `R` directory and `C:/TriNetX/bph_m_work` output tree before installing anything.
2. Copy the nine R files from `R/` into the project code directory **without overwriting any local version that differs**, until hashes/versions have been compared.
3. The included Stage 07 copy is supplied for methodological review only; there is already a Stage 07 script in the later `R(3).zip` collection.
4. **Do not source/rerun these historic stages merely to fill the reporting checklists.** Most scripts execute on `source()` and use `CREATE OR REPLACE TABLE` or overwrite CSV outputs in the same existing DuckDB and analysis directories. Rerunning upstream stages can invalidate the already reconciled Stage 15–20 outputs.
5. First compare source algorithms with final derived-table metadata and amend the RECORD-PE/TARGET Supplementary Methods in light of source evidence. If a from-scratch replication is desired, use a separate working database/output folder and a full, version-locked pipeline.

## Limitations

No R runtime or raw TriNetX extract is available in this execution environment, so this delivery is a **source recovery and static provenance audit**, not a verified full reproduction. Prior script packages document fixes made during development; original source dates and instructions should not be mischaracterized as prospective registration or genuine prespecification. The TriNetX data-sharing restrictions, full RECORD-PE/TARGET checklist, design choice justifications, and Stage 14 day-zero cohort flow remain distinct follow-up issues.
