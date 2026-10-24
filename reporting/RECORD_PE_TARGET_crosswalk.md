# RECORD-PE / TARGET reporting crosswalk

**Scope**: RECORD-PE applies to the active-comparator observational study, including RECORD and STROBE foundations. TARGET applies to the *day-zero target-trial-aligned sensitivity analyses* only. The checklist must not classify the full one-year-landmark analysis as a completed target-trial emulation.

| Reporting requirement | Current evidence and manuscript location | Status/action |
| --- | --- | --- |
| RECORD-PE: treatment vocabulary/classification | Stage 01 RxNorm dictionary; Supplementary Methods | **Source verified**; validate extract semantics |
| RECORD-PE: study design and time windows | Stages 02, 05, 07, 14; Figure S1 design schematic | **Source verified**; confirm displayed figure uses exact dates |
| RECORD-PE: treatment switching and sustained exposure | Stage 07 medication-date counts; Stages 11–12 | **Source verified**; label as recorded persistence, not filled adherence |
| RECORD-PE: BPH/comparator/new-user rule | Stages 01–02; original cohort flow | **Partly verified**; prior-use interpretation remains limited by captured history |
| RECORD-PE: proportional hazards, overlap and IPCW | Stages 06, 15; diagnostics CSVs | **Partly verified**; insert model-specific PH values |
| RECORD-PE: outcome algorithms and controls | Stage 14 codes; Table S4/S5 | **Source verified**; avoid claiming chart-adjudicated validation |
| STROBE/RECORD: participants, exclusions | Figure 1 and primary Stage 02/05 flow | **Present**; verify revised and final boundary denominator reconciliation |
| STROBE/RECORD: follow-up, missingness, incidence | Stage 20 Table S8, missingness report | **Present in internal draft**; author/data disclosure review |
| RECORD: reproducibility and data availability | Versioned R scripts, repository, sharing statement | **Conditional**; private until terms permit release |
| TARGET: hypothetical trial, strategy, contrast | Stage 14 trial-specification table (Table S1) | **Partial**; complete time-zero and assumptions descriptions |
| TARGET: observational implementation | Stage 14/15 code | **Present**; distinguish day-zero vs one-year lag |
| TARGET: day-zero cohort flow | Stage 21 `TARGET_day0_cumulative_flow.csv` | **Pending local run** |
| TARGET: day-zero baseline by arm | Stage 21 `TARGET_day0_baseline_by_arm.csv` | **Pending local run** |
| TARGET: day-zero follow-up and missingness | Stage 20 day-zero aggregates | **Present**; reconcile with Stage 21 |
| TARGET: outcomes and effect estimates | Stage 15 and Stage 17 Table 2 | **Present**; interpret as observational estimates |
| TARGET: registration/protocol and ethics | Declarations; RFD non-human-subjects determination | **Documented**; RFD STUDY00030063 (April 7, 2026), no written protocol; page mapping pending |

Do not mark a checklist item "complete" merely because an R script exists. The final journal checklists must map each item to actual page numbers **after** the manuscript is paginated. Keep source definitions and unresolved items distinguishable.


## Confirmed source and regulatory details (October 2026)

RFD STUDY00030063 (April 7, 2026) determined non-human-subjects research; no IRB oversight. No written protocol was prepared. No funding, conflicts, or additional disclosures were reported. The source Ex-US category denotes outside-US healthcare organizations. TriNetX EHR medication records may reflect prescriptions written or administrations depending on context; dispensing and adherence are not verified. Final submission checklist page numbers remain dependent on target-journal layout.
