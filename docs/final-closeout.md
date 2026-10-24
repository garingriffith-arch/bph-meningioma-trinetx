# Numerical reconciliation and reporting utilities

Stages 24–27 provide additional diagnostics and supporting outputs for the previously developed study. They are separate from the core cohort-construction and primary modeling stages.

| Stage | Purpose | Database writes |
| --- | --- | --- |
| 24 | Treatment-effect proportional-hazards diagnostics for selected boundary-corrected models | None; source database accessed read-only |
| 25 | Replay selected historical rare-outcome models and assess final-cohort denominator alignment | None; source database accessed read-only |
| 26 | Cross-file consistency checks and aggregate reporting register | None |
| 27 | Inventory or package locally prepared manuscript and checklist files | None |

## Execution

The `R/run_closeout.R` helper coordinates Stages 24–26 using the environment variables `BPH_R_DIR` (configured script directory), `BPH_WORK_DIR` (authorized local work directory), and `BPH_CLOSEOUT_MODE`.

- `BPH_CLOSEOUT_MODE=all`: run Stages 24–26 when expected outputs are not already present.
- `BPH_CLOSEOUT_MODE=24`, `25` or `26`: run the selected stage.
- `BPH_CLOSEOUT_MODE=audit`: run the aggregate Stage 26 review only.
- `BPH_PH_SCOPE=core` or `all`: choose the subset of Stage 16 fits inspected by Stage 24.

Each stage checks for prerequisite inputs. Stage 25 requires agreement with historical model summaries before accepting revised denominator-aligned fits and does not silently admit zero-day follow-up events. A non-significant proportional-hazards diagnostic is not proof that the proportional-hazards assumption holds.

Stage 27 is a separate document-packaging utility. It does not publish results or modify the source repository. Any output distribution still requires the applicable disclosure review.

## Interpretation

The primary cohort, day-zero sensitivity population and recorded medication-persistence subgroups answer different observational questions. Code-level reconciliation of aggregate estimates does not establish a causal effect or an independent clean-room reproduction from raw data. See [study methods](../reporting/methods-source-verified.md) and the [analysis stage guide](analysis-pipeline.md).
