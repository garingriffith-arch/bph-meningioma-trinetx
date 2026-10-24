# Analysis code provenance

The R scripts are retained as successive stages of the BPH–meningioma analysis. Their numbering documents the development of cohort construction, endpoint definitions, modeling, sensitivity analyses and reporting.

The original stages contain implementation-specific paths and comments that reflect the environment in which the analysis was developed. Later stages supersede parts of the earlier pipeline; retaining an earlier script does not imply that its results were used in the final analysis.

## Principal revisions

- **Stages 00–04:** Data validation, medication timeline, initial cohort and outcome analysis.
- **Stages 05–06:** Revision of follow-up eligibility and primary modeling.
- **Stages 07–12:** Recorded medication persistence, additional outcomes and exploratory models.
- **Stages 13–17:** Reporting revisions, day-zero sensitivity analyses and the exact-landmark outcome-boundary correction. Stage 17 supersedes Stage 13 for final reporting.
- **Stages 18–27:** Aggregate reconciliation, reporting diagnostics and supporting review utilities.

See [the analysis stage guide](../docs/analysis-pipeline.md) for inputs, outputs and execution precautions.

## Limits of reproducibility

These files constitute a versioned source-code record, not a complete computational archive. The restricted TriNetX extract, original database, exact package-version lockfile and generated outputs are not included. The public repository has not been independently rerun end to end from raw data. No file-hash manifest is distributed with this release.

Independent reproduction requires authorized data access, the applicable data dictionary, a compatible software environment and a separate working database. Running all historical scripts consecutively is not the prescribed procedure for recreating final results.
