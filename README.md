# 5α-reductase inhibitors and incident meningioma

Reproducibility materials for an active-comparator, new-user electronic health record cohort study of 5α-reductase inhibitors and incident meningioma among men treated for benign prostatic hyperplasia.

## Study and analysis

The primary analysis compared observed initiation of finasteride or dutasteride with tamsulosin, alfuzosin, or silodosin. Risk time began after a one-year outcome-free landmark. A separate day-zero analysis examined treatment-class assignment using information available at the index date. One- and three-year documented-persistence comparisons were exploratory and are not estimates of causal effects of sustained use.

The locked primary analytic population comprised 501,445 patients with 1,456 incident meningiomas. The overlap-weighted primary Cox estimate was HR 1.04 (95% CI 0.89–1.20). Full model outputs are not included because the underlying data and aggregate exports are subject to institutional disclosure requirements.

## Repository structure

- `R/00–04`: source validation, medication extraction, initial cohort and historical models.
- `R/05–06`: revised primary cohort and main model.
- `R/07–12`: persistence, alternative endpoints, sensitivity and exploratory analyses.
- `R/13–17`: historical reporting, target-trial-aligned analyses, censoring weights, final boundary correction, and tables/figures. Stage 13 is superseded by Stage 17.
- `R/18–27`: numerical reconciliation, RECORD-PE and TARGET reports, day-zero cohort description, proportional-hazards checks, rare-endpoint boundary reconciliation and release audit.
- `docs/analysis-pipeline.md`: dependencies, stage inputs/outputs, and execution policy.
- `reporting/`: reporting methodology and checklist mapping.
- `source-provenance/`: checksums and original-source recovery notes.

Stage numbering follows the original research development sequence. The historical files are retained for provenance; **running every script sequentially is not the prescribed final analysis**. The final locked results are obtained from the revised and boundary-corrected stages as specified in the pipeline documentation.

## Reproducibility

Reproduction requires licensed access to the relevant TriNetX Dataworks extract, the corresponding DuckDB database, and the specified R package environment. The archived scripts retain local Windows paths used for the original execution. Review `docs/setup.md` and `docs/analysis-pipeline.md`, configure the project paths, and run against an authorized copy of the data. Never execute historical cohort-construction scripts against the locked database without a backup.

The primary model, FDR family and rare-endpoint outputs were reconciled in independent follow-on checks. The final boundary-corrected treatment-effect proportional-hazards tests and rare-outcome model replay also reconciled. These numerical audits do **not** amount to a full clean-room rerun from raw data.

To check repository hygiene without data access:

```bash
python tests/check_repository.py
```

## Reporting and regulatory status

Reporting follows RECORD-PE for the observational study and TARGET for the day-zero sensitivity analyses. A formal written study protocol was not prepared, and post hoc analyses are not represented as prospectively prespecified.

An institutional Research Feasibility Determination (RFD; **STUDY00030063**, **April 7, 2026**) found that the research did not constitute human subjects research; IRB oversight was not required. No funding, conflicts of interest, or additional disclosures were reported.

The `Ex-US` field indicates records categorized as originating from healthcare organizations outside the United States. TriNetX medication events may reflect orders/prescriptions or administrations according to source context; repeated records cannot establish dispensing or ingestion. These data-source limitations are included in the manuscript and Supplementary Methods.

## Data and code access

This repository contains analysis code and documentation only. It does not contain patient-level records, raw exports, DuckDB files, model artifacts or unreviewed aggregate outputs. Data are not publicly redistributable. **Keep the repository private pending explicit confirmation of institutional and TriNetX code-disclosure requirements.** No open-source license is granted.

## References

- Langan SM, et al. RECORD-PE. *BMJ*. 2018;363:k3532.
- Cashin AG, et al. TARGET Statement. *JAMA*. 2025;334:1084–1093.
- Stein E, et al. TriNetX Dataworks-USA. *Pharmacoepidemiol Drug Saf*. 2025;34:e70198.
