# 5α-reductase inhibitors and incident meningioma

Analysis code for an active-comparator, new-user electronic health record study of 5α-reductase inhibitors and incident meningioma among men treated for benign prostatic hyperplasia (BPH).

## Study overview

The primary comparison was initiation of finasteride or dutasteride versus initiation of tamsulosin, alfuzosin, or silodosin. The primary cohort used a one-year outcome-free landmark. Separate day-zero initiation analyses and recorded medication-persistence analyses addressed different design questions.

The final primary analysis included 501,445 individuals and 1,456 incident meningiomas; the reported overlap-weighted Cox hazard ratio was 1.04 (95% CI 0.89–1.20). These observational estimates should not be interpreted as establishing a causal effect. Repeated electronic medication records do not establish dispensing or adherence.

## Repository contents

| Location | Contents |
| --- | --- |
| `R/00`–`R/06` | Export validation, medication history, cohort construction and revised primary analysis |
| `R/07`–`R/12` | Persistence, alternative endpoints and exploratory sensitivity analyses |
| `R/13`–`R/17` | Reporting, day-zero analyses, methodological checks and boundary-corrected models |
| `R/18`–`R/27` | Aggregate validation, reporting diagnostics and reconciliation |
| `R/final_analysis_helpers.R` | Shared propensity-score and survival-model utilities |
| `docs/analysis-pipeline.md` | Stage dependencies, superseded analyses and execution guidance |
| `docs/setup.md` | Local configuration and prerequisites |
| `docs/data-governance.md` | Data access and disclosure constraints |
| `reporting/` | Study definitions, declarations and reporting-guideline mapping |
| `source-provenance/` | Notes on the retained development history |
| `tests/` | Offline repository file-hygiene checks |

The numbered stages reflect the order in which the analysis was developed, **not a single command to rerun every script**. Some early models and reporting outputs were superseded. Stage 13 is superseded by Stage 17. See the [stage guide](docs/analysis-pipeline.md) before executing code.

## Reproducibility

The repository does not include source EHR records, the local DuckDB database, patient-level derivatives, model objects or unpublished aggregate result files. Reproduction requires authorized access to the corresponding TriNetX extract and data schema.

The R scripts preserve historical workstation paths, and some stages replace database tables or analytical outputs. Prepare a separate working copy with `scripts/configure_local_paths.py` and read the [setup instructions](docs/setup.md). Original package versions have not been fully recovered; no exact `renv.lock` is provided. Existing numerical reconciliation checks are not a substitute for an independent end-to-end reproduction.

Run the offline file-hygiene check:

```bash
python tests/check_repository.py
```

## Reproducibility audit and optional source-geography check

The [reproducibility audit](docs/reproducibility-audit.md) documents which parts of the pipeline are source-verifiable and which depend on controlled-access records or missing historical package versions.

On an authorized workstation, run `scripts/record_environment.R` to record the **currently installed** R and package versions in a local report. This is not a reconstruction of the original analysis environment.

An optional, **post hoc** source-geography sensitivity is available in `scripts/check_geography_sensitivity.R`. It reproduces the locked primary model before re-estimating that model after excluding HCOs coded `Ex-US`. It uses DuckDB read-only and does not change the primary study cohort. Any resulting estimate is a sensitivity check, not a replacement primary result.

## Reporting and data access

The observational study is described with reference to RECORD-PE; separate day-zero sensitivity analyses are described with reference to TARGET. The latter should not be mistaken for a randomized trial or a complete emulation of the primary landmark analysis.

An institutional Research Feasibility Determination classified this work as non-human-subjects research, so IRB oversight was not required. No written protocol was prepared before the analyses. See [study declarations](reporting/study-declarations.md).

The code repository is public, but **underlying TriNetX data are not publicly redistributable**. Public availability of this code does not confer access to source data or override applicable institutional and contractual restrictions. No open-source license has been granted.

## Reporting references

- Langan SM, et al. RECORD-PE. *BMJ*. 2018;363:k3532.
- Cashin AG, et al. TARGET Statement. *JAMA*. 2025;334:1084–1093.
- Stein E, et al. TriNetX Dataworks-USA. *Pharmacoepidemiology and Drug Safety*. 2025;34:e70198.
