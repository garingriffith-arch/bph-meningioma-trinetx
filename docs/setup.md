# Local setup

## Locations

Historical R scripts expect:

- R scripts: `C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/bph_m/R`
- Source export: `C:/Users/garin/OneDrive/OHSU/Research/NRSG/TriNetX/Data/BPH`
- Working area: `C:/TriNetX/bph_m_work`
- DuckDB: `C:/TriNetX/bph_m_work/duckdb/bph_m.duckdb`

These paths occur in the **original historical scripts**. Do not globally search-and-replace them in the preserved archive without tracking the change and validating downstream inputs. Stage 21 supports `BPH_WORK_DIR` as an environment override.

## Packages

R (the original local environment used R 4.x), with `DBI`, `duckdb`, `data.table`, `survival`, `splines` (base recommended), `ggplot2`, `scales`, `flextable`, `officer` and optional `speedglm`. Install only from your standard approved R package mirror. This repository does not invent an `renv.lock` file: package versions from the exact analysis execution environment have not been recovered. Capture `sessionInfo()` before any reproducibility rerun.

## Reproduction prerequisites

Stage 00 requires the licensed, locally exported TriNetX files. Stages 01–16 materialize tables in the same DuckDB database. Stages 15–21 assume prerequisite cohort tables and final analytical outputs already exist.

For current manuscript reporting, use the validated existing DuckDB. Run only read-only stages until authorized to reconstruct the entire cohort. Back up the `.duckdb` file and outputs if any model or source stage must be rerun.

## Data boundaries

Do not place study extracts or the working DuckDB inside the repository. Do not publish generated CSVs before reviewing institutional minimum-cell disclosure requirements. If a generated file contains IDs or dates of care, treat it as restricted even when deidentified at source.

## Configure paths without modifying the source archive

The historical scripts record the Windows locations used for their original execution. For a different authorized workstation, create a separate configured copy of the R scripts:

```bash
python scripts/configure_local_paths.py \
  --raw-dir /path/to/authorized/trinetx-csvs \
  --work-dir /path/to/analysis-work \
  --out-dir /path/to/local-r-scripts
```

This is a literal path substitution, not a reconstruction or scientific change. Review `CONFIGURATION_NOTES.txt` before running. The raw data schema, appropriate package versions and the full input history remain necessary for independent reproduction; the original `R/` code is not altered.
