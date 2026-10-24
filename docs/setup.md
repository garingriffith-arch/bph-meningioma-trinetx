# Local setup and reproducibility prerequisites

## Software

R 4.x and packages used across analytical stages include `DBI`, `duckdb`, `data.table`, `survival`, `splines`, `ggplot2`, `scales`, `flextable` and `officer`; some optional analyses use `speedglm`. Python 3 is used by the configuration and repository-check utilities.

An exact package lockfile was not recovered. Record the local R session and package versions before attempting a reproducibility run.

## Data and working directories

The R source was retained from the original Windows analysis environment. Several scripts contain absolute paths to the original code directory, restricted TriNetX export, and DuckDB database. These paths are historical configuration values, not publicly available resources.

To prepare a separate copy with workstation-specific paths, run from the repository root:

```bash
python scripts/configure_local_paths.py --raw-dir /path/to/authorized/trinetx-export --work-dir /path/to/separate-analysis-work --out-dir /path/to/configured-r-scripts
```

The configuration tool performs literal path substitution without modifying versioned source files. Review the generated `CONFIGURATION_NOTES.txt`. Confirm local data schemas and dependencies before executing any stage.

## Execution precautions

- Stages 00–17 include table-creation, cohort-construction or modeling operations and may overwrite an existing database or outputs.
- Historical stages are retained for provenance. Final estimates depend on revisions, including the Stage 16 boundary correction and Stage 17 reporting.
- Stages 18–27 support aggregate reporting, diagnostics and reconciliation. Each script documents dependencies.
- Do not rerun historical scripts against an existing validated database without a separate backup and output directory.

See the [stage guide](analysis-pipeline.md). A passing repository hygiene test is not a successful reproduction of the research.

## Repository check

```bash
python tests/check_repository.py
```

The check looks for expected scripts and common restricted-data file types. It does not validate statistical models or inspect every possible disclosure risk.
