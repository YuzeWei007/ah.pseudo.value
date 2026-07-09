# Reproducibility Code for Average Hazard PO and Stacked GMM

This repository reproduces Tables 2 and 3 for the paper on
pseudo-observation (PO) and stacked GMM estimation in average hazard
regression.

The primary workflow has three entry scripts:

1. `main.R` runs all simulation settings.
2. `make_table2.R` creates Table 2 for Model 1.
3. `make_table3.R` creates Table 3 for Model 2.

The statistical helper functions remain in the separate `utility.R`
file used by the original project structure.

## Repository structure

```text
.
├── README.md
├── main.R
├── utility.R
├── make_table2.R
├── make_table3.R
├── analysis/
│   └── run_poplar_analysis.R
└── .gitignore
```

Simulation output CSV files are not distributed in the repository.
The scripts automatically create `results/latest/` and `tables/` when
they run.

## Software requirements

The code requires R and the following packages:

```r
install.packages(c(
  "survival",
  "pseudo",
  "rootSolve",
  "eventglm",
  "openxlsx"
))
```

`openxlsx` is required only for the POPLAR analysis.

## Complete reproduction workflow

First, run the complete 1000-replication simulation:

```bash
Rscript main.R
```

The default settings are:

```text
seed = 2026
nsim = 1000
n = 300
tau = 7
```

The script saves timestamped diagnostic files and the stable file:

```text
results/latest/simulation_summary_long.csv
```

Second, after the simulation finishes, create Tables 2 and 3:

```bash
Rscript make_table2.R
Rscript make_table3.R
```

The table scripts read:

```text
results/latest/simulation_summary_long.csv
```

and create:

```text
tables/table2_model1.tex
tables/table3_model2.tex
```

The generated files are complete LaTeX table environments that can be
included in Overleaf. They require the LaTeX packages `booktabs`,
`graphicx`, and `pdflscape`.

## Quick smoke test

A full 1000-replication run can take substantial time. To check that
the installation and paths work, run a small test first:

```bash
AH_GMM_NSIM=10 Rscript main.R
Rscript make_table2.R
Rscript make_table3.R
```

The smoke-test tables will correctly report `B = 10`. Run the full
simulation afterward to restore the final `B = 1000` results.

## Changing the seed or number of replications

Use environment variables without editing the script:

```bash
AH_GMM_SEED=1234 AH_GMM_NSIM=100 Rscript main.R
```

If these variables are omitted, the script uses seed 2026 and 1000
replications.

## Changing the simulation output path

The default output directory is:

```text
results/latest/
```

To use another directory:

```bash
AH_GMM_OUTPUT_DIR="/absolute/path/to/results" \
  Rscript main.R
```

Paths containing spaces must be placed inside quotation marks.

## Using a different result file for Table 2 or Table 3

Each table script accepts an optional input CSV followed by an optional
output `.tex` path:

```bash
Rscript make_table2.R \
  "/absolute/path/to/simulation_summary_long.csv" \
  "/absolute/path/to/table2.tex"

Rscript make_table3.R \
  "/absolute/path/to/simulation_summary_long.csv" \
  "/absolute/path/to/table3.tex"
```

This is the preferred way to change table paths. No source-code edits
are required.

## Table definitions

Table 2 reports Model 1, with Weibull shape 1. Table 3 reports Model 2,
with Weibull shape 2. Both tables contain four censoring settings and
the following estimators:

- PO
- PO-GS
- PO-GS-QC
- AH-Cox
- Stacked GMM
- GMM-POGS
- GMM-POGS-QC

The reported columns are:

- `RB`: 100 times relative bias, reported as a percentage
- `ESD`: empirical standard deviation
- `MSE`: mean estimated standard error
- `CP`: empirical coverage probability of the 95% Wald interval

## POPLAR real-data analysis

The analysis code is provided in:

```text
analysis/run_poplar_analysis.R
```

The patient-level POPLAR workbook is not included because the data
cannot be distributed through this public repository. Consequently,
the real-data analysis is reproducible only for researchers who have
an authorized local copy of the workbook.

The default expected location is:

```text
data/41591_2018_134_MOESM3_ESM.xlsx
```

The `data/` directory and Excel files are excluded by `.gitignore`.
An authorized user may create that directory locally and run:

```bash
Rscript analysis/run_poplar_analysis.R
```

Alternatively, keep the workbook anywhere outside the repository and
provide its path:

```bash
POPLAR_DATA_FILE="/absolute/path/to/41591_2018_134_MOESM3_ESM.xlsx" \
POPLAR_OUTPUT_DIR="/absolute/path/to/poplar-results" \
  Rscript analysis/run_poplar_analysis.R
```

The script writes its derived outputs locally and does not print
patient-level rows.

## Reproducibility notes

- `main.R` calls `set.seed(2026)` by default.
- Two independently repeated 1000-replication runs were verified to be
  identical before preparing this repository.
- The table scripts validate that every required
  model-censoring-method-parameter combination occurs exactly once.
- Record `sessionInfo()` when archiving a final run so package versions
  are available for future audits.
