# Structural Roles of Asymptomatic Infection and a Downstream Bacterial State in a Cholera Model

This repository contains the R scripts used to generate the numerical illustrations and validation checks for the manuscript:

**Structural Roles of Asymptomatic Infection and a Downstream Bacterial State in a Cholera Model**

The scripts generate model-based numerical illustrations only. No individual-level, clinical, or surveillance dataset is analysed.

## Contents

- `simulations_paper1.R`: generates baseline simulations, threshold comparisons, sensitivity results, equilibrium branch calculations, numerical checks, and publication figures.
- `structural_checks_paper1.R`: generates structural model comparisons, identifiability transformation checks, hierarchy comparisons, and endemic eigenvalue calculations.

## Requirements

The scripts require R and the following package:

```r
deSolve
```
Install the required package with:
```r
install.packages("deSolve")
```
## How to run
From the repository directory, run:
```r
Rscript simulations_paper1.R
Rscript structural_checks_paper1.R
```
Alternatively, running:
```r
Rscript structural_checks_paper1.R
```
# Outputs
The scripts create the following folders:
- figures/: publication figures
- tables/: generated parameter values, trajectories, summaries, equilibrium values and structural comparison tables
- validation/: numerical validation checks and session information
## Reproducibility note
All parameter values used in the simulations are illustrative and are reported in the manuscript. The numerical experiments are deterministic model experiments and are not fitted to surveillance data.
## License
This repository is released under the MIT License.
