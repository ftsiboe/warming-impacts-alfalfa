# Temperature-Induced Yield Impacts and Geographic Shifts in US Alfalfa Production

Replication project for a study quantifying how warming temperatures affect
county-level alfalfa (and other hay) yields across the contiguous United States,
and how projected warming is likely to shift the geographic distribution of
alfalfa production.

**Author:** Francis Tsiboe (ftsiboe@hotmail.com), Agricultural Risk Policy Center,
North Dakota State University.

## What the study does

Using a county-level unbalanced panel (1951–2022; 2,573 counties) built from USDA
NASS alfalfa/hay production and livestock inventory data merged with PRISM daily
weather, the analysis estimates a piecewise linear degree-day yield-response model
(endogenous temperature thresholds at 14°C and 29°C). It then simulates uniform
warming scenarios (+0.5°C to +3.0°C) to produce county-level yield-shock maps and
assess regional shifts relative to livestock-dense areas.

## Project layout

```
data/         Inputs: processed NASS .rds panels, spatial representation, and
              raw PRISM climate/weather .rds (one file per year, 1950–2024).
scripts/      Numbered analysis pipeline (000–200) plus HPC .sbatch job files.
  helpers/    Reusable R functions (panel construction, delta method, price
              indices, plotting/theme, NASS processing, environment setup).
narrative/    Manuscript: section .Rmd files under sections/, source .docx drafts,
              and review_comments.Rmd (tracked reviewer comments + proposed fixes).
output/       Generated results: exhibits/ (figures + figure_data), summary/,
              releases/, and top-level model objects (.rds / .xlsx).
```

## Analysis pipeline (`scripts/`)

Run in order:

- `000_initialize_WarmingImpactsAlfalfa.R` — environment / paths.
- `001_WarmingImpactsAlfalfa_data_nass.R` — assemble NASS production & inventory.
- `002_WarmingImpactsAlfalfa_data_prism_weather.R` — county PRISM weather panel.
- `003_WarmingImpactsAlfalfa_knots.R` — endogenous degree-day thresholds.
- `004_WarmingImpactsAlfalfa_data_prism_climate.R` — climate baselines / scenarios.
- `005_WarmingImpactsAlfalfa_boots.R` — bootstrap estimation.
- `006_WarmingImpactsAlfalfa_summary.R` — summary objects (`output/summary/`).
- `100_WarmingImpactsAlfalfa_exhibits.R` — figures & tables (`output/exhibits/`).
- `200_WarmingImpactsAlfalfa_study_releases.R` — packaged release archive.

The `job00*.sbatch` files run the heavier steps on an HPC/SLURM cluster.

## Manuscript

The working paper lives in `narrative/`. Prose is split into
`narrative/sections/*.Rmd` (`00_abstract`, `01_introduction`, `02_data`,
`03_methods`, `04_results`, `05_discussion`, `08_conclusion`, `97_references`,
`98_tables_and_figures`). Open reviewer comments and their proposed fixes are
tracked in `narrative/review_comments.Rmd`.

## Acknowledgments

We gratefully acknowledge USDA NASS, the PRISM Climate Group, Professor Wolfram
Schlenker, NRCS, and the Beocat HPC cluster at Kansas State University. The views
expressed are solely those of the authors and do not reflect the official
positions or policies of their affiliated institutions.
