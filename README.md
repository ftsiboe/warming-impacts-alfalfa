<!-- README.md is generated from README.Rmd. Please edit that file, then knit. -->

# Temperature-Induced Yield Impacts and Geographic Shifts in US Alfalfa Production

Replication project for a study quantifying how warming temperatures affect
county-level alfalfa (and other hay) yields across the contiguous United States,
and how projected warming is likely to shift the geographic distribution of
alfalfa production.

**Authors:** Francis Tsiboe (ftsiboe@hotmail.com; corresponding), Agricultural Risk
Policy Center, North Dakota State University; Hannah Shear and Walker Davis, Department
of Agricultural Economics, Oklahoma State University; Jesse Tack, Department of
Agricultural Economics, Kansas State University; and Jisang Yu, Department of Food and
Resource Economics, Korea University.

## What the study does

Using a county-level unbalanced panel (1951–2022; 2,573 counties) built from USDA
NASS alfalfa/hay production and livestock inventory data merged with PRISM daily
weather, the analysis estimates a piecewise linear degree-day yield-response model
(endogenous temperature thresholds at 14°C and 29°C). It then simulates uniform
warming scenarios (+0.5°C to +3.0°C) to produce county-level yield-shock maps and
assess regional shifts relative to livestock-dense areas.

## Project layout

```
data/          Inputs: processed NASS .rds panels, spatial representation, and
               raw PRISM climate/weather .rds (one file per year, 1950–2024).
scripts/       Numbered analysis pipeline (000–200) + article layer (300–302,
               run_article.R) + HPC .sbatch job files.
  helpers/     Reusable R functions (panel construction, delta method, price
               indices, plotting/theme, NASS processing, environment setup).
narrative/     Manuscript: section .Rmd files under sections/, the master
               article_warming_effects_alfalfa.Rmd, references.bib, reference.docx,
               the generated article_objects.json, and review_comments.Rmd.
output/        Generated results: exhibits/ (figures + figure_data), summary/,
               releases/, and top-level model objects.
fastscratch/   Pipeline scratch (all intermediates; git-ignored, disposable).
```

`data/`, `output/`, and `fastscratch/` are **git-ignored**: research data and
generated outputs are distributed via GitHub Releases (piggyback), not committed.

## Analysis pipeline (`scripts/`)

Sourced from the repository root, in order:

- `000_initialize_WarmingImpactsAlfalfa.R`: sources helpers, creates the output
  tree, sets options/seed.
- `001_WarmingImpactsAlfalfa_data_nass.R`: assemble NASS production & inventory.
- `002_WarmingImpactsAlfalfa_data_prism_weather.R`: county PRISM weather panel.
- `003_WarmingImpactsAlfalfa_availability_associations.R`: geographically weighted
  cattle-forage availability & associations, computed once (boot-invariant).
- `004_WarmingImpactsAlfalfa_knots.R`: national endogenous degree-day thresholds.
- `005_WarmingImpactsAlfalfa_knots_cluster.R`: agro-climatic cluster county knots + slopes.
- `006_WarmingImpactsAlfalfa_data_prism_climate.R`: climate baselines / scenarios.
- `007_WarmingImpactsAlfalfa_boots_cluster.R`: bootstrap estimation of yield
  impacts (no GW; emits impacts only).
- `008_WarmingImpactsAlfalfa_summary.R`: summary objects for the yield-side
  quantities + associations pass-through (`output/summary/`).
- `009_WarmingImpactsAlfalfa_availability_warming.R`: warming-scenario alfalfa
  availability & cattle shifts, applying 003's fixed neighbourhood lag to every
  boot's impacts outside the loop (`output/summary/summary_availability.rds`,
  `summary_cattle.rds`).
- `100_WarmingImpactsAlfalfa_exhibits.R`: figures & tables (`output/exhibits/`).
- `200_WarmingImpactsAlfalfa_study_releases.R`: packaged release archive.

The `job00*.sbatch` files run the heavier steps on an HPC/SLURM cluster.

## Reproducing

To reproduce results from the shipped `data/` and `output/`, run
`003` → `004` → `005` → `007` → `008` → `009` → `100`, then build the article (below).
Scripts `001`, `002`, and `006` rebuild the processed panels / climate scenarios
from **external raw source archives** (NASS Quick Stats, PRISM, gSSURGO) on the
author's machine and are not required when the processed `data/` is present.

## Manuscript (article layer)

Every number quoted inline is computed once from the analysis outputs and read into
the prose from a single JSON file, so nothing is hard-coded:

- `scripts/300_article_helpers.R`: paths, constants, formatting helpers.
- `scripts/301_article_objects.R`: builds `narrative/article_objects.json`.
- `scripts/302_render_article.R`: knits the master to `.docx` + `.html`.
- `scripts/run_article.R`: runs 300 → 301 → 302.

Prose is split into `narrative/sections/*.Rmd` (`00_abstract`, `01_introduction`,
`02_data`, `03_methods`, `04_results`, `05_discussion`, `08_conclusion`,
`97_references`, `98_tables_and_figures`), stitched by
`narrative/article_warming_effects_alfalfa.Rmd`, with citations in
`narrative/references.bib`. Build with:

```r
source("scripts/run_article.R")
```

Reviewer comments and their resolution status are tracked in
`narrative/review_comments.Rmd`.

## Acknowledgments

We gratefully acknowledge USDA NASS, the PRISM Climate Group, Professor Wolfram
Schlenker, NRCS, and the Beocat HPC cluster at Kansas State University. The views
expressed are solely those of the authors and do not reflect the official
positions or policies of their affiliated institutions.
