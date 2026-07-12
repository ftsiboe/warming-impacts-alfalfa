# 005/006 gwkit consensus refactor — run & validation notes

## What changed

**`005_WarmingImpactsAlfalfa_boots.R`** was restructured so the SLURM array is
**by bootstrap** (one task = one boot; `--array=1-101`). Within a task, for each
analysis cell (crop × window × climate baseline) the spec-invariant yield model
is computed **once**, then the geographically-weighted quantities are estimated
under **all 50 gwkit distance-preset × kernel specifications** in parallel across
cores and reduced to a **per-county consensus** (median primary + mean retained)
**in-task**. Only the consensus is written — one file per cell per boot,
`consensus_<crop>_period<period>_<NAME>.rds`. No per-spec files are persisted.

gwkit is now the canonical GW engine:

| Step | gwkit tool |
|------|-----------|
| Distance grid (10 presets) | `gw_distance_metric_names()` |
| Availability m(zᵢ) — neighbour-weighted, self-excluded, multi-column | `estimate_gwlag_by_point()` **(new)** |
| Associations `cattle~avail00` and `cattle~prod00+prod00_LM` | `estimate_gwr_coefficients_by_point()` **(new)** |
| Across-spec consensus diagnostic | `gw_optimal_scalar_by_polygon()` |

Removed: the `optimal_gw` cross-validation pick and `output/optimal_gw.rds`
(consensus replaces the single-spec choice).

**`006_WarmingImpactsAlfalfa_summary.R`** now reads the per-cell consensus files,
drops all spec columns, and reports the boot-`"0000"` point estimate plus the
across-boot mean/sd/n keyed by `(crop, period, climate_base)`.

Two new gwkit functions live in `packages/gwkit/R/estimate_gwlag.R` and
`estimate_gwr_coefficients.R`, both faithful to the original pipeline (they use
`GWmodel::gw.dist` / `gw.weight` and, for the regressions, `stats::lm(weights=)`),
with tests in `tests/testthat/`.

## Run order

1. **gwkit**: `devtools::document("packages/gwkit"); devtools::test("packages/gwkit")`
   (all green), then `R CMD INSTALL` gwkit into the cluster R library, or rely on
   the `devtools::load_all(...)` call at the top of 005.
2. **boot list** (once): source the `boot list` `function(){...}` block in 005.
3. **bandwidth cache** (once): source the `Bandwidth cache` `function(){...}` block
   in 005 with several cores (`SLURM_CPUS_PER_TASK`) → `output/gw_bandwidths.rds`.
   This is the fixed-bandwidth reuse; without it every boot re-runs CV (valid but
   ~2× slower).
4. **array**: `sbatch scripts/job005_boots.sbatch` (101 tasks, 10 cpus, 60G, 3d).
5. **summary**: run `006`.

## Validate BEFORE committing cluster time

- **Single-boot smoke test**: interactively set `Sys.setenv(SLURM_ARRAY_TASK_ID=1)`
  (boot `"0000"`) and `SLURM_CPUS_PER_TASK=4`, source 005, and confirm one
  `consensus_*.rds` is written with finite `availability`/`associations`.
- **Consensus vs legacy**: run the preserved
  `005_WarmingImpactsAlfalfa_boots_LEGACY_preConsensus.R` for boot `"0000"`,
  preferred cell, and compare the legacy single-spec `avail00`/associations to the
  new consensus median. They should be close (consensus is more robust, not
  identical).

## Caveats / follow-ups

- **Downstream not yet updated.** `100_..._exhibits.R`, the v00 recovery script,
  and `301_article_objects.R` still expect the old per-spec schema
  (`optimal_gw[1,]`, `p/theta/kernel/specN`). They must be updated to read the new
  consensus `summary_*` (keyed by `crop/period/climate_base`, no spec columns).
- **Great Circle preset** uses the projected county centroids (as the original
  code did), so its distances are geographically unusual; it is one of 50 specs
  in the ensemble. Revisit if a true lon/lat metric is wanted.
- **Cost**: 50 specs is ~2.5× the old 20-spec preferred-cell work, but the yield
  model is computed once per cell and the specs run in parallel, so per-task wall
  time stays near the old envelope; the increase is CPU-hours.
