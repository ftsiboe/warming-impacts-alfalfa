# 006/007 gwkit consensus refactor — run & validation notes

## What changed

**`007_WarmingImpactsAlfalfa_boots.R`** was restructured so the SLURM array is
**by bootstrap** (one task = one boot; `--array=1-101` = boot `"0000"` + 100 draws).
Within a task, for each analysis cell (crop × window × climate baseline) the
spec-invariant yield model is computed **once**, then the geographically-weighted
quantities are estimated under **all 50 gwkit distance-preset × kernel
specifications** in parallel across cores and reduced to a **per-county consensus**
(median primary + mean retained) **in-task**. Only the consensus is written — one
file per cell per boot, `consensus_<crop>_period<period>_<NAME>.rds`. No per-spec
files are persisted.

gwkit is now the canonical GW engine (consolidated, class-detected point/polygon API):

| Step | gwkit tool |
|------|-----------|
| Distance grid (10 presets) | `gw_distance_metric_names()` |
| Availability m(zᵢ) — neighbour-weighted, self-excluded, multi-column | `estimate_gwlag()` |
| Associations `cattle~avail00` and `cattle~prod00+prod00_LM` | `estimate_gwr()` |
| County degree-day yield response (GWFE, panel/time) | `estimate_gwr(panel=, time=)` |
| Across-spec consensus (scalar / class) | `gw_consensus_scalar()` / `gw_consensus_class()` |

Removed: the `optimal_gw` cross-validation pick and `output/optimal_gw.rds`
(consensus replaces the single-spec choice).

**Full GWFE county yield response.** 006 now also estimates a per-county
degree-day yield response `b_DD1/b_DD2/b_DD3`, evaluated at each county's own knot
from `output/optimal_knots_gw.rds` (004), via `estimate_gwr(panel="fip",
time="year")` fit once per distinct selected knot pair and reduced to a per-county
median across the 50 specs. The county yield-impact projection uses these slopes,
falling back to the national `PWM.COEF` for aggregate rows / counties without a
county knot. Correspondingly, **`006_..._data_prism_climate.R`** stage 2 now trims
the `prism_climate` dday columns to the **union** of the national (003) and
county-GW (004) knot thresholds so those DD can be built downstream.

**`008_WarmingImpactsAlfalfa_summary.R`** reads the per-cell consensus files, drops
all spec columns, and reports the boot-`"0000"` point estimate plus the across-boot
mean/sd/n keyed by `(crop, period, climate_base)`.

The gwkit estimators live in `packages/gwkit/R/` (`estimate_gwlag.R`,
`estimate_gwr.R`, `estimate_gwss.R`, `gw_consensus.R`, shared engine
`gw_local_engine.R`). They use `GWmodel::gw.dist` / `gw.weight` and a shared local
WLS solve via `stats::lm.wfit`; tests are in `tests/testthat/`.

## Run order

1. **gwkit**: `devtools::document("packages/gwkit"); devtools::test("packages/gwkit")`
   (all green), then `R CMD INSTALL` gwkit into the cluster R library, or rely on
   the `devtools::load_all(...)` call at the top of 006.
2. **boot list** (once): source the `boot list` `function(){...}` block in 006.
3. **bandwidth cache** (once): source the `Bandwidth cache` `function(){...}` block
   in 006 with several cores (`SLURM_CPUS_PER_TASK`) → `output/gw_bandwidths.rds`.
   It caches **two** bandwidths per (cell × spec): `bw_lag` (fixed-distance, for
   `estimate_gwlag`) and `bw_gwr` (adaptive kNN, for `estimate_gwr`), reused across
   all boots. The GWFE yield block selects its own bandwidth once per spec
   (within-spec reuse; not cross-boot cached). Without the cache every boot re-runs
   CV (valid but slower).
4. **array**: `sbatch scripts/job006_boots.sbatch` (101 tasks).
5. **summary**: run `007`.

## Validate BEFORE committing cluster time

- **Single-boot smoke test**: interactively set `Sys.setenv(SLURM_ARRAY_TASK_ID=1)`
  (boot `"0000"`) and `SLURM_CPUS_PER_TASK=4`, source 006, and confirm one
  `consensus_*.rds` is written with finite `availability`/`associations` and a
  finite county `impact_yield`.
- **Consensus vs legacy**: run the preserved
  `006_WarmingImpactsAlfalfa_summary_LEGACY_preConsensus.R` for boot `"0000"`,
  preferred cell, and compare the legacy single-spec `avail00`/associations to the
  new consensus median. They should be close (consensus is more robust, not
  identical).

## Resolved

- **Downstream updated.** `100_..._exhibits.R` and `301_article_objects.R` now read
  the consensus `summary_*` objects (keyed by `crop/period/climate_base`, no spec
  columns); the `optimal_gw[1,]` picks were removed. The `v00_recovery` exhibits
  script and the `101_gw_scalar_consensus_map.R` demo were merged into the single
  `100_..._exhibits.R` (article Figures 4/5 keep their filenames; the recovered
  associations/availability/cattle and predicted-impact figures use
  `gw_consensus_scalar`).

## Caveats / follow-ups

- **Great Circle preset** uses the projected county centroids (as the original code
  did), so its distances are geographically unusual; it is one of 50 specs in the
  ensemble. Revisit if a true lon/lat metric is wanted.
- **Cost**: 50 specs is ~2.5× the old 20-spec preferred-cell work, and the GWFE
  yield block adds (distinct selected knot pairs) × 50 specs local fits per
  boot × cell. The yield model is computed once per cell and specs run in parallel,
  so per-task wall time stays near the old envelope; the increase is CPU-hours. If
  the GWFE block dominates, cache `bw_gwfe` across boots like `bw_lag`/`bw_gwr`.
- **Predicted-impact figures** (`predicted_county_impacts.png`,
  `predicted_impacts_mean.png`) are generated by 100 but not yet referenced by the
  article; the by-kernel / by-distance panels were dropped (no per-spec dimension
  under consensus).
