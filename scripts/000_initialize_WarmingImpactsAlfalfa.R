# 000_initialize_WarmingImpactsAlfalfa.R
# Self-contained project initializer for the reorganized (non-package) layout.
# Run from the repository ROOT. Sourcing this file loads all reusable helper
# functions and ensures the output directory tree exists. The analysis pipeline
# scripts (001-006, 100, 200) re-source the helpers themselves, so this file is a
# convenience for interactive setup and the single place documenting the project's
# directories. It carries NO package machinery (no devtools/roxygen/NAMESPACE).

# Hard reset of workspace
rm(list = ls(all = TRUE)); gc()

# --- Load reusable helper functions ---------------------------------------
# Every function the pipeline needs lives in scripts/helpers/*.R.
invisible(lapply(
  list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))

# --- Ensure the project directory tree ------------------------------------
project_dirs <- c(
  "data",
  "output",
  file.path("output", "summary"),
  file.path("output", "exhibits"),
  file.path("output", "exhibits", "figure_data"),
  file.path("output", "bootstraps"),
  file.path("output", "releases")
)
invisible(lapply(project_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# --- Global options + reproducibility seed --------------------------------
options(scipen = 999L)
options(future.globals.maxSize = 8 * 1024^3)   # ~8 GiB
options(dplyr.summarise.inform = FALSE)
set.seed(1980632)

message("Initialized: helpers sourced; data/ and output/ directories ensured.")
