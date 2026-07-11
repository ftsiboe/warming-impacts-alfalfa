rm(list = ls(all = TRUE));gc()

# Load reusable helper functions
invisible(lapply(list.files("scripts/helpers", pattern = "[.]R$", full.names = TRUE), source))
study_environment <- readRDS("data/study_environment.rds")

create_release_archive(
  source_directory = "output",
  output_directory = "output/releases",
  output_name = "ouputs"
)

# Verify auth first (nice sanity check)
if (requireNamespace("gh", quietly = TRUE)) try(gh::gh_whoami(), silent = TRUE)

# piggyback::pb_release_create(
#   repo = "ftsiboe/WarmingImpactsAlfalfa",
#   tag  = "data",
#   name = "Project Data",
#   body = "This release contains the cleaned and processed datasets used in the project."
# )
#
# piggyback::pb_release_create(
#   repo = "ftsiboe/WarmingImpactsAlfalfa",
#   tag  = "outputs",
#   name = "Project Outputs",
#   body = "This release contains the outputs from project."
# )


# Upload the assets
asset_list <- list.files("data", full.names = TRUE, recursive = TRUE)
piggyback::pb_upload(
  asset_list,
  repo  = "ftsiboe/WarmingImpactsAlfalfa",
  tag   = "data",
  overwrite = TRUE
)

piggyback::pb_upload(
  "output/releases/ouputs.zip",
  repo  = "ftsiboe/WarmingImpactsAlfalfa",
  tag   = "outputs",
  overwrite = TRUE
)





























