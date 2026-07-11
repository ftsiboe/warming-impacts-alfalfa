# run_article.R
# Build the alfalfa warming-effects journal article (narrative) end to end.
# Run with the working directory = the repository ROOT.
#
# Prerequisite: the analysis pipeline has been run and its outputs exist under
# output/summary/ and output/exhibits/ (built by scripts 000-006 and
# 100_WarmingImpactsAlfalfa_exhibits.R). This runner does NOT rebuild them; it computes
# the inline numbers (article_objects.json) and knits the document.
rm(list = ls(all = TRUE)); gc()

S <- "scripts"
source(file.path(S, "300_article_helpers.R"))   # paths, constants, formatting helpers
source(file.path(S, "301_article_objects.R"))   # -> narrative/article_objects.json
source(file.path(S, "302_render_article.R"))    # knit the .docx and .html

message("Done. See narrative/article_warming_effects_alfalfa.docx / .html")
