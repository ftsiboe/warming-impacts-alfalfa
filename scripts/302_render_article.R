# 302_render_article.R
# Knit the assembled article to Word AND HTML. The master
# narrative/article_warming_effects_alfalfa.Rmd pulls in sections/00..98 and reads
# article_objects.json for inline numbers; figures are pre-rendered into
# output/exhibits/ by 100_WarmingImpactsAlfalfa_exhibits.R. Knit from within NARRATIVE
# so child paths ("sections/.."), the object path ("article_objects.json"), and figure
# includes ("../output/exhibits/..") all resolve. Sourced after 301 by run_article.R.

# Allow standalone runs: load helper paths/constants if not already sourced.
if (!exists("NARRATIVE")) source(file.path("scripts", "300_article_helpers.R"))

master <- file.path(NARRATIVE, "article_warming_effects_alfalfa.Rmd")

rmarkdown::render(
  input         = master,
  output_format = "word_document",
  output_file   = "article_warming_effects_alfalfa.docx",
  knit_root_dir = normalizePath(NARRATIVE),
  quiet         = TRUE)

rmarkdown::render(
  input         = master,
  output_format = "html_document",
  output_file   = "article_warming_effects_alfalfa.html",
  knit_root_dir = normalizePath(NARRATIVE),
  quiet         = TRUE)

message("302_render_article: wrote ",
        file.path(NARRATIVE, "article_warming_effects_alfalfa.{docx,html}"))
