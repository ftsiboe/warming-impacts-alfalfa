#' Create a ZIP archive from a directory
#'
#' Creates a compressed ZIP archive by copying selected contents of
#' `source_directory` into a temporary staging directory and then zipping that
#' staged structure into `output_directory` using `output_name` as the base
#' filename.
#'
#' Files are selected using optional pattern matching. Only files whose full
#' paths match at least one of the supplied patterns are included in the archive.
#' Directory structure relative to `source_directory` is preserved in the ZIP.
#'
#' @param source_directory Character scalar. Path to the folder whose contents
#'   will be archived. Must exist.
#' @param output_directory Character scalar. Path to the folder where the ZIP
#'   archive will be written. Created if it does not exist.
#' @param output_name Character scalar. Base name for the archive (without
#'   `.zip`), e.g., `"my_archive"`.
#' @param patterns Optional character vector of regular expression patterns.
#'   Only files whose full paths match at least one pattern are included.
#'   If `NULL`, all files under `source_directory` are included.
#'
#' @details
#' The function stages files in a temporary directory before zipping in order to
#' control the internal directory structure of the archive. The working
#' directory is temporarily changed during ZIP creation and restored on exit.
#'
#' ZIP compression is performed using `utils::zip()` with recursive, maximum
#' compression flags (`-r9X`).
#'
#' @return Character scalar. Full file path to the created ZIP archive.
#'
#' @export
create_release_archive <- function(
    source_directory,
    output_directory,
    output_name,
    patterns = NULL
) {
  # basic checks
  stopifnot(is.character(source_directory), length(source_directory) == 1L, nzchar(source_directory))
  stopifnot(is.character(output_directory), length(output_directory) == 1L, nzchar(output_directory))
  stopifnot(is.character(output_name),      length(output_name)      == 1L, nzchar(output_name))
  
  folder <- normalizePath(source_directory, winslash = "/", mustWork = TRUE)
  
  out_dir <- normalizePath(output_directory, winslash = "/", mustWork = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # stage into temp to control zip structure
  tmp <- file.path(tempdir(), output_name)
  if (dir.exists(tmp)) unlink(tmp, recursive = TRUE, force = TRUE)
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  
  list_of_files <- list.files(folder, full.names = TRUE, all.files = TRUE, no.. = TRUE, recursive = TRUE)
  list_of_files <- list_of_files[grepl(paste0(patterns,collapse = "|"),list_of_files)]
  
  file.copy(
    from      = list_of_files,
    to        = tmp,
    recursive = TRUE
  )
  
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp)
  
  zip_path <- file.path(out_dir, paste0(output_name, ".zip"))
  if (file.exists(zip_path)) file.remove(zip_path)
  
  utils::zip(
    zipfile = zip_path,
    files   = ".",
    flags   = "-r9X"
  )
  
  zip_path
}
