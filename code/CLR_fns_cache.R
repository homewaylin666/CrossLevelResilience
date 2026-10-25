################################################################################
## Run a module once, then read its saved workspace on later runs
##
## Written by Hung-wei Lin Oct 2026
##
## run_cached(module) replaces source(module) in a runscript. After a module runs,
## the whole workspace is saved to '<cache_dir>/<module>.RData'. On a later run the
## module is skipped while its saved workspace is newer than both the module script and
## the step before it; the first step that needs a run (edited script, or an earlier
## step that ran again) starts from the saved workspace of the step before it, and every
## later step runs again. finish_cached() at the end loads the last saved workspace
## when every step was skipped.
## A saved workspace holds every object of the runscript so far, so the workspace in memory
## is cleared before one is loaded (two copies of the large tables would fill the memory).
## To force a full run: delete the cache folder (or set 'rerun_all <- TRUE' before the steps).
## Source data are outside this check (after a change in data.source, delete the cache).


  # Bookkeeping lives in a hidden environment, so it stays out of the saved workspaces
.cache <- new.env()
.cache$pending <- NULL        # saved workspace waiting to be loaded
.cache$last_time <- -Inf      # time stamp of the step before

  # Load a saved workspace in place of the current one; a file that fails to load is named,
  # with the way out (delete it, and the step that wrote it runs again)
load_cached <- function(file) {
  rm(list = setdiff(ls(globalenv()), c("run_cached", "finish_cached", "load_cached")), envir = globalenv())
  invisible(gc())
  ok <- tryCatch({ load(file, envir = globalenv()); TRUE }, error = function(e) conditionMessage(e))
  if (!isTRUE(ok))
    stop(sprintf("Could not read the saved workspace '%s' (%s). Restart R and run the runscript again; if the error comes back, delete this file (its step runs again).",
                 basename(file), ok), call. = FALSE)
  invisible(TRUE)
}

run_cached <- function(module, cache_dir) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  name  <- tools::file_path_sans_ext(basename(module))
  cache <- file.path(cache_dir, paste0(name, ".RData"))
  rerun <- exists("rerun_all", envir = globalenv()) && isTRUE(get("rerun_all", envir = globalenv()))
  up_to_date <- !rerun && file.exists(cache) &&
    file.mtime(cache) > max(file.mtime(module), .cache$last_time)
  if (up_to_date) {
    cat(sprintf("Up to date: %s (saved %s)\n", name, format(file.mtime(cache), "%d %b %H:%M")))
    .cache$pending <- cache
    .cache$last_time <- file.mtime(cache)
    return(invisible(FALSE))
  }
    # Start from the workspace of the step before, then run this step and save
  if (!is.null(.cache$pending)) load_cached(.cache$pending)
  .cache$pending <- NULL
  cat(sprintf("Run: %s\n", name))
  source(module)
  keep <- setdiff(ls(globalenv()), c("run_cached", "finish_cached", "load_cached"))
  save(list = keep, file = cache, envir = globalenv())
  .cache$last_time <- file.mtime(cache)
  invisible(TRUE)
}

finish_cached <- function() {
  if (!is.null(.cache$pending)) load_cached(.cache$pending)
  .cache$pending <- NULL
  invisible()
}
