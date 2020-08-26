# utils for executing in non-HGC supcom environment

#' get jobwatcher mode
#'
#' @export
#' @return a character. One of "not_uge", "uge", and "hgc".
get_jobwatcher_mode <- function() {
  jobwatcher_mode <- getOption("jobwatcher.mode")
  if (is.null(jobwatcher_mode) || is.na(jobwatcher_mode) || !jobwatcher_mode %in% c("not_uge", "uge", "hgc")) {
    if (all(check_commands_exist(c("qrecall", "qreport")))) jobwatcher_mode <- "hgc"
    else if (all(check_commands_exist(c("qsub", "qdel", "qstat", "qacct")))) jobwatcher_mode <- "uge"
    else jobwatcher_mode <- "not_uge"
    options(jobwatcher.mode = jobwatcher_mode)
  }
  jobwatcher_mode
}

warn_if_first_time_qsub <- function() {
  qsub_history <- getOption("jobwatcher.qsub_history")
  if (is.null(qsub_history) || is.na(qsub_history) || !qsub_history) {
    rlang::warn("This is not HGC supcom environment. 'qsub' command is disabled. (This warning is appeared only once.)")
    options(jobwatcher.qsub_history = TRUE)
  }
  invisible(TRUE)
}

check_commands_exist <- function(vec) {
  vec <- stringr::str_remove(vec, "[:space:].*$")
  if (is_command_exist_which("echo")) purrr::map_lgl(vec, is_command_exist_which)
  else if (is_command_exist_hash("echo")) purrr::map_lgl(vec, is_command_exist_hash)
  else {
    rlang::warn("cannot use neither 'which' nor 'hash' commands. return NAs.")
    rep(NA, length(vec))
  }
}

is_command_exist_which <- function(x) {
  is.null(attr(suppressWarnings(system(paste0("which ", x), intern = TRUE, ignore.stdout = TRUE)), "status"))
}

is_command_exist_hash <- function(x) {
  is.null(attr(suppressWarnings(system(paste0("hash ", x), intern = TRUE, ignore.stdout = TRUE)), "status"))
}
