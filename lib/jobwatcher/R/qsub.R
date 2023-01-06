#' make a file suitable for \emph{qsub}
#'
#' @param ... Your codes (default: \emph{bash} codes). Each argument should be a character vector. Multiple arguments and multiple elements will be separated with a line break.
#' @param name A character
#' @param first_line A character. It is written in the first line.
#' @param parallel A character
#' @param arrayjob A character
#' @param directory A character
#' @param use_bash_profile A logical. Whether \emph{source ~/.bash_profile} or not.
#' @param other_req A character. Other requirements for \emph{qsub}
#' @seealso \url{https://supcom.hgc.jp/internal/mediawiki/qsub_%E3%82%B3%E3%83%9E%E3%83%B3%E3%83%89}
#' @return qsub script as a character. In order to write this in a file, use \code{write} or \code{\link{write_qsubfile}}.
#' @export
make_qsubfile <- function(...,
                          name = NA_character_,
                          first_line = binbash(),
                          parallel = parallel_option(),
                          arrayjob = arrayjob_option(),
                          directory = directory_option(),
                          use_bash_profile = FALSE,
                          other_req = ""){#TODO docker file home directory check
  dots_parser(...) -> inputs
  list(name, first_line, parallel, arrayjob, directory, other_req) %>%
    purrr::walk(~ assertthat::assert_that(rlang::is_character(.x)))
  stringr::str_c(
    first_line,
    parallel,
    arrayjob,
    ifelse(is.na(name), "", resource("-N", name)) %>% character_1_0(),
    directory,
    ifelse(!is.na(use_bash_profile) && use_bash_profile, grov_env(), "") %>% character_1_0(),
    other_req,
    "##########",
    inputs,
    sep = "\n")
}

#' write a qsub file
#'
#' @param x Your qsub script.
#' @param path A character. The path to write a file.
#' @param recursive A logical. Whether make parent directory recursively when it does NOT exist.
#' @param add_time A logical. Whether add the time you execute this function to path for unique naming.
#' @return invisible. The path where you actually write your file.
#' @export
write_qsubfile <- function(x, path, recursive, add_time) {
  assertthat::assert_that(is.character(x))
  verify_path(path, recursive)
  time <- format(Sys.time(), "%Y%m%d%H%M")
  if (add_time) {
    ext <- fs::path_ext(path)
    if (ext == "") {
      path <- stringr::str_c(fs::path_ext_remove(path), "_", time)
    }else{
      path <- stringr::str_c(fs::path_ext_remove(path), "_", time, ".", fs::path_ext(path))
    }
  }
  write(x, path, append = F)
  invisible(path)
}

#' write out and \emph{qrecall}
#'
#' @param ... Paths to recall. Each argument should be a character vector. Multiple arguments and multiple elements will be separated with a line break.
#' @param path A character. The path to write a file for \emph{qrecall}. If path is a directory, prefix 'qrecall' will be added to the file name.
#' @param log_path A character (optional). The path to write a stdout of \emph{qrecall}. Default: home directory.
#' @param recursive A logical. Whether make parent directory recursively when it does NOT exist.
#' @param add_time A logical. Whether add the time you execute this function to path for unique naming.
#' @param watch A logical.  Whether trace your \emph{qsub}bed/\emph{qrecall}ed job by using \code{\link{watch}}.
#' @return Invisible. A list of Job ID, the path where you write your file, and the time you execute this function.
#' @export
write_and_qrecall <- function(..., path = fs::path_home(), log_path = NA_character_, recursive = FALSE, add_time = TRUE, watch = FALSE) {
  inputs <- dots_parser(...)
  time <- character()
  if (fs::is_dir(path)) path <- fs::path(path, "qrecall")
  path <- write_qsubfile(inputs, path, recursive, add_time)
  if (!is.null(log_path) && !is.na(log_path)) verify_path(log_path, recursive)
  qrecall(path, args = make_option(log_path, " -o"), watch = watch)
}


#' \emph{qsub} a file
#'
#' @param path A character. The path to a \emph{qsub/qrecall} file.
#' @param args A character. Additional arguments for \emph{qsub/qrecall}. Arguments written in the original file will be ignored.
#' @param watch A logical.  Whether trace your \emph{qsub}bed/\emph{qrecall}ed job by using \code{\link{watch}}.
#' @param verbose A logical.
#' @param ... Additional arguments passed directly to \code{\link{watch}} if watch is TRUE.
#' @return Invisible. A list of Job ID, the path you write your file to, the time you execute this function, and if you specified '-sync y' flag, the exit status.
#' @export
qsub <- function(path, args = NA, watch = FALSE, verbose = TRUE, ...){
  verify_scalar(path, args)
  path <- fs::path_abs(path)
  verify_file_exists(path)
  args <- .as_character(args)
  args <- ifelse(is.na(args) || args == "", "", paste0(" ", args))
  time <- format(Sys.time(), "%Y%m%d%H%M")
  if (get_jobwatcher_mode() %in% c("hgc", "uge")) {
    qsubres <- system(paste0("qsub ", path, args), intern = TRUE)
    rlang::inform(qsubres)
    exit_code <- stringr::str_subset(qsubres, "exited with exit code")
    id <- stringr::str_subset(qsubres, "Your .+ been submitted")
    id <- stringr::str_extract(id, "Your .+ been submitted")
    ID <- stringr::str_split(id, " ")[[1]][3]
    res <- list(ID = ID, path = path, time = time)
    if (length(exit_code) > 0) {
      exit_code <- stringr::str_extract(exit_code, "exit code [:digit:]+")
      exit_code <- stringr::str_remove(exit_code, "exit code ")
      exit_code <- as.integer(exit_code)
      res <- list(ID = ID, path = path, time = time, exit_code = exit_code)
    }
    if (watch) res <- watch(ID, path, time, ..., as_qrecall = FALSE)
  } else {# not-uge
    if (verbose && !is.na(args)) rlang::inform("This is not UGE environment. argument 'args' will be ignored.")
    if (fs::file_access(path, 'execute') && !is.na(read_shebang())) command <- path
    else {
      if (verbose) rlang::inform("Your file is not executable or shebang is not detected in your file. Your file will be executed on '/bin/bash'.")
      command <- paste0("/bin/bash ", path)
      }
    exit_code <- system(command, intern = FALSE, wait = TRUE)
    res <- list(ID = NA_character_, path = path, time = time, exit_code = exit_code)
  }
  invisible(res)
}

#' \emph{qrecall} a file
#'
#' @inheritParams qsub
#' @return Invisible. A list of Job ID, the path you write your file to, and the time you execute this function.
#' @export
qrecall <- function(path, args = NA, watch = FALSE, ...){
  verify_scalar(path, args)
  path <- fs::path_abs(path)
  verify_file_exists(path)
  args <- .as_character(args)
  args <- ifelse(is.na(args) || args == "", "", paste0(" ", args))
  time <- format(Sys.time(), "%Y%m%d%H%M")

  if (get_jobwatcher_mode() != "hgc") {
    rlang::inform("This is not HGC environment. 'qrecall' command will be ignored.")
    res <- list(ID = NA_character_, path = path, time = time)
  } else {
    qsubres <- system(paste0("qrecall -file ", path, args), intern = TRUE)
    rlang::inform(qsubres)
    ID <- stringr::str_split(qsubres, " ")[[1]][3]
    res <- list(ID = ID, path = path, time = time)
    if (watch) res <- watch(ID, path, time, ..., as_qrecall = TRUE)
  }
  invisible(res)
}
