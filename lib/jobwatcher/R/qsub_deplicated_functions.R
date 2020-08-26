# Deplicated. This is similar to \code{\link{write_qsubfile}}. (slightly different in "return".
write_job <- function(x, path, recursive, add_time) {
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
  invisible(list(path, time))
}

#' save and \emph{qsub}
#' 
#' @description Deplicated. We recommend to use \code{\link{qsub}}.
#'
#' @param x A character. contents of file.
#' @param script_path A character. The path to write a file.
#' @param script_dir A character. It will concatenated with file_path.
#' @param recursive A logical. Whether make parent directory recursively when it does NOT exist.
#' @param add_time A logical. Whether add the time you execute this function to path for unique naming.
#' @param qsub_args Additional arguments for \emph{qsub}.
#' @seealso \url{https://supcom.hgc.jp/internal/mediawiki/qsub_%E3%82%B3%E3%83%9E%E3%83%B3%E3%83%89}
#' @return Invisible. A list of Job ID, the path you write your file to, and the time you execute this function.
#' @export
save_and_qsub <- function(x, script_path, script_dir, recursive = FALSE,
                          add_time = TRUE, qsub_args = "") {
  time <- character()
  path <- ifelse(is.na(script_dir), script_path, fs::path(script_dir, script_path)) %>% fs::path_abs()
  c(path, time) %<-% write_qsubfile(x, path, recursive, add_time)
  qsubres <- try_system(paste0("qsub ", path, " ", qsub_args))
  rlang::inform(qsubres)
  ID <- stringr::str_split(qsubres, " ")[[1]][3]
  invisible(list(ID, path, time))
}

#' write and \emph{qsub}
#' 
#' @description Deplicated. We recommend to use \code{\link{qsub}}. shorthand of \code{\link{save_and_qsub}}(\code{\link{make_qsubfile}}())
#' @param ... Your codes (default: \emph{bash} codes). Each argument should be a character vector. Multiple arguments and multiple elements will be separated with a line break.
#' @param script_path A character. The path to write a file.
#' @param script_dir A character. It will concatenated with file_path..
#' @param name A character
#' @param first_line A character. It is written in the first line.
#' @param parallel A character
#' @param arrayjob A character
#' @param directory A character
#' @param use_bash_profile A logical. Whether \emph{source ~/.bash_profile} or not.
#' @param other_req A character. Other requirements for \emph{qsub}
#' @param recursive A logical. Whether make parent directory recursively when it does NOT exist.
#' @param add_time A logical. Whether add the time you execute this function to path for unique naming.
#' @param qsub_args Additional arguments for \emph{qsub}.
#' @seealso \url{https://supcom.hgc.jp/internal/mediawiki/qsub_%E3%82%B3%E3%83%9E%E3%83%B3%E3%83%89}
#' @return Invisible. A list of Job ID, the path you write your file to, and the time you execute this function.
#' @export
write_and_qsub <- function(...,
                           script_path, 
                           script_dir = NA_character_,
                           name = NA_character_,
                           first_line = binbash(),
                           parallel = parallel_option(),
                           arrayjob = arrayjob_option(),
                           directory = directory_option(),
                           use_bash_profile = TRUE,
                           other_req = character(0),
                           recursive = FALSE,
                           add_time = TRUE,
                           qsub_args = ""){
  NAME = FIRST_LINE = PARALLEL = ARRAYJOB = DIRECTORY = USE_BASH_PROFILE = OTHER_REQ = SCRIPT_PATH = SCRIPT_DIR = RECURSIVE = ADD_TIME = QSUB_ARGS = NA_character_
  c(NAME, FIRST_LINE, PARALLEL, ARRAYJOB, DIRECTORY, USE_BASH_PROFILE, OTHER_REQ, SCRIPT_PATH, SCRIPT_DIR, RECURSIVE, ADD_TIME, QSUB_ARGS) %<-% 
    list(name, first_line, parallel, arrayjob, directory, use_bash_profile, other_req, script_path, script_dir, recursive, add_time, qsub_args)
  make_qsubfile(..., name = NAME, first_line = FIRST_LINE, parallel = PARALLEL, arrayjob = ARRAYJOB, directory = DIRECTORY, use_bash_profile = USE_BASH_PROFILE, other_req = OTHER_REQ) %>% 
    save_and_qsub(script_path = SCRIPT_PATH, script_dir = SCRIPT_DIR, recursive = RECURSIVE, add_time = ADD_TIME, qsub_args = QSUB_ARGS)
}

#' watch a \emph{qsub}bed job by using \emph{qreport}
#'
#' @param x A list of your job ID, the path of your qsub file, and the time you execute qsub or time before that.
#' @param sys_sleep A numeric. \emph{qreport} interval in seconds.
#' @param max_repeat A integer. Total times of trying \emph{qsub} the same file.
#' @param qsub_args A character. Additional arguments for \emph{qsub/qrecall}.
#' @param modify_req A logical. When re-qsubbing, whether to add recommended requests to qsub_args
#' @param qrecall A logical. Whether use \emph{qrecall -file} instead of \emph{qsub} when re-subbing your job.
#' @param verbose A logical.
#' @param debug A logical.
#' @return Invisible. A list of your final job ID, the path of your qsub file, and the time of final qsub.
#' @export
jobwatch <- function(x, sys_sleep = 60L, max_repeat = 2L, qsub_args = "", modify_req = TRUE, qrecall = FALSE, verbose = FALSE, debug = FALSE){
  #perse ID and time
  if (debug) {verbose <- TRUE}
  x %>%
    purrr::walk(~ assertthat::assert_that(length(.x) == 1)) %>%
    purrr::map(vctrs::vec_cast, character()) -> x
  ID = path = time = NA_character_
  c(ID, path, time) %<-% x
  assertthat::assert_that(!is.na(ID))
  assertthat::assert_that(fs::file_exists(path))
  assertthat::assert_that(stringr::str_length(time) == 12)
  ID_vec <- stringr::str_split(ID, "\\.|-|:")[[1]] %>% as.integer()
  ID_body <- ID_vec[1]
  task <- ID_vec[2:4] %>% seq_int_chr()
  if (verbose) {
    todo(crayon::green(path))
    qsub_verbose(ID_body, task, time)
  }
  counter <- 0
  user = fs::path_home() %>% fs::path_file() %>% as.character()
  while (TRUE) {
    Sys.sleep(sys_sleep)
    rep <- qreport_tbl(ID_body, time, user)
    if (debug) {
      print("qreport: ")
      print(rep)
    }
    rep %>%
      tidyr::replace_na(list(failed_txt = "")) %>%
      dplyr::filter(!!sym("failed_txt") != "Rescheduling") -> rep_filt
    if (nrow(rep_filt) > 0) {
      rep_filt %>% dplyr::mutate_at(dplyr::vars("exit_status", "failed"), as.integer) -> rep_filt
      rep_filt %>% dplyr::filter(!!sym("taskid") %in% task) -> rep_filt
      if (debug) {
        print("filtered: ")
        print(rep_filt)
        print(paste0("setdiff: ", stringr::str_c(dplyr::setdiff(task, rep$taskid), collapse = ", ")))
        print(paste0("sum: ", sum(rep_filt$exit_status, rep_filt$failed)))
      }
      if (identical(dplyr::setdiff(task, rep$taskid), character(0))) {
        if (sum(rep_filt$exit_status, rep_filt$failed) == 0) {
          if (debug) as.data.frame(rep) %>% print()#debug
          message(paste0("'", path, "' has been done.")) #message->stderr, inform->stdout
          if (verbose) rlang::inform(done("'", crayon::cyan(path), "' has been done.")) #message and print
          break
        }else{
          counter <- counter + 1
          if (verbose) {
            fail("The job with",
                 "\n ID: ", crayon::cyan(ID_body),
                 "\n path: ", crayon::cyan(path),
                 "\nhas failed.")
            as.data.frame(rep) %>% print()
          }#debug
          if (counter < max_repeat) {
            qsub_args_new <- qsub_args
            if (modify_req) {
              qsub_args_new <- paste0(qsub_args, " ", rep_filt$recommended_option[1])
              if (qsub_args == "" || length(qsub_args) == 0) {
                readr::read_lines(path) %>% 
                  c(qsub_args_new) %>% 
                  write_job(path, recursive = TRUE, add_time = TRUE) %->% c(path, time)
                qsub_args_new <- qsub_args 
              }
            }
            c(ID, path, time) %<-% qsub(path, qsub_args_new, qrecall)
            if (modify_req) {
              message(paste0("#", counter, " resub: ", path, "\nadditional args: ", qsub_args_new))
            }else{
              message(paste0("#", counter, " resub: ", path))
            }
            ID_vec <- stringr::str_split(ID, "\\.|-|:")[[1]] %>% as.integer()
            ID_body <- ID_vec[1]
            task <- ID_vec[2:4] %>% seq_int_chr()
            if (verbose) {
              if (modify_req) {
                rlang::inform(todo("#", counter, " resub: ", crayon::cyan(path), "\nadditional args: ", qsub_args_new))
              }else{
                rlang::inform(todo("#", counter, " resub: ", crayon::cyan(path)))
              }
              qsub_verbose(ID_body, task, time)
            }
          }else{
            rlang::abort(paste0("'", path, "' has something to cause error or fail."), "qsub_contents_error")
          }
        }
      }
    }
  }
  invisible(list(ID, path, time))
}
