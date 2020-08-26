#' watch a \emph{qsub}bed job by using \emph{qreport}
#'
#' @param ID Your job ID or job name
#' @param path  A path of your qsub file (optional). If unspecified, max_repeat will be set as 0.
#' @param time A character of \strong{\%Y\%m\%d\%H\%M} format. The time you execute qsub or time before that (optional).
#' @param sys_sleep A numeric. \emph{qreport} interval in seconds.
#' @param max_repeat A integer. Total times of trying \emph{qsub} the same file.
#' @param give_up One of "error", "warning", "message". Default is "error".
#' @param qsub_args A character. Additional arguments for re-\emph{qsub}/re-\emph{qrecall}. Arguments written in the original file will be ignored.
#' @param modify_req A logical. When re-qsubbing, whether to add recommended requests to qsub_args
#' @param as_qrecall A logical. Whether use \emph{qrecall -file} instead of \emph{qsub} when re-subbing your job.
#' @param verbose A logical.
#' @param debug A logical.
#' @return Invisible. A list of your final job ID, the path of your qsub file, and the time of final qsub.
#' @export
watch <- function(ID, path = NA, time = NA,
                  sys_sleep = 60L, max_repeat = 2L, 
                  give_up = c("error", "warning", "message"),
                  qsub_args = "", modify_req = FALSE,
                  as_qrecall = FALSE,
                  verbose = FALSE, debug = FALSE){
  verify_hgc()
  #perse ID and time
  if (debug) {verbose <- TRUE}
  verify_scalar(ID, path, time)
  verify_no_na(ID)
  c(ID, path, time) %<-% .map_as_character(ID, path, time)
  if(max_repeat > 0L) verify_file_exists(path)
  if (is.na(path)) path <- "<Your job>"
  assertthat::assert_that(is.na(time) || stringr::str_length(time) == 12)
  assign_oneof_default_arg_chr("give_up")
  
  give_up_fun <- 
    switch(give_up,
           "error" = rlang::abort,
           "warning" = rlang::warn,
           "message" = rlang::inform,
           rlang::abort("function select error", "unexpected_error")
           )
  qsub_qrecall <- ifelse(as_qrecall, qrecall, qsub)
  
  ID_body = task = NULL
  c(ID_body, task) %<-% parse_id(ID)
  if (verbose) {
    rlang::inform(todo(crayon::green(path)))
    rlang::inform(qsub_verbose(ID_body, task, time))
  }
  
  counter <- 0
  user <- get_user_name()
  while (TRUE) {
    Sys.sleep(sys_sleep)
    rep <- qreport(ID_body, time, user, type = "tibble")
    if (debug) {
      print("qreport: ")
      print(rep)
    }
    rep_filt <- 
      rep %>%
      tidyr::replace_na(list(failed_txt = "")) %>%
      dplyr::filter(!!sym("failed_txt") != "Rescheduling")
    if (nrow(rep_filt) > 0) {
      rep_filt <-
        rep_filt %>% 
        dplyr::mutate_at(dplyr::vars("exit_status", "failed"), as.integer) %>% 
        dplyr::filter(!!sym("taskid") %in% task)
      if (debug) {
        print("filtered: ")
        print(rep_filt)
        print(paste0("setdiff: ", stringr::str_c(dplyr::setdiff(task, rep$taskid), collapse = ", ")))
        print(paste0("sum: ", sum(rep_filt$exit_status, rep_filt$failed)))
        }
      if (identical(dplyr::setdiff(task, rep$taskid), character(0))) {
        if (sum(rep_filt$exit_status, rep_filt$failed) == 0) {
          if (debug) print(as.data.frame(rep))#debug
          rlang::inform(done("'", crayon::cyan(path), "' has been done.")) #message->stderr, inform->stdout
          if (verbose) message(paste0("'", path, "' has been done.")) #message and print
          break
        }else{
          counter <- counter + 1
          if (verbose) {
            rlang::inform(fail(
              "The job with",
              "\n ID: ", crayon::cyan(ID_body),
              "\n path: ", crayon::cyan(path),
              "\nhas failed."
            ))
            as.data.frame(rep) %>% print()
            }#debug
          if (counter < max_repeat) {
            qsub_args_new <- qsub_args
            if (modify_req) {
              qsub_args_new <- paste0(qsub_args, " ", rep_filt$recommended_option[1])
              if (qsub_args == "" || length(qsub_args) == 0) {
                path <- write_qsubfile(c(readr::read_lines(path), qsub_args_new), path, recursive = TRUE, add_time = TRUE)
                qsub_args_new <- qsub_args 
               }
            }
            c(ID, path, time) %<-% qsub_qrecall(path, qsub_args_new, watch = FALSE)
            if (modify_req) {
              rlang::inform(paste0("#", counter, " resub: ", path, "\nadditional args: ", qsub_args_new))
            }else{
              rlang::inform(paste0("#", counter, " resub: ", path))
            }
            c(ID_body, task) %<-% parse_id(ID)
            if (verbose) {
              if (modify_req) {
                rlang::inform(todo("#", counter, " resub: ", crayon::cyan(path), "\nadditional args: ", qsub_args_new))
              }else{
                rlang::inform(todo("#", counter, " resub: ", crayon::cyan(path)))
              }
              rlang::inform(qsub_verbose(ID_body, task, time))
              }
          }else{
            give_up_fun(paste0("'", path, "' has something to cause error or fail."), "qsub_contents_error")
          }
        }
      }
    }
  }
  invisible(list(ID = ID, path = path, time = time))
}

#' make function for qsub a job and watch progress
#' 
#' @description Short hand of creating a function 
#'   doing \code{\link{qsub}} with \code{watch = TRUE}. 
#'   The created function has a dammy argument which has no effect.
#'
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
#' @param additional_args A list. Elements are passed to \code{\link{qsub}}
#' @seealso \url{https://supcom.hgc.jp/internal/mediawiki/qsub_%E3%82%B3%E3%83%9E%E3%83%B3%E3%83%89}
#' @return A function which has a dammy argument
#' @export
qsub_function <- function(...,
                          script_path,
                          script_dir = NA_character_,
                          name = NA_character_,
                          first_line = binbash(),
                          parallel = parallel_option(),
                          arrayjob = arrayjob_option(),
                          directory = directory_option(),
                          use_bash_profile = FALSE,
                          other_req = character(0),
                          recursive = FALSE,
                          add_time = TRUE,
                          qsub_args = "", 
                          additional_args = list()){
  force(list(name, first_line, parallel, arrayjob, directory, use_bash_profile, other_req, script_path, script_dir, recursive, add_time, qsub_args))
  function(dammy_arg){
    qsubfile <- make_qsubfile(..., name = name, first_line = first_line, parallel = parallel, arrayjob = arrayjob, directory = directory, use_bash_profile = use_bash_profile, other_req = other_req)
    path <- write_qsubfile(x = qsubfile, path = fs::path(script_dir, script_path), recursive = recursive, add_time = add_time)
    do.call(qsub, c(list(path = path, args = qsub_args), additional_args))
  }
}
