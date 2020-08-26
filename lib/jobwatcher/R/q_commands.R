####send_command####
send_command <- function(command, ID, ..., .hgc_mode = TRUE, .discard_std_error = TRUE) {
  dots <- c(...)
  dots <- dots[!purrr::map_lgl(dots, `==`, "")]
  dots_chr <- stringr::str_c(dots, collapse = " ")
  
  std_error_cmd <- if (.discard_std_error) "2> /dev/null" else character()
  
  if (is.na(ID)) {
    res <- 
      suppressWarnings(system(
        stringr::str_c(command, dots_chr, std_error_cmd, sep = " "),
        intern = TRUE
      ))
  } else {
    if (.hgc_mode) {
      if (!is_number(ID)) res <- ""
      else {
        res <- 
          suppressWarnings(system(
            stringr::str_c(command, "-j", ID, dots_chr, std_error_cmd, sep = " "),
            intern = TRUE
          ))
      }
      
      if (is.null(res) || length(res) == 0L || (length(res) == 1L && res == "")) {
        res <- 
          suppressWarnings(system(
            stringr::str_c(command, "-N", ID, dots_chr, std_error_cmd, sep = " "),
            intern = TRUE
          ))
      }
    } else {# UGE
      res <- 
        suppressWarnings(system(
          stringr::str_c(command, "-j", ID, dots_chr, std_error_cmd, sep = " "),
          intern = TRUE
        ))
    }
  }
  res
}

make_option <- function(x, prefix)  ifelse(is.na(x), "", paste0(prefix, " ", x))

####xml2tbl####
vacant_tbl <- function(colname){
  matrix(nrow = 0, ncol = length(colname)) %>%
    as.data.frame() %>% 
    `colnames<-`(colname) %>% 
    dplyr::as_tibble(.name_repair = "minimal") %>%
    dplyr::mutate_all(as.character)
}

try_xml_to_tbl_internal <- function(xml, command) {
  if (command == "qreport") {
    dplyr::as_tibble(XML::xmlToDataFrame(xml, stringsAsFactors = F), .name_repair = "minimal")
  } else if (command == "qstat") {
    job_list <- XML::xmlToList(xml)$queue_info
    purrr::map_dfr(job_list, ~ dplyr::as_tibble(purrr::compact(.x), .name_repair = "minimal"))
  }
}

try_xml_to_tbl <- function(xml, command = c("qreport", "qstat", "qacct")){
  assign_oneof_default_arg_chr("command")
  tbl_colnames <- 
    switch(command,
           "qreport" = 
             c("JB_owner","JB_job_number",
               "taskid","slots","JB_pe_id","granted_pe",
               "exit_status","failed","queue_name","host_name",
               "JB_name","JAT_qsub_time","JAT_start_time","JAT_end_time",
               "ru_wallclock","cpu","memory","maxvmem","r_mem","r_q","r_cpu",
               "qdel","failed_txt","recommended_queue","recommended_memory","recommended_option"),
           "qstat" = 
             c("JB_job_number", "JAT_prio", "JB_name", "JB_owner", "state", "JAT_start_time", "queue_name", "slots", ".attrs"),
           "qacct" = 
             c("qname", "hostname", "group", "owner", "project", "department", "jobname", "jobnumber", "taskid",
               "account", "priority", "qsub_time", "start_time", "end_time", "granted_pe", "slots",
               "failed", "exit_status", "ru_wallclock", "ru_utime", "ru_stime", "ru_maxrss", "ru_ixrss", "ru_ismrss", "ru_idrss", "ru_isrss", "ru_minflt", "ru_majflt",
               "ru_nswap", "ru_inblock", "ru_outblock", "ru_msgsnd", "ru_msgrcv", "ru_nsignals", "ru_nvcsw", "ru_nivcsw",
               "cpu", "mem", "io", "iow", "maxvmem", "arid", "r_mem", "r_q", "r_cpu"),
           rlang::abort("unknown command was specified in try_xml_to_tbl.")
    )
  tryCatch(
    {
      if ("tbl" %in% class(xml)) tbl <- xml
      else tbl <- try_xml_to_tbl_internal(xml, command)
      missing_colnames <- dplyr::setdiff(tbl_colnames, colnames(tbl))
      tbl <- purrr::reduce(missing_colnames, ~ dplyr::mutate(.x, !!.y := NA_character_), .init = tbl)
      dplyr::select(tbl, !!!tbl_colnames)
    },
    error = function(e) vacant_tbl(tbl_colnames),
    warning = function(e) vacant_tbl(tbl_colnames)
  )
}

####qreport####
qreport_xml_txt <- function(ID, begin = NA, user = NA, end = NA, type) {
  verify_scalar(ID, begin, user, end)
  c(ID, begin, user, end) %<-% .map_as_character(ID, begin, user, end)
  begin_option <- make_option(begin, "-b")
  user_option <- make_option(user, "-o")
  end_option <- make_option(end, "-e")
  xml_option <- ifelse(type == "xml", "-x", "")
  send_command("qreport", ID, user_option, begin_option, end_option, xml_option)
}

#' get \emph{qreport} results in xml fromat
#'
#' @param ID Job ID or Job Name.
#' @param begin A character of \strong{\%Y\%m\%d\%H\%M} format. (optional)
#' @param user Your user ID. (optional)
#' @param end A character of \strong{\%Y\%m\%d\%H\%M} format. (optional)
#' @export
qreport_xml <- function(ID, begin = NA, user = NA, end = NA) qreport_xml_txt(ID, begin, user, end, type = "xml")

#' get \emph{qreport} results as a tibble
#'
#' @inheritParams qreport_xml
#' @export
qreport_tbl <- function(ID, begin = NA, user = NA, end = NA) try_xml_to_tbl(qreport_xml(ID, begin, user, end), "qreport")

#' get \emph{qreport} results
#'
#' @inheritParams qreport_xml
#' @param type result type.
#' @export
qreport <- function(ID, begin = NA, user = NA, end = NA, type = c("tibble", "xml", "text")){
  verify_hgc()
  assign_oneof_default_arg_chr("type")
  if (type == "text") {
    res <- qreport_xml_txt(ID, begin, user, end, type = "text")
  } else {
    res <- qreport_xml_txt(ID, begin, user, end, type = "xml")
    if (type == "tibble") res <- try_xml_to_tbl(res, "qreport")
  }
  res
}

####qacct####
qacct_txt <- function(ID, begin = NA, user = NA, end = NA) {
  verify_scalar(ID, begin, user, end)
  c(ID, begin, user, end) %<-% .map_as_character(ID, begin, user, end)
  begin_option <- make_option(begin, "-b")
  user_option <- make_option(user, "-o")
  end_option <- make_option(end, "-e")
  if (get_jobwatcher_mode() == "hgc") {
    send_command("qreport -c", ID, user_option, begin_option, end_option)
  } else {
    send_command("qacct", ID, user_option, begin_option, end_option, .hgc_mode = FALSE)
  }
}

txt_to_tbl <- function(x) {
  is_spacer_lines <- (x == "" | stringr::str_detect(x, "^(.)\\1*$"))
  is_first_spacer_lines <- is_spacer_lines & !dplyr::lag(is_spacer_lines, default = FALSE)
  info_idx <- cumsum(is_first_spacer_lines)
  info_idx[is_spacer_lines] <- 0L
  purrr::map_dfr(seq_len(sum(is_first_spacer_lines)), ~ {
    vec <- x[info_idx == ..1]
    mat <- stringr::str_split_fixed(stringr::str_remove(vec, "[:space:]*$"), "[:space:]+", 2L)
    df <- data.frame(t(mat[,2L]), stringsAsFactors = FALSE)
    colnames(df) <- mat[,1L]
    dplyr::as_tibble(df, .name_repair = "minimal")
  })
}

# qacct_to_tbl <- function(x) {
#   #TODO
# }


#' get \emph{qacct} results
#' 
#' As for HGC super computer, \emph{qreport --compatible} is called internally. 
#' Unlike \emph{qreport}, \emph{qacct} does NOT have an option to return results in xml format.
#'
#' @inheritParams qreport_xml
#' @param type result type.
#' @export
qacct <- function(ID, begin = NA, user = NA, end = NA, type = c("tibble", "text")) {
  verify_uge()
  assign_oneof_default_arg_chr("type")
  res <- qacct_txt(ID, begin, user, end)
  if (type == "tibble") res <- try_xml_to_tbl(txt_to_tbl(res), "qacct")
  res
}

####qstat####
qstat_xml_txt <- function(ID = NA, user = NA, type) {
  verify_scalar(user)
  c(ID, user) %<-% .map_as_character(ID, user)
  user_option <- make_option(user, "-u")
  xml_option <- ifelse(type == "xml", "-xml", "")
  send_command("qstat", ID, user_option, xml_option, .hgc_mode = FALSE)
}

#' get \emph{qstat} results
#' 
#' @param ID Job ID or Job Name. (optional)
#' @param user Your user ID. (optional)
#' @param type result type.
#' @export
qstat <- function(ID = NA, user = NA, type = c("tibble", "xml", "txt")) { #TODO verify array job or options (e.g. -f)
  verify_uge()
  assign_oneof_default_arg_chr("type")
  if (type == "xml") res <- qstat_xml_txt(ID = NA, user = NA, "xml")
  else {
    res <- qstat_xml_txt(ID = NA, user = NA, "txt")
    if (type == "tibble") res <- try_xml_to_tbl(res, "qstat")
  }
  res
}