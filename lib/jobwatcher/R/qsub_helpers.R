###define functions for qsub###

#' shortcut for memory requirement
#' @param x A numeric scalar.
#' @export
#' @examples
#' mem(5)
#' resource("-l", mem(5.3))
mem <- function(x) {
  verify_scalar(x)
  x <- vctrs::vec_cast(x, numeric())
  paste0("s_vmem=", x, "G,mem_req=", x, "G")
}


#' resource requirement
#' @param ... Requirements. Each argument should be a character vector. Multiple arguments and multiple elements will be separated with a space.
#' @export
#' @examples
#' resource("-l", c("def_slot", "1"))
resource <- function(...) {
  inputs <- dots_parser(..., sep_collapse = " ")
  paste0("#$ ", inputs)
}

character_1_0 <- function(x) if (x == "") character(0) else x

#' Arrayjob requirement
#' @description Generates "-t" argument for "qsub".
#' @param n An integer which represents number of job-array. If n is 1, arrayjob is not required.
#' @param tc An integer. Max mumber of tasks executed simultaneously.
#' @param stepsize An integer(option). SGE_TASK_STEPSIZE option of "qsub -t".
#' @examples
#' arrayjob_option(10L)
#' @export
arrayjob_option <- function(n = 1L, tc = 100L, stepsize = NULL) {
  verify_scalar(n, tc)
  n <- vctrs::vec_cast(n, integer()) %>% as.character()
  tc <- vctrs::vec_cast(tc, integer()) %>% as.character()
  one2n <- paste0("1-", n)
  if (!is.null(stepsize)) {
    verify_scalar(stepsize)
    stepsize <- vctrs::vec_cast(stepsize, integer())
    one2n <- paste0(one2n, ":", stepsize)
  }
  resource("-t", one2n, "-tc", tc)
}

## @usage parallel_option(env = c("def_slot", "mpi", "mpi-fillup", "mpi_4", "mpi_8", "mpi_16", "mpi_24"), slot = 1L, memory = 5.3, master_memory = NULL, ljob = FALSE, no_rerun = TRUE, special_que = c(NULL, "cp", "docker", "knl", "gpu", "groupname"), docker_images = NA_character_)

#' Parallel job requirement
#' @description Generates parallel-job related arguments for "qsub".
#' @param env A character. Choose one of "def_slot", "mpi", "mpi-fillup", "mpi_4", "mpi_8", "mpi_16", "mpi_24".
#' @param slot An integer. The number of slots. If env is either "mpi" or "mpi-fillup", 2-length integer vector representing minimun and maximun nubmer of slots is also accepted.
#' @param memory A double. Memory requirement(Gb).
#' @param master_memory A double (option). Memory requirement for the master que. If \emph{slot} is 1 or \emph{master_memory} is equal to \emph{memory}, this argument will be ignored.
#' @param ljob A logical. Whether need to run more than 2 days. if you require more than 128Gb in total, this option is automatically set as FALSE.
#' @param no_rerun A logical. Whether allow to run on rerun ques.
#' @param special_que A character (option). Choose one of "cp", "docker", "knl", "gpu", "\emph{groupname}", "exclusive"(equivalent to groupname). If specified, ljob option will be ignored.
#' @param docker_images A character (option). Valid only if special_que == "docker".
#' @examples
#' parallel_option(slot = 4L, memory = 10, master_memory = 5, ljob = TRUE)
#' @export
parallel_option <- function(env = c("def_slot", "mpi", "mpi-fillup", "mpi_4", "mpi_8", "mpi_16", "mpi_24"),
                            slot = 1L, memory = 5.3, master_memory = NULL,
                            ljob = FALSE, no_rerun = TRUE, special_que = NULL, docker_images = NA_character_){
  resource_df <- 
    dplyr::tibble(
      ENV = c("def_slot", "mpi", "mpi-fillup", "mpi_4", "mpi_8", "mpi_16", "mpi_24"),
      BASE_SLOT = c(1L, 1L, 1L, 4L, 8L, 16L, 24L)
    )
  if (!is.null(special_que)) {
    verify_scalar(special_que)
    verify_in(special_que, c("cp", "docker", "knl", "gpu", "groupname", "exclusive"))
    if (special_que == "docker") {
      verify_scalar(docker_images)
      verify_no_na(docker_images)
      docker_images <- .as_character(docker_images)
    }
  }
  
  assign_oneof_default_arg_chr("env")
  if (env %in% c("mpi", "mpi-fillup")) {
    verify_in(length(slot), c(1, 2))
    if (length(slot) == 1) slot <- rep(slot, 2)
    max_slot <- ifelse(env == "mpi", 1L, vctrs::vec_cast(slot[2], integer()))
  }else{
    verify_scalar(slot)
    max_slot <- slot
  }
  slot <- vctrs::vec_cast(slot, integer())
  resource_df <- 
    resource_df %>%
    dplyr::mutate(SLOT_NODE = c(max_slot, max_slot, max_slot, 4L, 8L, 16L, 24L))
  verify_scalar(memory)
  memory <- vctrs::vec_cast(memory, numeric())
  if (is.null(master_memory)) {
    master_memory <- memory
  }else{
    verify_scalar(master_memory)
    master_memory <- vctrs::vec_cast(master_memory, numeric())
  }
  ljob <- vctrs::vec_cast(ljob, logical())
  no_rerun <- vctrs::vec_cast(no_rerun, logical())

  df <- 
    resource_df %>%
    dplyr::filter(!!dplyr::sym("ENV") == env)
  slot_node <- df$SLOT_NODE
  base_slot <- df$BASE_SLOT
  assertthat::assert_that(slot[1] %% base_slot == 0)
  total_memory <- slot_node * memory + ifelse(slot_node > 1L, master_memory - memory, 0)
  if (stringr::str_detect(env, "^mpi")) total_memory <- max(total_memory, slot_node * memory)

  lmem <- FALSE
  special_resource <- ""

  if (!is.null(special_que)) {
    if (ljob)  {
      ljob <- FALSE
      rlang::warn("Special que request. ljob option was set as FALSE.", "requirement_resource_warning")
    }
    if (special_que == "cp") {
      max_slot <- 12L
      max_memory <- 128L
    }else if (special_que == "docker") {
      max_slot <- 12L
      max_memory <- 128L
      special_que <- paste0('docker,docker_images="', docker_images, '"')
    }else if (special_que == "knl") {
      max_slot <- 64L
      max_memory <- 112L
      if (env != "def_slot") rlang::abort("env must be def_slot when knl que is specified.", "requirement_resource_error")
    }else if (special_que == "gpu") {
      max_slot <- Inf #FIXME
      max_memory <- 1000L
      if (env != "def_slot") rlang::abort("env must be def_slot when gpu que is specified.", "requirement_resource_error")
    }else if (special_que %in% c("groupname", "exclusive")) {
      max_slot <- 24L
      max_memory <- 128L
      special_que <- "exclusive"
      # system("id", intern = TRUE) -> id
      # assertthat::assert_that(id$status == 0)
      # stringr::str_split(id$stdout, ("\\)|\\("))[[1]][4] %>% paste0(".q") -> groupque
      # special_resource <- resource("-q", groupque)
    }
    if (slot_node > max_slot) rlang::abort("Large number of slot request.", "requirement_resource_error")
    if (total_memory > max_memory) rlang::abort("Large memory request.", "requirement_resource_error")
    special_resource <- resource("-l", special_que)
    result <- 
      stringr::str_c(
        resource("-pe", env, slot),
        special_resource,
        resource("-l", mem(memory)),
        ifelse(slot > 1L && master_memory != memory,
               resource("-masterl", mem(master_memory)),
               "") %>% character_1_0(),
        sep = "\n"
      )
  }else{
    if (slot_node > 39L) {
      rlang::abort("number of slots must be equal to or less than 24 per node.", "requirement_resource_error")
    }
    if (slot_node > 24L) {
      if (ljob) {
        ljob <- FALSE
        rlang::warn("Large number of slots request. lmem option was selected instead of ljob option. Running time is allowed up to 2 weeks.", "requirement_resource_warning")
      }else{
        rlang::inform("Large number of slots request. lmem option was selected.", "requirement_resource_message")
      }
      lmem <- TRUE
    }
    if (total_memory > 2000L) {
      rlang::abort("total memory must be equal to or less than 2Tb per node.", "requirement_resource_error")
    }
    if (total_memory > 128L) {
      if (ljob) {
        ljob <- FALSE
        rlang::warn("Large memory request. lmem option was selected instead of ljob option. Running time is allowed up to 2 weeks.", "requirement_resource_warning")
      }else{
        rlang::inform("Large memory request. lmem option was selected.", "requirement_resource_message")
      }
      lmem <- TRUE
    }

    result <- 
      stringr::str_c(
        resource("-pe", env, stringr::str_c(slot, collapse = "-")),
        resource("-l", mem(memory)),
        ifelse(!is.na(ljob) && ljob, resource("-l", "ljob"), "") %>% character_1_0(),
        ifelse(!is.na(no_rerun) && no_rerun, resource("-q", "'!mjobs_rerun.q'"), "") %>% character_1_0(),
        ifelse(!is.na(lmem) && lmem, resource("-l", "lmem"), "") %>% character_1_0(),
        ifelse(max_slot > 1L && master_memory != memory,
               resource("-masterl", mem(master_memory)),
               "") %>% character_1_0(),
        ifelse(max_slot > 1L && master_memory != memory && ljob,
               resource("-masterl", "ljob"),
               "") %>% character_1_0(),
        ifelse(max_slot > 1L && master_memory != memory && no_rerun,
               resource("-masterq", "'!mjobs_rerun.q'"),
               "") %>% character_1_0(),
        ifelse(max_slot > 1L && master_memory != memory && lmem,
               resource("-masterl", "lmem"),
               "") %>% character_1_0(),
        sep = "\n"
      )
  }
  result
}

#' shebang and -S specification for /bin/bash
#' @export
binbash <- function() "#!/bin/bash\n#$ -S /bin/bash"

#' use ~/.bash_profile
#' @export
grov_env <- function() paste0("source ", fs::path_expand("~/.bash_profile"))

convert_to_array <- function(x) {
  x[is.na(x)] <- "" #escape NA in order not to return NA
  stringr::str_c("[", seq_along(x), ']="', x, '"', collapse = " ")
}

is_bash_name <- function(x) {
  is.character(x) &
    stringr::str_detect(x, "^([:alnum:]|_)+$") &
    !stringr::str_detect(x, "^[:digit:]")
}

#' convert lists, vectors, tibbles into \emph{bash-array}
#' @description Vectors or each column of tibbles will be interpretted as a bash-array when your output is read by bash. NA will be changed to "".
#' @param ... Lists, vectors, or tibbles. Elements with the same name will be overwritten by the last one.
#' @param option An option for declare function of bash. This argument is used for all arguments.
#' @export
as_bash_array <- function(..., option = "-a") {
  verify_scalar(option)
  option <- .as_character(option)
  dots <- rlang::list2(...)
  dots %>%
    rlist::list.flatten() %>%
    purrr::iwalk(~ assertthat::assert_that(is_bash_name(.y))) %>%
    purrr::imap_chr(~ stringr::str_c("declare ", option, " ", .y, "=(", convert_to_array(.x), ")")) %>%
    stringr::str_c(collapse = "\n")
}
#TODO always require names. e.g. hello = c(1,2);as_bash_array(hello = hello).

#' Directory requirements
#' @param cwd A logical. Whether set the directory where you run your code as current working directory. Otherwise, your home directory is set as current working directory.
#' @param out Path to write stdout
#' @param err Path to write stderr (option). If unspecified, \emph{out} argument will be used instead.
#' @export
directory_option <- function(cwd = FALSE, out = fs::path_home(), err = NA_character_) {
  verify_scalar(cwd, out, err)
  if (is.na(err)) err <- out
  out <- resource("-o", as.character(fs::path_abs(out)))
  err <- resource("-e", as.character(fs::path_abs(err)))
  cwd <- ifelse(!is.na(cwd) && cwd, "#$ -cwd", "")
  stringr::str_c(cwd, out, err, sep = "\n")
}
