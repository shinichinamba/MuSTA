p_load_chr <- function(pkg) {
  paste0(
    ifelse(base::requireNamespace(pkg, quietly = TRUE),
           "",
           paste0("if (!requireNamespace('", pkg, "', quietly = TRUE)) install.packages('", pkg, "', )\n")
    ),
    "library('", pkg, "')"
  )
}

p_load_github_chr <- function(pkg) {
  paste0(
    ifelse(
      base::requireNamespace("remotes", quietly = TRUE),
      "",
      "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes')\n"
    ),
    ifelse(
      base::requireNamespace(pkg, quietly = TRUE),
      "",
      paste0("if (!requireNamespace('", pkg, "', quietly = TRUE)) remotes::install_github('", pkg, "', )\n")
    ),
    "library('", pkg, "')"
  )
}

p_ruler <- function(x) paste0("#", x, "=================================")

pipeline_preset <- function(pipe_name, pipe_dir, n_parallel, pipe_memory) {
  pipe_dir <-
    fs::path_abs(pipe_dir)
  log_qrecall <-
    fs::path(pipe_dir, "log", "file_management")
  # file_qrecall <-
  #   fs::path(pipe_dir, "config", "qrecall.txt")
  dir_output <-
    fs::path(pipe_dir, "output")
  dir_script <- 
    fs::path(pipe_dir, "script")
  log_output <-
    fs::path(pipe_dir, "log", "hello")
  file_helloworld <-
    fs::path(dir_output, "helloworld", ext = "txt")
  log_pipeline <-
    fs::path(pipe_dir, "log")

  pipe_Rfile <-
    stringr::str_glue(
      '
      {p_ruler("attach libraries")}
      {p_load_chr("ggplot2")}
      {p_load_chr("drake")}
      {p_load_github_chr("jobwatcher")}
      
      ',{p_ruler("declare variables")},'
      Hello <- c("H", "e", "l", "l", "o")
      
      ',{p_ruler("parse config")},'
      # If you use config yaml files, we recommend fascinating {{config}} and {{rlist}} packages to parse them.
      
      ',{p_ruler("make directories")},'
      fs::dir_create("{dir_output}")
      fs::dir_create("{log_output}")
      
      ',{p_ruler("summarize in/out paths")},'
      your_qrecall_objects_1 <- "/archive/data/hgc1043/snamba/.jobwatch/for_qrecall.txt"
      your_qrecall_objects_2 <- "/archive/data/hgc1043/snamba/.jobwatch/for_qrecall_2.txt"
      
      ',{p_ruler("intermediate file paths")},'
      
      ',{p_ruler("make functions of qscript files")},'
      dir_opt <- directory_option(
        out = "{log_output}",
        err = "{log_output}"
      )
      pl_makefile <- qsub_function(
        "touch {file_helloworld}",
        script_path = "makefile",
        script_dir = "{dir_script}",
        directory = "dir_opt",
        other_req = "#$ -l os7"
      )
      pl_hello <- qsub_function(
        as_bash_array(Hello = Hello),
        "echo ${{Hello[$SGE_TASK_ID]}} >> {file_helloworld}",
        script_path = "hello",
        script_dir = "{dir_script}",
        directory = dir_opt,
        other_req = "#$ -l os7"
      )
      pl_world <- qsub_function(
        "echo World >> {file_helloworld}",
        script_path = "world",
        script_dir = "{dir_script}",
        directory = dir_opt,
        other_req = "#$ -l os7"
      )
      pl_helloworld <- qsub_function(
        "echo HelloWorld >> {file_helloworld}",
        script_path = "helloworld",
        script_dir = "{dir_script}",
        directory = dir_opt,
        other_req = "#$ -l os7"
      )
      
      ',{p_ruler("qrecall")},'
      job_recall <- write_and_qrecall(
        your_qrecall_objects_1, your_qrecall_objects_2,
        path = "{log_pipeline}",
        log_path = "{log_qrecall}",
        watch = TRUE
      )
      
      ',{p_ruler("make pipeline")},'
      pipeline <- 
        drake::drake_plan(
          makefile = pl_makefile(),
          hello = pl_hello(makefile),
          world = pl_world(makefile),
          helloworld = pl_helloworld(c(hello, world))
        )
        
      ',{p_ruler("run pipeline")},'
      ggsave_pipeline(pipeline, "{fs::path(pipe_dir, "log", "pipeline_pre.pdf")}", width = 30, height = 10)
      drake::make(pipeline, jobs = {n_parallel}L)
      ggsave_pipeline(pipeline, "{fs::path(pipe_dir, "log", "pipeline_post.pdf")}", width = 30, height = 10)
      ', .sep = "\n"
    )

  pipe_qsubfile <-
    make_qsubfile(
      paste0("Rscript ", fs::path(pipe_dir, pipe_name), ".R"),
      name = pipe_name,
      parallel = parallel_option(slot = n_parallel, memory = pipe_memory, ljob = TRUE),
      arrayjob = arrayjob_option(),
      directory = directory_option(out = log_pipeline, err = log_pipeline),
      use_bash_profile = TRUE,
      other_req = "#$ -l os7"
    )
  list(pipe_qsubfile, pipe_Rfile)
}

#' save a plan image
#'
#' @description wrapper of \code{drake::\link[drake]{drake_ggraph}} and \code{ggplot2::\link[ggplot2]{ggsave}}
#' @param plan A plan made with \code{drake::\link[drake]{make}}
#' @param path A path to write image
#' @param ... Additional arguments for \code{ggplot2::\link[ggplot2]{ggsave}}
#' @export
ggsave_pipeline <- function(plan, path, ...){
  if (!requireNamespace("ggplot2", quietly = TRUE)) rlang::abort("ggplot2 is required. Please install the package.")
  g <- plan %>% drake::drake_config() %>% drake::drake_ggraph(targets_only = TRUE, mode = "out", label_nodes = TRUE) +
    ggplot2::theme_void() + 
    ggplot2::coord_flip()
  suppressMessages(
    g <- g + ggplot2::scale_x_continuous(trans = "reverse")
  )
  ggplot2::ggsave(path, g, ...)

}

#' build a pre-build pipeline
#'
#' @param pipe_name A character. Your pipeline name.
#' @param pipe_dir A directory path where you intend to build a pipeline.
#'  If it does not exist, it will be made recursively.
#' @param n_parallel A integer. The number you run elements of your pipeline simultaneously.
#' @param pipe_memory A numeric. Request of memory size par slot for pipeline manager.
#' @param force Alogical. Whether to run this function when \emph{pipe_dir} exists and \emph{make_subdir} is TRUE.
#' @param make_subdir A logical. Whether to make minimum subdirectories under \emph{pipe_dir}
#' @export
build_pipeline <- function(pipe_name, pipe_dir, n_parallel = 2L, pipe_memory = 3L, force = FALSE, make_subdir = TRUE){
  file_list <- pipeline_preset(pipe_name, pipe_dir, n_parallel, pipe_memory)
  dir_exist <- fs::dir_exists(pipe_dir)
  if ((!force) && make_subdir && dir_exist) {
    rlang::abort(paste0(pipe_dir, "exists. If you'd like to continue, set force = TRUE."))
  }
  if (!dir_exist) {
    fs::dir_create(pipe_dir)
    rlang::inform(done("Directory '", crayon::cyan(pipe_dir), "' has been created."))
  }
  if (make_subdir) {
    purrr::map(list("log", "script", "config"), ~ fs::path(pipe_dir, .x)) %>%
      purrr::walk(fs::dir_create) %>%
      purrr::walk(~ rlang::inform(done("Directory '", crayon::cyan(.x), "' has been created.")))
  }
  path_list <- purrr::map(c("sh", "R"), ~ fs::path(pipe_dir, pipe_name, ext = .x))
  purrr::walk2(path_list, file_list, ~ write(.y, .x, append = FALSE)) %>%
    purrr::walk(~ rlang::inform(done("File '", crayon::cyan(.x), "' has been written.")))
  rlang::inform(todo("Please edit '", crayon::cyan(path_list[[2]]), "' for your own pipeline."))
  rlang::inform(todo("Then, run ", crayon::green(paste0('qsub("', path_list[[1]], '")'))))
}

