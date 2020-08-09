#!/usr/bin/env Rscript

####functions####
# logged.error <- function(dir.path) {
#   # file.name <- paste0("pipeline_", format(Sys.time(), "%Y%m%d%H%M%S"), ".dump")
#   # dump.frames(dumpto = file.name, to.file = TRUE, include.GlobalEnv = TRUE)
#   print(rlang::last_error())
#   quit(save = "no", status = 1, runLast = TRUE)
# }
# .opt_error <- getOption("error")
# if (!interactive()) options(error = logged.error)

this_file <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    res <- sub(needle, "", cmdArgs[match])
  } else if (identical(cmdArgs, c("RStudio", "--interactive"))) {
    # Rstudio interactive mode
    res <- rstudioapi::getActiveDocumentContext()$path
  } else {
    # 'source'd via R console
    res <- sys.frames()[[1]]$ofile
  }
  if (res == "") "." else res
}

is_fun <- function(obj_name) is.function(eval(parse(text = obj_name)))

####check dependencies & load packages####
.this_dir <- dirname(this_file()) # not detected by ls()
source(paste0(.this_dir, "/src/version.R"))
source(paste0(.this_dir, "/src/check_dependencies.R"))

suppressPackageStartupMessages({library(argparse, quietly = TRUE)})
source(paste0(.this_dir, "/src/parse_arg.R"))

suppressPackageStartupMessages({library(rlang, quietly = TRUE)})
if (!identical(base::grep("/", args$prefix), integer(0))) abort("The value of '--prefix' must NOT contain '/'", "requirement error")

suppressPackageStartupMessages({
  library(stringr, quietly = TRUE)
  library(readr, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(tidyr, quietly = TRUE)
  library(purrr, quietly = TRUE)
  library(fs, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
  library(drake, quietly = TRUE)
})

obj <- ls()
obj <- setdiff(obj[!map_lgl(obj, is_fun)], "args")
rm(list = c(obj, "obj"))

####source files####
source(paste0(.this_dir, "/src/file_location.R"))  # written as relative paths from THIS file
pipeline_dir <- path_abs(.this_dir) # get the directory path of THIS file
lib <- map(lib, ~ path(pipeline_dir, .)) # as absolute paths
src <- map(src, ~ path(pipeline_dir, .)) # as absolute paths

source(src$utils.R)
load_jobwatcher_chr <- str_glue("suppressPackageStartupMessages(pkgload::load_all('{lib$jobwatcher}', export_all = FALSE, quiet = TRUE))")
eval(parse(text = load_jobwatcher_chr))
jobwatcher_mode <- get_jobwatcher_mode()

####print info####
cat(paste0("software: ", .name, "\n"))
cat(paste0("version: ", .version, "\n"))
cat(print_list(args))
cat("\n================================\n")
if (args$dry_run) {
  print(if (requireNamespace("sessioninfo", quietly = TRUE)) sessioninfo::session_info() else sessionInfo())
  cat("================================\n")
}
####parse config####
source(src$parse_config.R)

####make dir####
output_subdir <- 
  c("log", "script", "report", "genome", "samples", "merge", "result", "fusion",
  "merge/salmon", "merge/sqanti_filter_supplement") %>% 
  path(output_dir, .)
dir <- 
  map(output_subdir, force) %>% 
  setNames(str_remove(output_subdir, "^.*/"))

sample_dirs <- path(dir$samples, samples)

dir_log_path <- function(softname, o_dir = output_dir) path(o_dir, "log", softname)
log_categories <- c("file_management", "interleave", "lordec", "minimap2", "samtools", "merge", "salmon", "sqanti", "fusion", "others")
dir_logs <- 
  map(log_categories, dir_log_path) %>% 
  setNames(log_categories)

if (!args$dry_run) {
  walk(dir, dir_create, mode = "0775")
  walk(sample_dirs, dir_create, mode = "0775")
  walk(dir_logs, dir_create, mode = "0775")
  inform(paste0("The folder \'", output_dir, "\' has been created."))
}

####intermediate file location####
source(src$intermediate_files.R)

####tidy IO files####
source(src$update_input_files.R)
IO_summary_file <- path(dir$report, "IO.summary.csv")
if (!args$dry_run) readr::write_csv(input_files, IO_summary_file)

####Qscripts####
source(src$make_qscript.R)

####plan####
source(src$plan.R)
eval(parse(text = drake_txt))

if (args$force && !args$dry_run) clean(garbage_collection = TRUE, jobs = args$thread, purge = TRUE)

g <- drake_ggraph(drake_config(df_plan), targets_only = TRUE, mode = "out", label_nodes = TRUE) +
  theme_void() +
  coord_flip() + 
  ggtitle("MuSTA workflow")
suppressMessages({
  g <- g + scale_x_continuous(trans = "reverse")
})
g$layers <- map(g$layers, ~ {if (is(.x[["geom"]], "GeomText")) .x[["aes_params"]][["size"]] <- 3L;.x})

ggsave(if (!args$dry_run) path(dir$report, "plan_pre_run.pdf") else "plan_pre_run.pdf", g, width = 20, height = 30, units = "cm")

####check file existence####
qrecall_files <- c(input_files$long_read_hq, input_files$short_read_1, input_files$short_read_2, ref_gtf, genome_fa)
if (args$use_lq) qrecall_files <- c(qrecall_files, input_files$long_read_lq)
not_exist_files <- qrecall_files[!fs::file_exists(qrecall_files)]
if (length(not_exist_files) > 0L) abort(paste0("These files are not found: ", str_c(not_exist_files, collapse = ", ")), "requirement error")

####args$dry_run####
if (args$dry_run) {
  cat(print_list(file))
  cat("\n================================\n")
  cat(drake_txt)
  cat("\n")
  quit(save = "no", status = 0, runLast = FALSE)
}

####qrecall (SHIROKANE super computer only)####
# if (jobwatcher_mode == "hgc") {
#   write_and_qrecall(
#     qrecall_files,
#     path = path(dir$report, "recall_list"),
#     log_path = dir_logs[["file_management"]],
#     recursive = TRUE,
#     watch = TRUE
#   )
# }

####run pipeline####
future::plan(future::multiprocess, workers = args$thread)

make(df_plan, jobs = args$thread, parallelism = "future", prework = load_jobwatcher_chr)

cat("Completed.\n")

####save report####
g <- drake_ggraph(drake_config(df_plan), targets_only = TRUE, mode = "out", label_nodes = TRUE) +
  theme_void() + 
  coord_flip() + 
  ggtitle("MuSTA workflow")
suppressMessages({
  g <- g + scale_x_continuous(trans = "reverse")
})
g$layers <- map(g$layers, ~ {if (is(.x[["geom"]], "GeomText")) .x[["aes_params"]][["size"]] <- 3L;.x})

ggsave(path(dir$report, "plan_post_run.pdf"), g, width = 20, height = 30, units = "cm")
