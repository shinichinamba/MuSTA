args <- commandArgs(TRUE)

isocan <- args[1]# "../isocan"
genome_pkg <- args[2]
input_files <- args[3]
variant_colname <- args[4]
mut_colname <- args[5]
correlation_result <- args[6]
variant_filtered_out <- args[7] #output
n_workers <- as.integer(args[8]) #thread


suppressPackageStartupMessages({
  pkgname <- fs::path_file(genome_pkg)
  pkgload::load_all(genome_pkg, export_all = FALSE)
  genome <- eval(parse(text = paste0(pkgname, "::", pkgname)))
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(furrr)
  pkgload::load_all(isocan, export_all = FALSE)
})

plan(multicore, workers = n_workers)
input_files <- read_csv(input_files, col_types = cols(.default = col_character()))
correlation_result <- read_tsv(correlation_result, col_types = cols(.default = col_character()))

#functions----
tbl_0 <- 
  tibble(
    chr = character(),
    start = integer(), 
    end = integer(), 
    ref = character(),
    alt = character(),
    type = character(),
    .name_repair = "minimal"
  )

read_genomon_var_skip_na <- function(path) {
  if (is.na(path)) tbl_0
  else read_genomon_variant(path)
}

filter_variant_par_sample <- function(samplename, mut_path, variant_path, correlation) {
  message(paste0("processing: ", samplename))
  genomon_var <- read_genomon_var_skip_na(mut_path) %>% mutate(sample = samplename)
  retained_seq_id <- correlation %>% filter(sample == samplename) %>% `[[`("seq_id")
  pb_var <- read_tsv(variant_path, col_types = "ccciicciicc") %>% filter(seq_id %in% retained_seq_id) %>% filter_mono_exon() %>% mutate(sample = samplename)
  res <- future_map_dfr(split_df_group(pb_var, 200L, "seq_id"), isocan::filter_variant, genomon_var, genome, .progress = TRUE)
  message(paste0("      done: ", samplename))
  res
}

filter_mono_exon <- function(variant_df) {
  multi_exon_id <- 
    variant_df %>% 
    filter(type == "~") %>%
    distinct(seq_id) %>%
    `[[`("seq_id")
  variant_df %>%
    filter(seq_id %in% multi_exon_id)
}

split_df_group <- function(df, n, group_colname) {
  if (nrow(df) == 0L) df
  else {
    identifier <- unique(df[[group_colname]])
    len_id <- length(identifier)
    splt <- 
      tibble(
        idx = identifier,
        .name_repair = "minimal"
      ) %>%
      mutate(splitter = row_number() %>% `-`(1L) %>% `/`(len_id) %>% `*`(min(n, len_id)) %>% floor()) %>%
      right_join(
        tibble(idx = df[[group_colname]],
               .name_repair = "minimal"),
        by = "idx"
      ) %>%
      `[[`("splitter")
    split.data.frame(df, splt)
    }
}


message("filtering")
pmap_dfr(
  list(input_files$sample, input_files[[mut_colname]], input_files[[variant_colname]]),
  filter_variant_par_sample,
  correlation_result
) %>% 
  write_tsv(variant_filtered_out)
