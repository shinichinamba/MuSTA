args <- commandArgs(TRUE)
input_files <- args[1]
chimera_colname <- args[2]
path_out <- args[3]

suppressPackageStartupMessages({
  library(readr)
})

input_files <- read_csv(input_files, col_types = cols(.default = col_character()))
res <- purrr::map2_dfr(input_files[[chimera_colname]], input_files$sample, 
                       ~ dplyr::mutate(read_tsv(..1, col_types = "ciiicciic"), sample = ..2))
write_tsv(res, path_out)
