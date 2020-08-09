args <- commandArgs(TRUE)

args_IO_summary <- args[1]
args_IO_colname <- args[2]
args_salmon_colname <- args[3]
args_outname <- args[4]

suppressPackageStartupMessages({
  library(rlang)
  library(dplyr)
  library(readr)
})

IO_summary <- 
  read_csv(args_IO_summary, col_types = cols(.default = col_character())) %>% 
  select(c("sample", args_IO_colname)) %>%
  rename(chain_salmon_iso := !!sym(args_IO_colname))

bind_cols(
  read_tsv(IO_summary$chain_salmon_iso[1L], col_types = cols_only(Name = col_character())),
  purrr::map2_dfc(
    IO_summary$sample,
    IO_summary$chain_salmon_iso,
    ~ read_tsv(.y, col_types = do.call(cols_only, list2(!!sym(args_salmon_colname) := col_character()))) %>% `colnames<-`(.x)
  )
) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  tibble::column_to_rownames("Name") %>% 
  write.table(args_outname, sep = "\t", quote = F)