args <- commandArgs(TRUE)

input_files <- args[1]
final_corr <- args[2]
variant_file_colname <- args[3]
original_multiex_out <- args[4]
sample_name_colname <- "sample"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
})

final_corr <- read_tsv(final_corr, col_types = cols(.default = col_character()))
input_files <- read_csv(input_files, col_types = cols(.default = col_character()))

isoseq_tx_multiex <- 
  map2_dfr(input_files[[variant_file_colname]], input_files[[sample_name_colname]], ~ {
    read_tsv(..1, col_types = "ccciicciicc") %>% 
      filter(type == "~") %>% 
      distinct(seq_id, chr, tx_strand, tx_start, tx_end) %>% 
      mutate(sample = ..2)
  }) %>% 
  mutate(tx_start = tx_start + 1L) %>% #0-based -> 1-based
  rename(strand = "tx_strand", start = "tx_start", end = "tx_end")

inner_join(
  final_corr %>% select(gene_id, transcript_id, sample, seq_id, type),
  isoseq_tx_multiex,
  by = c("seq_id", "sample")
) %>% 
  select(gene_id, transcript_id, sample, seq_id, type, chr, strand, start, end) %>% 
  write_tsv(original_multiex_out)
