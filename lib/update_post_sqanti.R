#parse args
args <- commandArgs(TRUE)

isocan <- args[1]

input_files <- args[2] #input
original_gtf <- args[3] #input
original_corr <- args[4]
original_fa <- args[5]
filter_res <- args[6] #sqanti result
gtf_out <- args[7] #update out
corr_out <- args[8]
fa_out <- args[9]
count_unique_out <- args[10] #used
existence_out <- args[11]

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(readr)
  library(purrr)
  pkgload::load_all(isocan, export_all = FALSE)
})

input_files <- read_csv(input_files, col_types = cols(.default = col_character()))

#filter
message("read sqanti filter result")
selected_isoform <- 
  filter_res %>% 
  read_tsv(
    col_types = cols(
      isoform = col_character(),
      strand = col_character(),
      structural_category = col_character(),
      associated_gene = col_character(),
      associated_transcript = col_character(),
      subcategory = col_character(),
      RTS_stage = col_logical(),
      all_canonical = col_character(),
      min_sample_cov = col_logical(),
      min_cov = col_logical(),
      sd_cov = col_logical(),
      n_indels = col_logical(),
      bite = col_logical(),
      FSM_class = col_character(),
      coding = col_character(),
      ORF_length = col_logical(),
      CDS_length = col_logical(),
      CDS_start = col_logical(),
      CDS_end = col_logical(),
      ML_classifier = col_logical()
    )
  ) %>% 
  filter(SQANTI_filter == "Isoform") %>% 
  `[[`("isoform")

#update gtf
message("update gtf")
updated_gtf <- 
  original_gtf %>% 
  isocan::read_gtf() %>% 
  filter(transcript_id %in% selected_isoform)
multi_exon_isoform <- 
  updated_gtf %>% 
  count(transcript_id) %>% 
  filter(n >= 2L) %>% 
  `[[`("transcript_id")
isocan::write_gtf(updated_gtf, gtf_out)

#update corr
message("update correlation.txt")
corr <- 
  original_corr %>% 
  read_tsv(col_types = cols(.default = col_character())) %>% 
  filter(transcript_id %in% selected_isoform)
write_tsv(corr, corr_out)

message("calculate FL count")

unique_corr <- 
  anti_join(
    corr, 
    corr %>% select(seq_id, sample) %>% mutate(dupli = duplicated(.)) %>% filter(dupli),
    by = c("seq_id", "sample")
  )

add_pbcount <- function(df, cluster_report_path, id_colname = "seq_id") {
  extract_count <- function(df, id_colname) {
    count_matrix <-
      df[[id_colname]] %>% 
      str_replace("^.+/f([:digit:]+)p([:digit:]+).+$", "\\1,\\2") %>% 
      str_split_fixed(",", 2)
    FLvec = as.integer(count_matrix[, 1])
    nFLvec = as.integer(count_matrix[, 2])
    if (any(is.na(FLvec)) || any(is.na(nFLvec))) stop("Failed to extract FL, nFL information from cluster_ids (this method cannot apply to Isoseq3). Please specify valid cruster_report.csv")
    df$FL <- FLvec
    df$nFL <- nFLvec
    df
  }
  if (is.null(cluster_report_path) || is.na(cluster_report_path)) {
    extract_count(df, id_colname)
  } else {
    count_df <- 
      cluster_report_path %>% 
      read_csv(col_types = cols(.default = col_character())) %>% 
      group_by(cluster_id) %>% 
      summarise(FL = sum(read_type == "FL"), nFL = sum(read_type == "nFL"))
    res <- left_join(df, count_df, by = `names<-`("cluster_id", id_colname))
    if (any(is.na(res$FL))) {
      warning("Not all cluster_ids are appeared in cluster_report.csv. Trying to extract FL, nFL information from cluster_ids...")
      extract_count(df, id_colname)
    } else {res}
  }
}
add_pbcount_map <- function(df, ifs = input_files) {
  df <- df %>% 
    mutate(sample = factor(sample, levels = ifs$sample))
  df %>% 
    group_split(sample) %>% 
    map2_dfr(
      ifs$cluster_report[sort(as.integer(unique(df$sample)))], 
      add_pbcount) %>% 
    mutate(sample = as.character(sample))
}
get_count <- function(corr_df) {
  corr_df %>% 
    add_pbcount_map() %>% 
    group_by(gene_id, transcript_id, sample) %>% 
    summarise(count_fl = sum(FL)) %>% 
    ungroup() %>% 
    spread(key = "sample", value = "count_fl", fill = 0L)
}


write_tsv(get_count(unique_corr), count_unique_out)

message("judge FL existence")
bind_rows(
  corr %>% filter(type == "intronic_match_include" & transcript_id %in% multi_exon_isoform),
  unique_corr %>% filter(!transcript_id %in% multi_exon_isoform)
) %>% 
  mutate(existence = TRUE) %>% 
  group_by(gene_id, transcript_id, sample) %>% 
  summarise(existence = any(existence)) %>% 
  ungroup() %>% 
  spread(key = "sample", value = "existence", fill = FALSE) %>% 
  write_tsv(existence_out)

#update fa
message("update fasta")
bs_fa <- Biostrings::readDNAStringSet(original_fa)
bs_fa[names(bs_fa) %in% selected_isoform] %>% Biostrings::writeXStringSet(fa_out)
message("done.")
