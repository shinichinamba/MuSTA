args <- commandArgs(TRUE)

isocan <- args[1]# "~/long_read/isocan"

input_files <- args[2]
chimera <- args[3]
ref_gtf <- args[4]
chimera_class <- args[5]

summary_out <- args[6]
gtf_out <- args[7]
sv_out <- args[8]

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  pkgload::load_all(isocan, export_all = FALSE)
})

input_files <- read_csv(input_files, col_types = cols(.default = col_character()))

# fun ----
separate_seqid_idx <- function(df) separate(df, seq_id, c("seq_id", "idx"), sep = "___") %>% mutate_at(vars("idx"), as.integer)
# add_count_cols <- function(df) {
#   df %>%
#     mutate(.count_col = seq_id %>% str_remove("^.+\\|.+/f") %>% str_remove("/[:digit:]+$")) %>% 
#     separate(.count_col, c("full", "partial"), sep = "p") %>% 
#     mutate_at(vars("full", "partial"), as.integer)
# }
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

# input ----
chimera <- readRDS(chimera)
ref_gtf <- isocan::read_gtf(ref_gtf)

chimera_class <- 
  read_tsv(
    chimera_class,
    col_types = cols(
      isoform = col_character(),
      length = col_integer(),
      exons = col_integer(),
      associated_gene = col_character(),
      associated_transcript = col_character(),
      structural_category = col_character(),
      subcategory = col_character(),
      .default = col_skip()
    )
  ) %>% 
  rename(seq_id = "isoform") %>% 
  separate_seqid_idx()

# main ----
chimera_sv <- 
  bind_rows(
    chimera$sv_uni_match$sv,
    chimera$sv_multi_match$sv,
    chimera$sv_partial_match$sv
  )
chimera_gtf <- 
  bind_rows(
    chimera$sv_uni_match$gtf,
    chimera$sv_multi_match$gtf,
    chimera$sv_partial_match$gtf
  )
chimera_id <- 
  c(
    chimera$sv_uni_match$id,
    chimera$sv_multi_match$id,
    chimera$sv_partial_match$id
  )

chimera_block_info <- 
  if (all(c("gene_id", "gene_name") %in% colnames(ref_gtf))) {
    left_join(chimera_class, 
              ref_gtf %>% distinct(gene_id, gene_name),
              by = c("associated_gene" = "gene_id")
    ) %>% arrange(seq_id, idx)
  } else chimera_class %>% arrange(seq_id, idx)

chimera_gene_combination <-
  chimera_block_info %>% 
  replace_na(list(gene_name = "novel")) %>% 
  group_by(seq_id) %>% 
  summarise(gene_combination = str_c(gene_name, collapse = "_")) %>% 
  ungroup() %>% 
  left_join(
    chimera_sv %>% 
      distinct(seq_id, sample),
    by = "seq_id"
  ) %>% 
  arrange(gene_combination, sample, seq_id)


chimera_summary <- 
  chimera_gene_combination %>% 
  filter(seq_id %in% chimera_id) %>% 
  add_pbcount_map() %>% 
  dplyr::rename(full = "FL", partial = "nFL") %>% 
  mutate(n_SV = str_count(gene_combination, "_")) %>% 
  arrange(desc(n_SV))

chimera_gtf_summary <- 
  chimera_gtf %>% 
  inner_join(chimera_summary %>% select(seq_id, gene_id = "gene_combination"), by = "seq_id") %>% 
  arrange(sample, seq_id, exon_number) %>% 
  rename(transcript_id = "seq_id", feature = "type", seqname = "chr", source = "sample")

chimera_sv_summary <- 
  chimera_sv %>% 
  arrange(sample, seq_id, sv_number) %>% 
  rename(up_pos_genomic = "Pos_up", down_pos_genomic = "Pos_down", inserted_seq_genomic = "Inserted_Seq", sv_crossing_type = "sv_type", sv_type = "Type")

write_tsv(chimera_summary, summary_out)
write_gtf(chimera_gtf_summary, gtf_out)
write_tsv(chimera_sv_summary, sv_out)
