#common functions for merge_iso_intersamples.R and merge_iso_intrasamples.R
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rlang)
  library(readr)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  pkgload::load_all(isocan, export_all = FALSE)
})

plan(multicore, workers = n_worker)
options(future.globals.maxSize = 3000L * (1024L^2L)) #allow 3GB. for intermerge


split_df <- function(df, n) if (nrow(df) == 0L) df else split.data.frame(df, floor((seq_len(nrow(df)) - 1) * min(n, nrow(df)) / nrow(df)))
split_df_group <- function(df, n, group_colname) {
  if (nrow(df) == 0L) df 
  else {
    identifier <- unique(df[[group_colname]])
    splt <- 
      tibble(idx = identifier, .name_repair = "minimal") %>% 
      mutate(splitter = row_number() %>% `-`(1L) %>% `/`(length(identifier)) %>% `*`(min(n, length(identifier))) %>% floor()) %>% 
      right_join(
        tibble(idx = df[[group_colname]], .name_repair = "minimal"),
        by = "idx"
      ) %>% 
      `[[`("splitter")
    split.data.frame(df, splt)
  }
} 

gtf_longer <- function(gtf_nested) {
  gtf_nested <- gtf_nested %>% mutate(.index = row_number())
  gtf_nested %>% 
    mutate(exon_structure = paste0(tx_start, exon_structure, tx_end) %>%
             stringi::stri_split_fixed(" ")) %>% 
    select(.index, exon_structure) %>% 
    unnest() %>% 
    separate(col = "exon_structure", into = c("start", "end"), sep = "-", convert = TRUE) %>% 
    left_join(
      gtf_nested %>% select(-c(tx_start, exon_structure, tx_end)),
      by = ".index"
    ) %>% 
    select(-.index)
}

get_lastexon <- function(gtf_longer, identifier_colname) {
  get_all_first <- function(df) {
    if (nrow(df) > 0L) {
      unnest(summarise_all(df, dplyr::first)) 
      # dplyr > 1.0.0 automatically unnest list columns...
    } else {
      unnest(df)
    }
  }
  bind_rows(
    gtf_longer %>% 
      filter(strand == "-") %>% 
      group_by(!!sym(identifier_colname)) %>% 
      arrange(start) %>% 
      get_all_first() %>% 
      ungroup(),
    gtf_longer %>% 
      filter(strand == "+") %>% 
      group_by(!!sym(identifier_colname)) %>% 
      arrange(desc(start)) %>% 
      get_all_first() %>% 
      ungroup()
  )
}

df_query_subject <- function(df) {
  if (nrow(df) == 0 && ncol(df) == 0) tibble(query = integer(), subject = integer(), .name_repair = "minimal")
  else df
}

collapse_mono <- function(gtf_mono, index_colname) {
  ex_gr <- GRanges(gtf_mono)
  ex_grl <- split(ex_gr, gtf_mono[[index_colname]])
  split_rows <- function(x) split.default(x, seq_along(x)) %>% GRangesList()
  ex_grl_collapsed <- ex_gr %>% GenomicRanges::reduce() %>% split_rows()
  gtf_mono_collapsed <- 
    ex_grl_collapsed %>%
    as_tibble() %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate(group = -group,
           exon_structure = "-") %>% 
    select(collapsed_idx = "group", seqname = "seqnames", strand, tx_start = "start", tx_end = "end", exon_structure)
  gene_index_df_match <- 
    inner_join(
      gtf_mono %>% select(!!index_colname, seqname, strand, tx_start, tx_end),
      gtf_mono_collapsed,
      by = c("seqname", "strand", "tx_start", "tx_end")) %>% 
    select(!!index_colname, collapsed_idx) %>% 
    `colnames<-`(c("old_idx", "new_idx"))
  df <- 
    GenomicRanges::findOverlaps(
      ex_grl, 
      ex_grl_collapsed,
      type = "within",
      select = "all") %>% 
    as_tibble()
  gene_index_df <- 
    tibble(
      new_idx = -df$subjectHits,
      index = names(ex_grl)[df$queryHits],
      .name_repair = "minimal"
    ) %>% 
    mutate_all(as.integer)
  gene_index_df_include <- anti_join(gene_index_df, gene_index_df_match, by = c("old_idx", "new_idx"))
  stopifnot(nrow(anti_join(gene_index_df_match, gene_index_df, by = c("old_idx", "new_idx"))) == 0L)
  list(gene_index_df_match, gene_index_df_include, gtf_mono_collapsed)
}
