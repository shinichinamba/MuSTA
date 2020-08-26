args <- commandArgs(TRUE)

isocan <- args[1]# "../isocan"
merge_iso_commonpart.R <- args[2]

paf_hq <- args[3]# "hq.lordec.qfilt.paf"
paf_lq <- args[4]# "lq.lordec.qfilt.paf" | "NULL"

corr_out <- args[5]
gtf_out <- args[6]
variant_out <- args[7]
paf_fusion_out <- args[8]
n_worker <- as.integer(args[9])
mono_multi_merge <- as.logical(args[10])
collapse_mono_range <- as.logical(args[11])

source(merge_iso_commonpart.R)

message("read inputs")
paf <- read_paf(paf_hq, .query_position = TRUE)
if (paf_lq != "NULL") paf <- bind_rows(paf, read_paf(paf_lq, .query_position = TRUE)) 
paf_fusion <- paf %>% filter(seq_id %in% paf$seq_id[duplicated(paf$seq_id)])
paf <- anti_join(paf, paf_fusion, by = colnames(paf))

#exclude fusion
message("extract variants")
placeholder = ""
variant <- 
  paf %>% 
  mutate(sample = placeholder) %>% 
  split_df(n_worker) %>% 
  future_map_dfr(~ {
    suppressPackageStartupMessages(pkgload::load_all(isocan, export_all = FALSE, quiet = TRUE))
    isocan::convert_paf_variant(..1, sample_name = placeholder)
  })

message("extract gtf")
gtf <- 
  paf %>% 
  mutate(sample = placeholder) %>% 
  isocan::retrieve_gtf(sample_name = placeholder, paf_variants = variant) %>% 
  select(seq_id, strand, seqname = "chr", start, end)
gtf_nested <- 
  gtf %>% 
  isocan::gather_exon() %>% 
  isocan::separate_ends()
message("separate multi-/mono-exon isoforms")

#multi-exon: merge tx whose splicing junctions are completely same.
message("multi-exon nest")
gtf_nested_multi <- 
  gtf_nested %>% 
  filter(exon_structure != "-") %>% 
  group_by(seqname, strand, exon_structure) %>% 
  summarise(seq_id = list(seq_id),
            tx_start = min(tx_start),
            tx_end = max(tx_end)) %>% 
  ungroup() %>% 
  mutate(index = row_number())

#mono-exon: included tx (not overlapping tx) are merged
message("mono-exon nest")
gtf_nested_mono <- 
  gtf_nested %>% 
  filter(exon_structure == "-") %>% 
  group_by(seqname, strand, tx_start, tx_end, exon_structure) %>%
  summarise(seq_id = list(seq_id)) %>%
  ungroup() %>%
  mutate(index = -row_number())

#included in other tx?

include_5prime_multiex_1 <- function(gtf_nested) {
  match_3prime_list <- match_3prime_exonstart_multiex(gtf_nested)
  future_map_dfr(match_3prime_list, ~{
    gtf_tmp <- gtf_nested %>% filter(index %in% .x)
    include_table_nest_core_multiex(gtf_tmp, gtf_tmp)
  },
  .progress = TRUE) %>% 
    unnest() %>% 
    df_query_subject() %>% 
    filter(query != subject) %>% 
    future_include_table_edge_multiex(gtf_nested, gtf_nested) %>% 
    filter(!subject %in% query)
}

match_3prime_exonstart_multiex <- function(gtf_nested) {
  match_list_minus <- 
    gtf_nested %>% 
    filter(strand == "-") %>% 
    mutate(last_exonstart = exon_structure %>% str_extract("^-[:digit:]+ ")) %>% 
    select(seqname, strand, last_exonstart, index) %>% 
    group_by(seqname, strand, last_exonstart) %>% 
    summarise(index = list(index)) %>% 
    ungroup() %>% 
    `[[`("index")
  match_list_plus <- 
    gtf_nested %>% 
    filter(strand == "+") %>% 
    mutate(last_exonstart = exon_structure %>% str_extract(" [:digit:]+-$")) %>% 
    select(seqname, strand, last_exonstart, index) %>% 
    group_by(seqname, strand, last_exonstart) %>% 
    summarise(index = list(index)) %>% 
    ungroup() %>% 
    `[[`("index")
  match_list <- c(match_list_plus, match_list_minus)
  match_list[map_int(match_list, length) > 1L]
}

include_nest_multiex_1 <- function(gtf_nested) {
  gtf_nested %>% 
    group_split(seqname, strand) %>% 
    future_map_dfr(~ include_table_nest_core_multiex(.x, .x), .progress = TRUE) %>% 
    unnest() %>% 
    df_query_subject() %>% 
    filter(query != subject) %>% 
    future_include_table_edge_multiex(gtf_nested, gtf_nested) %>% 
    filter(!subject %in% query)
}

include_table_nest_core_multiex <- function(gtf_nest_1, gtf_nest_2) {
  tibble(
    query = map(gtf_nest_2[["exon_structure"]], ~ {
      gtf_nest_1[["index"]][stringi::stri_detect_fixed(.x, gtf_nest_1[["exon_structure"]])]
    }),
    subject = gtf_nest_2[["index"]]
  ) 
}

future_include_table_edge_multiex <- function(include_table, gtf_nest_1, gtf_nest_2) {
  if (nrow(include_table) == 0L) {
    return(include_table)
  } else {
    gtf_nest_2_mod <- 
      gtf_nest_2 %>% 
      mutate(exon_structure_2 = paste0(tx_start, exon_structure, tx_end)) %>% 
      select(index, exon_structure_2) %>% 
      filter(index %in% include_table$subject)
    gtf_nest_1_mod <- 
      gtf_nest_1 %>% 
      select(index, tx_start, tx_end, exon_structure) %>% 
      filter(index %in% include_table$query)
    main_table <- 
      include_table %>% 
      left_join(gtf_nest_1_mod, by = c("query" = "index")) %>% 
      left_join(gtf_nest_2_mod, by = c("subject" = "index"))
    filter_by_edge <- function(main_table) {
      edge_table <- 
        map2(main_table$exon_structure_2, main_table$exon_structure, stringi::stri_split_fixed, n = 2L) %>% 
        simplify_all() %>% 
        transpose(.names = c("smaller", "larger")) %>% 
        as_tibble() %>% 
        unnest() %>% 
        mutate(smaller = str_extract(smaller, "[:digit:]+$"),
               larger = str_extract(larger, "^[:digit:]+")) %>% 
        mutate_all(as.integer)
      bind_cols(
        main_table,
        edge_table
      ) %>% 
        filter(smaller <= tx_start & tx_end <= larger) %>% 
        select(query, subject)
    }
    main_table %>% 
      split_df(n_worker) %>% 
      future_map_dfr(filter_by_edge, .progress = TRUE)
  }
}
  
include_nest_mono_any_2 <- function(gtf_nested_mono, gtf_nested_any, index_colname = "index") {
  combination_table <- 
    inner_join(
      gtf_nested_mono %>% dplyr::rename(query := !!sym(index_colname)) %>% select(seqname, strand, query),
      gtf_nested_any %>% dplyr::rename(subject := !!sym(index_colname)) %>% select(seqname, strand, subject),
      by = c("seqname", "strand")
    )
  if (nrow(combination_table) == 0L) {
    combination_table
  }else {
    gtf_long_any <- 
      gtf_nested_any %>% 
      split_df(n_worker) %>% 
      future_map_dfr(gtf_longer, .progress = TRUE) %>% 
      get_lastexon(index_colname)
    combination_table %>% 
      split_df(n_worker) %>% 
      future_map_dfr(filter_by_include, gtf_nested_mono, gtf_long_any, index_colname, .progress = TRUE)
  }
}

include_nest_mono_1 <- function(gtf_nested_mono, index_colname = "index") {
  gtf_long_any <- 
    gtf_nested_mono %>% 
    split_df(n_worker) %>% 
    future_map_dfr(gtf_longer, .progress = TRUE)
  combination_table <- 
    inner_join(
      gtf_nested_mono %>% dplyr::rename(query := !!sym(index_colname)) %>% select(seqname, strand, query),
      gtf_nested_mono %>% dplyr::rename(subject := !!sym(index_colname)) %>% select(seqname, strand, subject),
      by = c("seqname", "strand")
    ) %>% 
    filter(query != subject)
  if (nrow(combination_table) == 0L) {
    combination_table
  }else {
    combination_table %>% 
      split_df(n_worker) %>% 
      future_map_dfr(filter_by_include, gtf_nested_mono, gtf_long_any, index_colname, .progress = TRUE) %>% 
      filter(!subject %in% query)
  }
}

filter_by_include <- function(combination_table, gtf_nested_mono, gtf_long_any, index_colname = "index") {
  if (nrow(combination_table) == 0L) {
    return(combination_table)
  } else {
    combination_table %>% 
      left_join(gtf_nested_mono %>% select(index_colname, "tx_start", "tx_end"), by = c("query" = index_colname)) %>% 
      left_join(gtf_long_any %>% select(index_colname, "start", "end"), by = c("subject" = index_colname)) %>% 
      filter(start <= tx_start & tx_end <= end) %>% 
      select(query, subject)
  }
}


certify_seq_id <- function(df) {# In case unnest drop columns with empty lists
  if (!"seq_id" %in% colnames(df)) df %>% mutate(seq_id = character())
  else df
}


df_dammy_res <- tibble(index = integer(), seq_id = character(), query = integer(), .name_repair = "minimal")

#multi vs multi
message("multi-exon merge")
include_multi <- 
  gtf_nested_multi %>% 
  include_5prime_multiex_1() %>% 
  select(query, index = subject) %>% 
  left_join(gtf_nested_multi %>% select(index, seq_id), by = c("query" = "index")) %>% 
  unnest() %>% certify_seq_id()

intronic_match_include_multi <- 
  gtf_nested_multi %>%
  select(seq_id, index) %>% 
  filter(!index %in% include_multi$query) %>% 
  unnest() %>% certify_seq_id()

gtf_nested_multi_collapsed <- gtf_nested_multi %>% filter(!index %in% include_multi$query)

#mono vs multi
message("mono-exon merge")
#leave only non-overlapped tx
include_mono_multi <- 
  if (mono_multi_merge == TRUE) {
    include_nest_mono_any_2(gtf_nested_mono, gtf_nested_multi_collapsed) %>% 
      select(query, index = "subject") %>% 
      left_join(gtf_nested_mono %>% select(index, seq_id), by = c("query" = "index")) %>% 
      unnest() %>% certify_seq_id()
  } else df_dammy_res

gtf_nested_mono_include <- gtf_nested_mono %>% filter(index %in% include_mono_multi$query)
gtf_nested_mono_noninclude <- gtf_nested_mono %>% filter(!index %in% include_mono_multi$query)

#mono vs mono
if (collapse_mono_range == TRUE) {
  collapse_mono_range_res <- collapse_mono(gtf_nested_mono_noninclude, "index")
  match_mono <- 
    collapse_mono_range_res[[1]] %>% 
    left_join(gtf_nested_mono %>% select(index, seq_id), by = c("old_idx" = "index")) %>% 
    unnest() %>% certify_seq_id() %>% 
    select(index = new_idx, seq_id)
  include_mono <- 
    collapse_mono_range_res[[2]] %>% 
    left_join(gtf_nested_mono %>% select(index, seq_id), by = c("old_idx" = "index")) %>% 
    unnest() %>% certify_seq_id() %>% 
    select(index = new_idx, seq_id)
  gtf_nested_mono_collapsed <- 
    collapse_mono_range_res[[3]] %>% 
    dplyr::rename(index = "collapsed_idx")
} else {
  include_mono <- 
    include_nest_mono_1(gtf_nested_mono_noninclude) %>% 
    select(query, index = "subject") %>% 
    left_join(gtf_nested_mono %>% select(index, seq_id), by = c("query" = "index")) %>% 
    unnest() %>% certify_seq_id()
  gtf_nested_mono_collapsed <- gtf_nested_mono_noninclude %>% filter(!index %in% include_mono$query)
  match_mono <- 
    gtf_nested_mono_collapsed %>% 
    select(index) %>%
    left_join(gtf_nested_mono %>% select(index, seq_id), by = "index") %>%
    unnest() %>% certify_seq_id()
}

#mono(included in multi) vs mono
include_mono_mono <- 
  include_nest_mono_any_2(gtf_nested_mono_include, gtf_nested_mono_collapsed) %>% 
  dplyr::rename(index = "subject") %>% 
  left_join(gtf_nested_mono %>% select(index, seq_id), by = c("query" = "index")) %>% 
  unnest() %>% certify_seq_id()

#summarize correlation
message("summarize")
tidy_corr <- function(df, type = "include", info = "") {
  message(paste0("  ", info))
  if (nrow(df) == 0L) df <- df %>% mutate(seq_id = character())
  df %>% select(index, seq_id) %>% mutate(type = type)
}
corr <- 
  bind_rows(
    include_multi %>% tidy_corr(info = "multi"),
    include_mono %>% tidy_corr(info = "mono"),
    include_mono_multi %>% tidy_corr(info = "mono-multi"),
    include_mono_mono %>% tidy_corr(info = "mono-mono"),
    intronic_match_include_multi %>% tidy_corr("intronic_match_include", info = "mch-inc-multi"),
    match_mono %>% tidy_corr("monoex_match", info = "mch-mono")
  )

#result
message("get gtf")
gtf_collapsed <- 
  bind_rows(gtf_nested_multi_collapsed, gtf_nested_mono_collapsed) %>% 
  select(index, seqname, strand, tx_start, tx_end, exon_structure) %>% 
  split_df(n_worker) %>% 
  future_map_dfr(gtf_longer)
ex_grl <- 
  gtf_collapsed %>% 
  GRanges() %>% 
  split(.$index)
tx_gr <- range(ex_grl)

split_rows <- function(x) split.default(x, seq_along(x)) %>% GRangesList()

include_2_mod <- function(grl_x, grl_y) {
  include_table(grl_x, grl_y) %>% 
    retrieve_id_mod(grl_x, grl_y)
}
include_table <- function(grl_x, grl_y) {
  GenomicRanges::findOverlaps(grl_x, grl_y, type = "within", select = "all") %>% 
    as_tibble()
}
retrieve_id_mod <- function(df, grl_x, grl_y) {
  tibble(
    super_index = names(grl_y)[df$subjectHits],
    index = names(grl_x)[df$queryHits],
    .name_repair = "minimal"
  ) %>% 
    mutate_all(as.integer)
}
message("PB style ID")
gene_index_df <- 
  include_2_mod(
    tx_gr,
    tx_gr %>% unlist() %>% GenomicRanges::reduce() %>% split_rows()
  )
gene_index_df <- 
  gene_index_df %>% 
  mutate(gene_id = group_indices(gene_index_df, super_index) %>% paste0("PB.", .)) %>% 
  select(index, gene_id) %>% 
  group_by(gene_id) %>% 
  mutate(transcript_id = row_number()) %>% 
  ungroup() %>% 
  mutate(transcript_id = paste0(gene_id, ".", transcript_id))

corr_result <- 
  full_join(corr, gene_index_df, by = "index") %>% 
  select(gene_id, transcript_id, seq_id, type)
gtf_result <- 
  gtf_collapsed %>% 
  select(index, seqname, start, end, strand) %>% 
  full_join(gene_index_df, by = "index") %>% 
  select(-index) %>% 
  mutate(feature = "exon")

#write out
message("write out")
write_tsv(corr_result, corr_out)
isocan::write_gtf(gtf_result, gtf_out)
write_tsv(variant, variant_out)
write_tsv(paf_fusion, paf_fusion_out)
