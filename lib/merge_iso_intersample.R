args <- commandArgs(TRUE)

isocan <- args[1]# "../isocan"
merge_iso_commonpart.R <- args[2]

input_files <- args[3]# "hgc_results/pipeline/report/IO.summary.csv"
intra_merge_correlation_colname <- args[4]# "intra_merge_correlation"
intra_merge_gtf_colname <- args[5]# "intra_merge_gtf"

corr_out <- args[6]
gtf_out <- args[7]
count_out_colname <- args[8]
n_worker <- as.integer(args[9])
mono_multi_merge <- as.logical(args[10])
collapse_mono_range <- as.logical(args[11])

source(merge_iso_commonpart.R)

message("===read inputs===")
input_files <- read_csv(input_files, col_types = cols(.default = col_character()))
message("gtf")
gtf_sample <- map2_dfr(input_files[[intra_merge_gtf_colname]], input_files$sample, ~ read_gtf(..1) %>% mutate(sample = ..2))
message("correlation files")
corr_sample <- map2_dfr(
  input_files[[intra_merge_correlation_colname]],
  input_files$sample,
  ~ {
    read_tsv(
      ..1,
      col_types = cols(
        gene_id = col_character(),
        transcript_id = col_character(),
        seq_id = col_character(),
        type = col_character()
      )
    ) %>% mutate(sample = ..2)
  })

# functions----
message("===declare functions===")

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

collapse_sj_match <- function(gtf_nested_multi, add_minus = FALSE) {
  gtf_nested_res <-
    gtf_nested_multi %>% 
    group_by(seqname, strand, exon_structure) %>% 
    summarise(tx_start = min(tx_start),
              tx_end = max(tx_end),
              index = list(index)) %>% 
    ungroup() %>% 
    mutate(new_idx = row_number())
  if (add_minus) gtf_nested_res <- gtf_nested_res %>% mutate(new_idx = new_idx * -1L)
  list(gtf = gtf_nested_res %>% select(-index),
       tbl = gtf_nested_res %>% select(new_idx, index) %>% unnest())
}

collapse_tx_match <- function(gtf_nested, add_minus = FALSE) {
  gtf_nested_res <-
    gtf_nested %>% 
    group_by(seqname, strand, exon_structure, tx_start, tx_end) %>% 
    summarise(index = list(index)) %>% 
    ungroup() %>% 
    mutate(new_idx = row_number())
  if (add_minus) gtf_nested_res <- gtf_nested_res %>% mutate(new_idx = new_idx * -1L)
  list(gtf = gtf_nested_res %>% select(-index),
       tbl = gtf_nested_res %>% select(new_idx, index) %>% unnest())
}

include_nest_mono_any_2 <- function(gtf_nested_mono, gtf_nested_any, index_colname = "index") {
  message("include_nest_mono_any_2")
  if (nrow(gtf_nested_mono) == 0L || nrow(gtf_nested_any) == 0L) {
    return(tibble(query = integer(), subject = integer(), .name_repair = "minimal"))
  } else {
    message("gtf_longer")
    gtf_long_any <- 
      gtf_nested_any %>% 
      split_df(n_worker) %>% 
      future_map_dfr(gtf_longer, .progress = TRUE) %>% 
      get_lastexon(index_colname)
    message("findOverlaps")
    nest_table <- 
      inner_join(
        gtf_nested_mono %>% nest(-c(seqname, strand), .key = "nest_query"),
        gtf_long_any %>% nest(-c(seqname, strand), .key = "nest_subject"),
        by = c("seqname", "strand")
      )
    future_map2_dfr(nest_table$nest_query, nest_table$nest_subject, ~{
      within_matrix <- 
        IRanges::findOverlaps(
          .x %>% as_GRangesList("tx_start", "tx_end", index_colname),
          .y %>% as_GRangesList("start", "end", index_colname),
          type = "within",
          select = "all"
        ) %>% 
        as.matrix() %>% 
        as.data.frame()
      tibble(
        query = .x[[index_colname]][within_matrix$queryHits],
        subject = .y[[index_colname]][within_matrix$subjectHits],
        .name_repair = "minimal"
      )
    }) %>% 
      df_query_subject()
  }
}

as_GRangesList <- function(gtf, start_col, end_col, index_col) {
  GRanges(seqnames = "1", strand = "+", ranges = IRanges(gtf[[start_col]], gtf[[end_col]])) %>% split(gtf[[index_col]])
}

include_nest_mono_1 <- function(gtf_nested_mono, index_colname = "index") {
  message("include_nest_mono_1")
  if (nrow(gtf_nested_mono) == 0L) {
    return(tibble(query = integer(), subject = integer(), .name_repair = "minimal"))
  } else {
    message("gtf_longer")
    gtf_long_any <-
      gtf_nested_mono %>% 
      split_df(n_worker) %>% 
      future_map_dfr(gtf_longer, .progress = TRUE)
    message("findOverlaps")
    nest_table <- 
      inner_join(
        gtf_nested_mono %>% nest(-c(seqname, strand), .key = "nest_query"),
        gtf_long_any %>% nest(-c(seqname, strand), .key = "nest_subject"),
        by = c("seqname", "strand")
      )
    future_map2_dfr(nest_table$nest_query, nest_table$nest_subject, ~{
      within_matrix <- 
        GenomicRanges::findOverlaps(
          .x %>% as_GRangesList("tx_start", "tx_end", index_colname),
          .y %>% as_GRangesList("start", "end", index_colname),
          type = "within",
          select = "all"
        ) %>% 
        as.matrix() %>% 
        as.data.frame(stringsAsFactors = FALSE)
      tibble(
        query = .x[[index_colname]][within_matrix$queryHits],
        subject = .y[[index_colname]][within_matrix$subjectHits],
        .name_repair = "minimal"
      ) %>% 
        filter(subject != query) %>% 
        filter(!subject %in% query)
    }) %>% 
      df_query_subject()
  }
}

certify_index <- function(df) {# In case unnest drop columns with empty lists
  if (!"index" %in% colnames(df)) df %>% mutate(index = integer())
  else df
}

df_dammy_res_inter <- tibble(index = integer(), new_idx = integer(), query = integer(), .name_repair = "minimal")

is_fun <- function(obj_name) is.function(eval(parse(text = obj_name)))

filter_function <- function(names) names[!purrr::map_lgl(names, is_fun)]


# preparation----
message("===preparation===")
message("indexing")
gtf_sample <- gtf_sample %>% mutate(index = group_indices(gtf_sample, sample, transcript_id))
index_tbl <- gtf_sample %>% distinct(index, sample, gene_id, transcript_id)

message("gather exon")
gtf_nested_sample <-
  gtf_sample %>% 
  split_df_group(n_worker, "index") %>%
  future_map_dfr(isocan::gather_exon) %>% 
  isocan::separate_ends()
message("separate multi/mono")
gtf_nested_multi_sample <- gtf_nested_sample %>% dplyr::filter(exon_structure != "-")
gtf_nested_mono_sample <- gtf_nested_sample %>% dplyr::filter(exon_structure == "-")
rm(gtf_sample)
rm(gtf_nested_sample)
# Multi. merge sj-match isoforms----
message("===merge isoforms===")
message("merge sj-match isoforms")
message("  multi")
nest_multi <- collapse_sj_match(gtf_nested_multi_sample)
message("  mono")
nest_mono <- collapse_tx_match(gtf_nested_mono_sample, add_minus = TRUE)
rm(gtf_nested_multi_sample)
rm(gtf_nested_mono_sample)
rm(list = setdiff(ls(all.names = TRUE), ls(all.names = FALSE)))

# Mono - Multi inclusion----
message("Mono - Multi inclusion")
if (mono_multi_merge == TRUE) {
  inc_mono_multi_res <- include_nest_mono_any_2(nest_mono$gtf, nest_multi$gtf, index_colname = "new_idx")
  # print(inc_mono_multi_res)
  include_mono_multi_sample <- 
    inc_mono_multi_res %>% 
    rename(new_idx = "subject") %>% 
    left_join(nest_mono$tbl, by = c("query" = "new_idx")) %>% 
    unnest() %>% certify_index()
  message("  update")
} else {
  include_mono_multi_sample <- df_dammy_res_inter
}
gtf_nest_mono_include <- nest_mono$gtf %>% filter(new_idx %in% include_mono_multi_sample$query)
gtf_nest_mono_noninclude <- nest_mono$gtf %>% filter(!new_idx %in% include_mono_multi_sample$query)

rm(list = setdiff(ls(all.names = TRUE), ls(all.names = FALSE)))

# Mono - Mono inclusion----
message("Mono - Mono inclusion")
if (collapse_mono_range == TRUE) {
  collapse_mono_range_res <- collapse_mono_inter(gtf_nest_mono_noninclude, "new_idx")
  match_mono_sample <- 
    collapse_mono_range_res[[1]] %>% 
    left_join(nest_mono$tbl, by = c("old_idx" = "new_idx")) %>% 
    unnest() %>% certify_index()
  include_mono_sample <- 
    collapse_mono_range_res[[2]] %>% 
    left_join(nest_mono$tbl, by = c("old_idx" = "new_idx")) %>% 
    unnest() %>% certify_index()
  gtf_nest_mono_collapsed <- 
    collapse_mono_range_res[[3]] %>% 
    rename(new_idx = "collapsed_idx")
} else {
  include_mono_sample <- 
    include_nest_mono_1(gtf_nest_mono_noninclude, index_colname = "new_idx") %>% 
    rename(new_idx = "subject") %>% 
    left_join(nest_mono$tbl, by = c("query" = "new_idx")) %>% 
    unnest() %>% certify_index()
  
  gtf_nest_mono_collapsed <- 
    gtf_nest_mono_noninclude %>% filter(!new_idx %in% include_mono_sample$query)
  match_mono_sample <- 
    gtf_nest_mono_collapsed %>% 
    select(new_idx) %>%
    left_join(nest_mono$tbl,
              by = "new_idx") %>%
    unnest() %>% certify_index()
}

rm(list = setdiff(ls(all.names = TRUE), ls(all.names = FALSE)))
#included-Mono - Mono inclusion----
message("included-Mono - Mono inclusion")
include_mono_mono_sample <- 
  include_nest_mono_any_2(gtf_nest_mono_include, gtf_nest_mono_collapsed, index_colname = "new_idx") %>% 
  rename(new_idx = "subject") %>% 
  left_join(nest_mono$tbl, by = c("query" = "new_idx")) %>% 
  unnest() %>% certify_index()

tidy_corr_sample <- function(df, type = "include", info = "") {
  message(paste0(info, ". nrow: ", nrow(df)))
  # print(df)
  if (nrow(df) == 0L) df <- df %>% mutate(index = integer())
  df %>% select(new_idx, index) %>% mutate(type = type)
}
corr_merge <- 
  bind_rows(
    nest_multi$tbl %>% tidy_corr_sample("intronic_match_include", info = "mch_inc_multi"),
    include_mono_multi_sample %>% tidy_corr_sample(info = "mono_multi"),
    include_mono_sample %>% tidy_corr_sample(info = "mono"),
    match_mono_sample %>% tidy_corr_sample("monoex_match", info = "mch_mono"),
    include_mono_mono_sample %>% tidy_corr_sample(info = "mono_mono")
  )

#result----
message("===summarize===")
gtf_collapsed <- 
  bind_rows(nest_multi$gtf, gtf_nest_mono_collapsed) %>% 
  select(new_idx, seqname, strand, tx_start, tx_end, exon_structure) %>% 
  split_df(n_worker) %>% 
  future_map_dfr(gtf_longer)
ex_grl_merge <- 
  gtf_collapsed %>% 
  GRanges() %>% 
  split(.$new_idx)
tx_gr_merge <- range(ex_grl_merge)

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

gene_index_df <- 
  include_2_mod(
    tx_gr_merge,
    tx_gr_merge %>% unlist() %>% GenomicRanges::reduce() %>% split_rows()
  )
gene_index_df <- 
  gene_index_df %>% 
  mutate(gene_id = group_indices(gene_index_df, super_index) %>% paste0("PB.", .)) %>% 
  select(index, gene_id)
gene_index_df <- 
  gene_index_df %>% 
  group_by(gene_id) %>% 
  mutate(transcript_id = row_number()) %>% 
  ungroup() %>% 
  mutate(transcript_id = paste0(gene_id, ".", transcript_id)) %>% 
  rename(new_idx = "index")

corr_result_sample <- 
  full_join(corr_merge, gene_index_df, by = "new_idx") %>% 
  full_join(index_tbl %>% 
              rename(tx_id_sample = "transcript_id", g_id_sample = "gene_id"),
            by = "index") %>% 
  select(type, gene_id, transcript_id, sample, tx_id_sample, g_id_sample)

classify <- function(intra_type, inter_type) {
  allowed_types = c("intronic_match_include", "monoex_match", "include")
  stopifnot(all(intra_type %in% allowed_types))
  stopifnot(all(inter_type %in% allowed_types))
  if_else(intra_type == inter_type, intra_type, "include")
}

corr_final <- 
  full_join(
    corr_result_sample %>% 
      rename(inter_type = "type"),
    corr_sample %>% 
      rename(intra_type = "type", g_id_sample = "gene_id", tx_id_sample = "transcript_id"),
    by = c("sample", "g_id_sample", "tx_id_sample")
  ) %>% 
  mutate(type = classify(intra_type, inter_type)) %>% 
  select(gene_id, transcript_id, sample, g_id_sample, tx_id_sample, seq_id, type, intra_type, inter_type)

gtf_final <- 
  ex_grl_merge %>% 
  as_tibble() %>% 
  mutate(new_idx = as.integer(group_name)) %>% 
  rename(seqname = "seqnames") %>% 
  select(new_idx, seqname, start, end, strand) %>% 
  full_join(gene_index_df, by = "new_idx") %>% 
  select(-new_idx) %>% 
  mutate(feature = "exon")

#count full read----
#only unique-match CCS
unique_corr_final <- 
  corr_final %>% 
  anti_join(corr_final %>% select(seq_id, sample) %>% mutate(dupli = duplicated(.)) %>% filter(dupli),
            by = c("seq_id", "sample"))
count_final <- 
  unique_corr_final %>% 
  add_pbcount_map() %>% 
  group_by(transcript_id, sample) %>%
  summarise(count_fl = sum(FL)) %>%
  ungroup() %>%
  rename(pbid = "transcript_id") %>%
  mutate(
    count_nfl = 0,
    count_nfl_amb = 0,
    norm_fl = 0,
    norm_nfl = 0,
    norm_nfl_amb = 0)
count_final_nest <- 
  count_final %>%
  nest(-sample) %>% 
  full_join(input_files[c("sample", count_out_colname)], by = "sample")

#write out----
write_tsv(corr_final, corr_out)
isocan::write_gtf(gtf_final, gtf_out)
map2(count_final_nest[["data"]], count_final_nest[[count_out_colname]], write_tsv)
