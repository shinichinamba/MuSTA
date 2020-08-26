#' Filtering 'variant' slot using 'exon_ranges' slot and provided true variant tibble.
#'
#' @param isodata An isodata class object with 'variant' and 'exon_ranges' slots.
#' @param df A tibble containing variant information to be retained. It is usually generated from WES or WGS data.
#' If unspecified, filtering is conducted using only 'exon_ranges' slot.
#' @param ref_genome A BSgenome class object. if specified, spliced introns that exist on exons in exon_range slot will be treated as deletion.
#' If unspecified, spliced introns will be ignored.
#' @param .restrict_on_exon A logical. Whether filtering variants to retain only those on exons of correlated transcripts.
#' @export
filter_variant_isodata <- function(isodata, df = NULL, ref_genome = NULL, .restrict_on_exon = T){
  if(is_isodata(isodata) == F) stop("The first argument must be isodata class.")
  if(!is.null(df) && !"variant_without_seq_id" %in% distinguish_slot(df)) stop("The second argument must satisfy isodata@variant slot requirement. ('seq_id' column is unnecessary.)")
  is_null_ref <- is.null(ref_genome)

  isodata <- if_else(.restrict_on_exon, filter_on_exon(isodata, is_null_ref), isodata@variant)
  if(!is_null_ref) variant <- fill_na_w_ref(isodata@variant, ref_genome)
  if(!is.null(df)) variant <- filter_variant_w_df(variant, df)

  isodata@variant <- variant
  return(isodata)
}

#' Filtering 'variant' slot using 'exon_ranges' slot and provided true variant tibble.
#'
#' @param variant A tibble
#' @param true_variant A tibble containing variant information to be retained. It is usually generated from WES or WGS data.
#' @param ref_genome A BSgenome class object.
#' @export
filter_variant <- function(variant, true_variant, ref_genome){
  if (!"variant" %in% distinguish_slot(variant)) stop("variant requirement failure")
  if (!"variant_without_seq_id" %in% distinguish_slot(true_variant)) stop("true_variant requirement failure")
  if (nrow(true_variant) == 0L) {
    tibble::tribble(
      ~chr, ~type, ~start, ~end, ~ref, ~alt, ~sample, ~status,
      character(1), character(1), integer(1), integer(1), character(1), character(1), character(1), character(1)
    ) %>% dplyr::slice(-1)
  } else {
    variant %>%
      add_another_start_repeptive(ref_genome) %>%
      add_another_start_sj() %>%
      filter_variant_w_df(true_variant)
  }
}

# Fit the slot 'variant' for the slot 'exon_ranges'.
# splice sites are handled as the sub class of deletions, and variant outside of exon_ranges will be filtered.
filter_on_exon <- function(isodata, is_null_ref, ref_genome){
  if (is_isodata(isodata) == F) stop("The first argument must be isodata class.")
  isodata@variant %>%
    add_another_start_repeptive(ref_genome) %>%
    add_another_start_sj() %>% #must be ahead of masking splice sites
    convert_splice_del(is_null_ref) %>%
    filter_outside_variant(isodata@exon_ranges, isodata@tx_seq_correspondence) -> isodata@variant
  return(isodata)
}

add_another_start_sj <- function(df) {
  df %>% arrange(seq_id, sample, chr, start, end) -> df
  if (!"another_start" %in% colnames(df)) mutate(df, another_start = NA_integer_) -> df

  df %>%
    select(seq_id, start, end, sample, chr, type) %>%
    mutate(lag_seq_id = lag(seq_id),
           lag_start = lag(start),
           lag_end = lag(end),
           lag_sample = lag(sample),
           lag_chr = lag(chr),
           lag_type = lag(type)) -> tmp

  (tmp$lag_seq_id == tmp$seq_id &
      tmp$lag_sample == tmp$sample &
      tmp$lag_chr == tmp$chr &
      tmp$lag_start == tmp$start &
      tmp$lag_type == "+" &
      tmp$type == "~") -> is_splice_w_last_ins

  (tmp$lag_seq_id == tmp$seq_id &
      tmp$lag_sample == tmp$sample &
      tmp$lag_chr == tmp$chr &
      tmp$lag_end == tmp$end &
      tmp$lag_type == "~" &
      tmp$type == "+") -> is_ins_w_last_splice

  df[tidyr::replace_na(lead(is_splice_w_last_ins), F), "another_start"] <- df[tidyr::replace_na(is_splice_w_last_ins, F), "end"]
  df[tidyr::replace_na(is_ins_w_last_splice, F), "another_start"] <- df[tidyr::replace_na(lead(is_ins_w_last_splice), F), "start"]
  return(df)
}

get_repeat_unit <- function(chr_vec) stringr::str_replace(chr_vec, "^(.+)\\1+$", "\\1")

add_another_start_repeptive <- function(df, ref_genome) {
  if (!is(ref_genome, "BSgenome")) stop("The second argument must be a BSgenome object.")
  if (!"another_start" %in% colnames(df)) mutate(df, another_start = NA_integer_) -> df
  df %>%
    mutate(.index = row_number()) %>%
    select(seq_id, sample, chr, tx_start, tx_end, start, end, type, ref, alt, .index) %>%
    arrange(seq_id, sample, chr, tx_start, tx_end, start) %>%
    mutate(last_end_gap = lag(end) - start,
           next_start_gap = lead(start) - end) -> df_idx

  fun_equal <- function(colname, df, fun) df[[colname]] == fun(df[[colname]])
  na_to_false <- function(vec) `[<-`(vec, is.na(vec), FALSE)

  is_last_same_seq <- purrr::map(c("seq_id", "sample", "chr", "tx_start", "tx_end"), fun_equal, df_idx, lag) %>% reduce(`&`) %>% na_to_false()
  is_next_same_seq <- purrr::map(c("seq_id", "sample", "chr", "tx_start", "tx_end"), fun_equal, df_idx, lead) %>% reduce(`&`) %>% na_to_false()
  df_idx[!is_last_same_seq, "last_end_gap"] <- c(df_idx[!is_last_same_seq, "tx_start"] - df_idx[!is_last_same_seq, "start"])
  df_idx[!is_next_same_seq, "next_start_gap"] <- c(df_idx[!is_last_same_seq, "tx_end"] - df_idx[!is_last_same_seq, "end"])

  bind_rows(
    df_idx %>%
      filter(type == "+") %>%
      mutate(search_seq = get_repeat_unit(alt)),
    df_idx %>%
      filter(type == "-") %>%
      mutate(search_seq = get_repeat_unit(ref))
  ) %>%
    mutate(len_search = str_length(search_seq)) -> df_idx
  df_idx %>%
    mutate(
      last_ref = retrieve_ref_seq_2(chr, "+", start - len_search + 1L, start, ref_genome),
      next_ref = retrieve_ref_seq_2(chr, "+", end + 1, end + len_search, ref_genome),
      last_match = (last_ref == search_seq),
      next_match = (next_ref == search_seq)
    ) -> df_idx
  if (any(df_idx$len_search > df_idx$last_end_gap)) df_idx[df_idx$len_search > df_idx$last_end_gap, "last_match"] <- FALSE
  if (any(df_idx$len_search > df_idx$next_start_gap)) df_idx[df_idx$len_search > df_idx$next_start_gap, "next_match"] <- FALSE
  df_idx %>%
    filter(last_match | next_match) -> df_idx
  if (any(df_idx$last_match & df_idx$next_match)) {
    print(as.data.frame(df_idx[df_idx$last_match & df_idx$next_match, ]))
    stop("Bad alignment tool for repeptive region.")
    }

  get_all_repeat_len_to3prime <- function(chr, position, rep_pattern, ref_genome, .detected_n_rep = 0L) {
    rep_len <- str_length(rep_pattern)
    if (retrieve_ref_seq(chr, position + 1L, position + rep_len, ref_genome) == rep_pattern) {
      .detected_n_rep <- get_all_repeat_len_to3prime(chr, position + rep_len, rep_pattern, ref_genome, .detected_n_rep + 1L)
    }
    return(.detected_n_rep)
  }
  get_all_repeat_len_to5prime <- function(chr, position, rep_pattern, ref_genome, .detected_n_rep = 0L) {
    rep_len <- str_length(rep_pattern)
    if (retrieve_ref_seq(chr, position - rep_len + 1L, position, ref_genome) == rep_pattern) {
      .detected_n_rep <- get_all_repeat_len_to5prime(chr, position - rep_len, rep_pattern, ref_genome, .detected_n_rep + 1L)
    }
    return(.detected_n_rep)
  }

  bind_rows(
    df_idx %>%
      filter(last_match) %>%
      mutate(n_rep = purrr::pmap_int(list(chr, start, search_seq), get_all_repeat_len_to5prime, ref_genome),
             another_start = start - str_length(search_seq) * n_rep) %>%
      select(.index, another_start),
    df_idx %>%
      filter(next_match) %>%
      mutate(n_rep = purrr::pmap_int(list(chr, end, search_seq), get_all_repeat_len_to3prime, ref_genome),
             another_start = start + str_length(search_seq) * n_rep) %>%
      select(.index, another_start)
  ) -> another_pos_df
  df[another_pos_df$.index, "another_start"] <- another_pos_df$another_start
  df
}

convert_splice_del <- function(df, .drop = F){
  if (.drop) {
    df %>% filter(type != "~") -> df
  }else{
  df[df$type == "~", "ref"] <- NA
  df[df$type == "~", "type"] <- "-"
  }
  return(df)
}

filter_outside_variant <- function(variant, exon_ranges, tx_seq_correspondence){
  variant %>%
    inner_join(tx_seq_correspondence %>% select(seq_id, tx_id) %>% rename(.tx_id = tx_id),
               by = "seq_id",
               na_matches = "never") %>%
    inner_join(exon_ranges %>%
                 mutate(seqnames = as.character(seqnames)) %>%
                 select(transcript_id, seqnames, start, end) %>%
                 rename(.tx_id = transcript_id,
                        chr = seqnames,
                        .start = start,
                        .end = end),
               by = c(".tx_id", "chr"),
               na_matches = "never") -> variant
  variant %>%
    filter(.start <= start & start <= .end & .start <= end & end <= .end) -> tmpOK
  variant %>% anti_join(tmpOK) -> variant
  tmpOK[tidyr::replace_na(tmpOK$another_start < tmpOK$.start | tmpOK$.end < tmpOK$another_start, F),
        "another_start"] <- NA

  variant %>%
    filter(.start <= another_start & another_start <= .end) -> tmp
  variant %>% anti_join(tmp) -> variant
  tmp[,"start"] <- tmp$another_start
  tmp[,"end"] <- tmp$another_start
  tmp[,"another_start"] <- rep(NA, nrow(tmp))

  variant %>% #rescue del.
    filter(type == "-") -> variant
  variant %>%
    mutate(.max_start = pmax(start, .start),
           .min_end = pmax(end, .end)) %>%
    filter(.min_end > .max_start) %>%
    mutate(.relative_start = .max_start - start + 1,
           .relative_end = .min_end - start + 1,
           ref = str_sub(ref, .relative_start, .relative_end)) -> tmp2

  bind_rows(tmpOK %>% select(-starts_with(".")),
            tmp   %>% select(-starts_with(".")),
            tmp2  %>% select(-starts_with(".")))
}

filter_variant_w_df <- function(variant, df, .allow_inclusion = TRUE){
  variant %>% filter(type != "~") %>% mutate(.index = row_number()) -> variant
  rlang::inform("===filter variants===")
  varOK <- semi_join(variant, df, by = c("chr", "type", "start", "end", "ref", "alt", "sample"))

  is_another_pos = ("another_start" %in% colnames(variant))

  #partial match -> ambiguous status. matchしなかった部分を取り出
  if (is_another_pos) {
    rlang::inform("column 'another_start' is detected.")
    change_inssite_another <- function(df){
      df %>%
        filter(!is.na(another_start)) %>%
        mutate(start = another_start,
               end = another_start)
    }
    change_delsite_another <- function(df){
      df %>%
        filter(!is.na(another_start)) %>%
        mutate(end = end - start + another_start,
               start = another_start)
    }

    variant %>%
      filter(type == "+") %>%
      anti_join(varOK) -> var_left_ins
    df %>%
      filter(type == "+") -> df_ins

    var_left_ins %>%
      change_inssite_another() %>%
      semi_join(df_ins, by = c("chr", "start", "end", "ref", "alt", "sample")
                ) -> varOK_alt_ins

    variant %>%
      filter(type == "-") %>%
      anti_join(varOK) -> var_left_del
    df %>%
      filter(type == "-") -> df_del

    var_left_del %>%
      change_delsite_another() %>%
      semi_join(df_del, by = c("chr", "start", "end", "ref", "alt", "sample")
      ) -> varOK_alt_del

    varOK <- bind_rows(varOK, varOK_alt_ins, varOK_alt_del)
  }
  variant %>%
    anti_join(varOK, by = ".index") -> var_left
  varOK <- varOK %>% mutate(status = "match")

  if (.allow_inclusion) {
    rlang::inform("check included variants\ninsertion")
    #ins
    var_left %>%
      filter(type == "+") -> var_left_ins
    df %>%
      filter(type == "+") -> df_ins
    df_ins %>%
      rename(.alt_df = alt) %>%
      select(chr, type, start, end, .alt_df, sample) -> df_ins

    find_inclusion_ins <- function(var_left_ins, df_ins){
      var_left_ins %>%
        inner_join(df_ins,
                   by = c("chr", "type", "start", "end", "sample"),
                   na_matches = "never") %>%
        filter(str_detect(alt, .alt_df)) %>%
        mutate(alt = .alt_df) %>%
        select(-.alt_df)
    }

    find_inclusion_ins(var_left_ins, df_ins) -> varOK_inclusion_ins

    if (is_another_pos) {
      rlang::inform("  another_position")
      var_left_ins %>%
        anti_join(varOK_inclusion_ins, by = ".index") -> var_left_ins

      var_left_ins %>%
        change_inssite_another() %>%
        find_inclusion_ins(df_ins) %>%
        bind_rows(varOK_inclusion_ins) -> varOK_inclusion_ins
    }

    #del
    rlang::inform("deletion")
    var_left %>%
      filter(type == "-") %>%
      anti_join(varOK) -> var_left_del
    df %>%
      filter(type == "-") -> df_del
    df_del %>%
      rename(.start_df = start,
             .end_df = end,
             .ref_df = ref) %>%
      select(chr, type, .start_df, .end_df, .ref_df, sample) -> df_del

    rlang::inform("  nest and join")
    inner_join(# join wo/ nest will require TOO MANY MEMORY.
      var_left_del %>%
        tidyr::nest(-c(chr, type, alt, sample), .key = "var"),
      df_del %>%
        tidyr::nest(-c(chr, type, sample), .key = "df"),
      by = c("chr", "type", "sample"),
      na_matches = "never"
    ) -> nested_var_df

    find_inclusion_del <- function(var_nested, df_nested, info = NULL){
      if (!is.null(info)) rlang::inform(info)
      inner_join(
        var_nested %>% mutate(.dammy = 0L),#全組み合わせでjoin
        df_nested %>% mutate(.dammy = 0L),
        by = ".dammy"
      ) %>%
        filter(start <= .start_df & .end_df <= end) %>%
        mutate(.ref_var = str_sub(ref, .start_df - start + 1L, .end_df - start)) %>%
        filter(!is.na(.ref_var) & !is.na(.ref_df)) %>%
        filter(.ref_var == .ref_df) -> inclusion_del

      inclusion_del %>%
        mutate(start = .start_df,
               end = .end_df,
               ref = .ref_df) %>%
        select(-c(".dammy", ".ref_var", ".ref_df" ,".start_df", ".end_df"))
    }

    nested_var_df %>%
      mutate(
        .index_2 = row_number(),
        res = purrr::pmap(list(var, df, .index_2), find_inclusion_del)) %>%
      select(-c(df, var, .index_2)) %>%
      tidyr::unnest() -> varOK_inclusion_del

    if (is_another_pos) {
      rlang::inform("  another_position")
      if (nrow(varOK_inclusion_del) > 0L) var_left_del %>% anti_join(varOK_inclusion_del, by = ".index") -> var_left_del

     inner_join(# join wo/ nest will require TOO MANY MEMORY.
       var_left_del %>%
         change_delsite_another() %>%
          tidyr::nest(-c(chr, type, alt, sample), .key = "var"),
        df_del %>%
          tidyr::nest(-c(chr, type, sample), .key = "df"),
        by = c("chr", "type", "sample"),
        na_matches = "never"
      ) %>%
       mutate(
         .index_2 = row_number(),
         res = purrr::pmap(list(var, df, .index_2), find_inclusion_del)) %>%
       select(-c(df, var, .index_2)) %>%
       tidyr::unnest() %>%
        bind_rows(varOK_inclusion_del) -> varOK_inclusion_del
    }

    #summarise
    rlang::inform("summarize")
    bind_rows(varOK_inclusion_ins, varOK_inclusion_del) %>%
      mutate(status = "ambiguous") -> varOK_inclusion
    bind_rows(
      varOK,
      varOK_inclusion
    ) -> varOK
  }

  if (is_another_pos) varOK %>% select(-another_start) -> varOK
  varOK %>%
    arrange(.index) %>%
    select(-.index) -> varOK
  rlang::inform("===done===")
  return(varOK)
}

flag_variant_w_df <- function(variant, df, .allow_inclusion = T){
  #
}


# TOO SLOW! Should not map. Change to retrieve once par chr.
fill_na_w_ref <- function(df, ref_genome){
  if (!is(ref_genome, "BSgenome")) stop("The second argument must be a BSgenome object.")
  #when df$ref is NA, fill it using ref_genome
  #BSgenome objects are 1-based.
  df %>% filter(is.na(ref)) -> na_df
  anti_join(df, na_df) -> rest_df
  na_df %>%
    mutate(
      ref = retrieve_ref_seq_2(chr, "+", start + 1L, end, ref_genome)
      ) -> na_df

  bind_rows(rest_df, na_df)
}

retrieve_ref_seq <- function(seqname, start, end, ref_genome){
  # if (!str_detect(seqname, "^chr")) seqname <- paste0("chr", seqname)
  as.character(ref_genome[[seqname]][start:end]) %>% str_to_lower()
}

retrieve_ref_seq_2 <- function(seqname, strand, start, end, ref_genome) {
  if (length(strand) == 1L) strand <- rep(strand, length(seqname))
  gtf_gr <-
    GenomicRanges::makeGRangesListFromFeatureFragments(
      seqnames = seqname,
      fragmentStarts = purrr::map(start, force),
      fragmentEnds = purrr::map(end, force),
      strand = strand
    )
  as.character(GenomicFeatures::extractTranscriptSeqs(ref_genome, gtf_gr))
}
