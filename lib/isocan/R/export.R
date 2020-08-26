#' export isodata as SAM format.
#'
#' You MUST confirm that all variants lie ONLY on exons.
#' Some arguments are set to fixed values:
#' MAPQ = 255, RNEXT = asterisk, PNEXT = 0, TLEN = 0, QUAL = asterisk.
#' @param isodata An isodata class object.
#' @param path A file path to create.
#' @param sample_name A sample name to export.
#' @param ref_genome A BSgenome class object. It is used to construct SEQ segment and SQ tags.
#' @export
export_sam <- function(isodata, path, sample_name, ref_genome){
  #first, should check the format. e.g. do not contain "*" in strand
  #Must do filter_on_exon before execute this function.

  isodata@variant %>%
    filter(sample == sample_name) %>%
    .$chr %>%
    unique() -> chr_variant
  isodata@exon_ranges$seqnames %>% unique() %>% as.character() -> chr_exon_ranges
  dplyr::intersect(chr_variant, chr_exon_ranges) -> chr_iso
  dplyr::intersect(chr_iso, names(ref_genome)) -> chr_common
  dplyr::setdiff(chr_iso, chr_common) -> chr_drop

    if(length(chr_drop) > 0) warning(str_c(
    length(chr_drop),
    " chromosomes do NOT exist in ref_genome. Omitted:\n",
    str_c(chr_drop, collapse = ", ")))

  SQ <- get_sq(ref_genome, chr_common)
  PG <- "@PG\tID:isocan\tPN:isocan"
  BODY <- get_sam(isodata, ref_genome, sample_name, chr_common)

  readr::write_tsv(SQ, path, col_names = F)
  cat(PG, "\n", file = path, append = T)
  readr::write_tsv(BODY, path, append = T, col_names = F)
}

get_sq <- function(ref_genome, chr_common){
  tibble(SN = chr_common) %>%
    mutate(SQ = "@SQ",
           LN = purrr::map_chr(SN, ~ ref_genome[[.x]] %>% length())
           ) %>%
    mutate(SN = str_c("SN:", SN),
           LN = str_c("LN:", LN)) %>%
    select(SQ, SN, LN)
}

get_sam <- function(isodata, ref_genome, sample_name, chr_common){
  isodata@variant %>%
    filter(sample == sample_name &
             chr %in% chr_common) %>%
    select(seq_id, chr, type, start, end, alt) -> seq_df #do not check "ref" column is consistent to ref_genome


  iso@tx_seq_correspondence %>%
    select(tx_id, seq_id) %>%
    left_join(seq_df, by = "seq_id", na_matches = "never") %>% #export all seqences listed in tx_seq_correspondence
    tidyr::nest(-c(seq_id, tx_id), .key = "var") -> seq_df


  iso@exon_ranges %>%
    select(seqnames, start, end, strand, transcript_id) %>%
    rename(tx_id = transcript_id) %>%
    filter(tx_id %in% seq_df$tx_id) %>%
    arrange(tx_id, start) -> exon_df

  exon_df %>%
    group_by(tx_id) %>%
    summarize(
      FLAG = if_else(strand[1] == "+", 0, 16),
      RNAME = seqnames[1] %>% as.character(),
      POS = start[1],
      seq_pos = alternate_colon(start, end) %>% list()
    ) %>%
    mutate(seq = purrr::map2(RNAME, seq_pos, ~ ref_genome[[.x]][.y] %>% as.character()))-> exon_df_sum

  exon_df %>%
    tidyr::nest(-tx_id, .key = "exon") %>%
    left_join(exon_df_sum, by = "tx_id") -> exon_df

  inner_join(seq_df, exon_df,
               by = "tx_id",
               na_matches = "never") -> seq_df
  #seq_df has 8 columns: seq_id, tx_id, var, exon, FLAG, RNAME, POS, seq

  seq_df %>%
    mutate(cigar_seq = purrr::pmap(list(var, exon, seq, POS), construct_CIGAR_SEQ, ref_genome)) %>%
    select(seq_id, cigar_seq, FLAG, RNAME, POS) %>%
    rename(QNAME = seq_id) %>%
    tidyr::unnest() %>%
    mutate(MAPQ = 255,
           RNEXT = "*",
           PNEXT = 0,
           TLEN = 0,
           QUAL = "*") %>%
    select(QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL)
}

construct_CIGAR_SEQ <- function(var, exon, seq, POS, ref_genome){
  #return CIGAR, SEQ
  var %>% mutate(start_relative = start - POS + 1,
                 end_relative = end - POS + 1) -> var
  exon %>% mutate(start_relative = start - POS + 1,
                  end_relative = end - POS + 1) -> exon

  cigar = construct_CIGAR(var, exon)
  seq = construct_SEQ(var, exon, seq, ref_genome)

  tibble(CIGAR = cigar,
         SEQ = seq)
}

construct_CIGAR <- function(var, exon){
  exon_width <- exon$end_relative - exon$start_relative + 1
  intron_width <- exon$start_relative - lag(exon$end_relative) - 1
  exon_intron_width <- alternate_c(intron_width, exon_width)[-1]
  cigar <- S4Vectors::Rle(
    rep(c("N", "M"), length(exon_width))[-1],
    exon_intron_width)

  #del
  var %>%
    filter(type == "-") -> var_del
  alternate_colon(var_del$start_relative + 1, # 0-based -> 1-based
                  var_del$end_relative) -> range_del
  cigar[range_del] <- "D"

  #ins : must be the last
  var %>%
    filter(type == "+") %>%
    mutate(ins_len = str_length(alt)) -> var_ins
    purrr::reduce2(var_ins$start_relative, var_ins$ins_len,
                   function(a, x, y) c(a[1:x], rep("I", y), a[(x + 1) : length(a)]),
                   .init = cigar
                   ) -> cigar

    str_c(S4Vectors::runLength(cigar), S4Vectors::runValue(cigar), collapse = "")
}

alternate_c <- function(x, y){
  purrr::map2(x, y, ~ c(.x, .y)) %>% purrr::flatten_dbl()
}

alternate_colon <- function(x, y){
  purrr::map2(x, y, ~ .x:.y) %>% purrr::flatten_dbl()
}

construct_SEQ <- function(var, exon, seq, ref_genome){
  relative_exon_vector <- tibble(pos_relative = alternate_colon(exon$start_relative, exon$end_relative))
  relative_exon_vector %>%
    mutate(index_relative = 1:length(relative_exon_vector)) -> relative_exon_vector
  var %>%
    left_join(relative_exon_vector,
              by = c("start_relative" = "pos_relative")) %>%
    rename(start_index = index_relative) %>%
    left_join(relative_exon_vector,
              by = c("end_relative" = "pos_relative")) %>%
    rename(end_index = index_relative) -> var

  var %>% filter(type == "*") -> var_mut
  var %>% filter(type == "-") -> var_del
  var %>% filter(type == "+") -> var_ins

  #mut
  var_mut %>% mutate(start_end_index = purrr::map2(start_index, end_index, ~ c(.x, .y))
                     ) -> var_mut
  purrr::reduce2(var_mut$start_end_index, var_mut$alt,
                 function(a, x, y) `str_sub<-`(a, x[1], x[2], value = y),
                 .init = seq) -> seq

  #del
  purrr::reduce2(var_mut$start_index, var_mut$end_index,
                 function(a, x, y) `str_sub<-`(a, x, y, value = paste0(rep("_", y - x + 1))),
                 .init = seq) -> seq #"_" is a spaceholder.

  #ins
  var_ins %>% arrange(desc(start_index)) -> var_ins # not to change index on inserting.
  purrr::reduce2(var_ins$start_index, var_ins$alt,
                 function(a, x, y) str_c(str_sub(a, end = x), y, str_sub(a, start = x+1)),
                 .init = seq) -> seq

  #delete spaceholders
  seq %>% str_replace_all("_", "")
}
