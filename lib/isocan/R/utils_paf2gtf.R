#' retreive gtf_oriented tibble from paf
#'
#' @description Retreive granges object from paf data and paf-based variant data
#' @param paf A tbl. A paf file generated by minimap2 with --cs argument and read by read_paf.
#' @param sample_name A one-length character vector. Used to fill the "sample" column.
#' @param withcare Logical. If true, stop when a transcript's end position is NOT consistent between described one and calculated one.
#' @param verbose Logical.
#' @param parallel Logical. If true, use furrr::future_map instead of purrr::map
#' @param paf_variants A tibble generated by convert_paf_variant function. If specified, common procedure with convert_paf_variant will be skipped.
#' @param retreive_type One of 'exon', 'intron', 'tx'(='transcript').
#' @return A gtf_oriented tibble object. 1-based coordination.
#' @export
retrieve_gtf <- function(paf, sample_name, retreive_type = "exon", withcare = F, verbose = T, parallel = T, paf_variants = NULL){
  if(length(sample_name) != 1) stop("The length of sample_name must be 1.")
  if(!retreive_type %in% c("tx", "transcript", "exon", "intron")) stop("'type' must be one of 'exon', 'intron', 'tx'(='transcript').")

  if(retreive_type %in% c("tx", "transcript")){
    paf %>%
      mutate(
        sample = sample_name,
        type = "transcript",
        start = start + 1) %>% #1-based
      arrange(seq_id, strand, chr, start, end) %>% return()
  }else{
    #retrieve introns
    if(is.null(paf_variants)){
      #convert_paf_variant like function
      if (requireNamespace("furrr")) {
        map_choice <- furrr::future_map
        if (verbose) map_choice <- purrr::partial(furrr::future_map, .progress = TRUE)
      } else {
        map_choice <- purrr::map
      }

      paf %>%
        mutate(variant = map_choice(paf$cs, convert_cs_variant)) %>%
        tidyr::unnest() %>%
        mutate(var_start = var_start + start,
               var_end = var_end + start) -> paf2

      if(withcare){
        paf2 %>%
          filter(type == "") %>%
          mutate(success = (end == var_end)) %>%
          filter(success == F) -> fail
        if(nrow(fail) > 0) stop(str_c("Some transcript(s) was failed to convert: ", str_c(fail$seq_id, collapse = ", ")))
      }

      paf2 %>%
        filter(type == "~") %>%
        select("seq_id", "chr", "var_start", "var_end", "start", "end", "strand") -> introns
    }else{
      paf_variants %>%
        filter(sample == sample_name & type == "~") %>%
        rename("var_start" = "start",
               "var_end" = "end",
               "start" = "tx_start",
               "end" = "tx_end",
               "strand" = "tx_strand") %>%
        select("seq_id", "chr", "var_start", "var_end", "start", "end", "strand") -> introns
    }
    #filter isoforms with "intron-insertion-intron" structure (They might be generated for some technological reasons, including mismapping or errornous sequencing)
    introns %>%
      arrange(seq_id, strand, chr, start, end, var_start) %>%
      filter(seq_id == lag(seq_id) &
               strand == lag(strand) &
               chr == lag(chr) &
               start == lag(start) &
               end == lag(end) &
               var_start == lag(var_end)) %>%
      rename(position = "var_start") %>%
      select(seq_id, strand, chr, start, end, position) -> int_ins_int
    if(nrow(int_ins_int) > 0L) {
      warning('"intron-insertion-intron" structure was detected. Those isoforms with this structure were filtered.')
      print('"intron-insertion-intron" structure was detected. Those isoforms with this structure were filtered.')
      print(as.data.frame(int_ins_int))
      introns %>%
        anti_join(int_ins_int, by = c("seq_id", "strand", "chr", "start", "end")) -> introns
    }

    if(retreive_type == "intron"){
      introns %>%
        mutate(
          sample = sample_name,
          type = "intron",
          start = var_start + 1, #1-based
          end = var_end) %>%
        select(seq_id, strand, chr, start, end, sample, type) %>%
        arrange(seq_id, strand, chr, start, end) %>% return()
      }else{
      #intron2exon
      introns %>%
        arrange(seq_id, strand, chr, start, end, var_start) %>%
        distinct(seq_id, strand, chr, start, end, .keep_all = T) %>%
        mutate(start = start + 1, #1-based
               end = var_start) %>%
        select(seq_id, strand, chr, start, end) -> first_exons

      introns %>%
        arrange(seq_id, strand, chr, start, end, desc(var_start)) %>%
        distinct(seq_id, strand, chr, start, end, .keep_all = T) %>%
        mutate(start = var_end + 1) %>% #1-based
        select(seq_id, strand, chr, start, end) -> last_exons

      introns %>%
        arrange(seq_id, strand, chr, start, end, var_start) %>%
        mutate(exon_start = lag(var_end) + 1) %>%
        filter(seq_id == lag(seq_id) &
                 strand == lag(strand) &
                 chr == lag(chr) &
                 start == lag(start) &
                 end == lag(end)) %>%
        select(-c(end, start)) %>%
        rename("end" = "var_start",
               "start" = "exon_start") %>%
        select(seq_id, strand, chr, start, end) -> inner_exons

      paf %>%
        anti_join(introns, by = c("seq_id", "strand", "chr", "start", "end")) %>%
        select(!!! syms(c("seq_id", "strand", "chr", "start", "end"))) %>%
        mutate(
          start = start + 1) -> one_exon_txs #1-based

      bind_rows(first_exons, last_exons, inner_exons, one_exon_txs) %>%
        mutate(sample = sample_name,
               type = "exon") %>%
        arrange(seq_id, strand, chr, start, end)
    }
  }
}

#' convert gtf_tbl to granges object
#'
#' @param gtf_tbl A gtf_oriented tibble object.
#' @return A granges object.
#' @export
convert_gtf_granges <- function(gtf_tbl){
  min_granges_colnames <- c("seq_id", "strand", "chr", "start", "end")
  if (!is(gtf_tbl, "tibble") || all(min_gtf_colnames %in% colnames(gtf_tbl))) {
    stop("Input must be a tibble with columns named seq_id, strand, chr, start, and end.")
  }
  gtf_tbl %>%
    with(GenomicRanges::GRanges(chr, IRanges(start, end), strand, type, gene_id = NA_character_, transcript_id = seq_id))
}
