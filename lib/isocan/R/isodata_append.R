#tmp <- rtracklayer::import("/home/sinnhazime/Desktop/projects/long_read/output/merged/all_samples.chained.5merged.gff")

#' import a GRanges object into an isodata class object.
#'
#' @inheritParams iso_append
isodata_append.GRanges <- function(obj, isodata, append){
  suppressWarnings(obj <- tibble::as_tibble(obj))
  # check chimera
  chimera <- find_chimera(obj)
  if(nrow(chimera) > 0) stop(str_c("chimeric transcripts are not supported in the current version: ", str_c(chimera$transcript_id, collapse = ", ")))
  # import exonic data
  obj %>%
    filter_('type == "exon"') %>%
    select(-"type") -> obj
  if(append == F){
    isodata@exon_ranges <- obj %>% select(-"gene_id")
    isodata@tx_attr <- obj %>% distinct_("transcript_id", "gene_id")
  }else {
    dupli_exon <- intersect(isodata@exon_ranges$transcript_id, obj$transcript_id)
    dupli_tx <- intersect(isodata@tx_attr$transcript_id, obj$transcript_id)
    if(length(dupli_exon) > 0) stop(str_c("transcript_ids are duplicated between newly-added data and old data in the exon_ranges level (and parhaps also in the tx_attr level): ", str_c(dupli_exon, collapse = ", ")))
    if(length(dupli_tx) > 0) stop(str_c("transcript_ids are duplicated between newly-added data and old data in the tx_attr level: ", str_c(dupli_tx, collapse = ", ")))
    isodata@exon_ranges <- bind_rows(isodata@exon_ranges, obj %>% select(-"gene_id"))
    isodata@tx_attr <- bind_rows(isodata@tx_attr, obj %>% distinct_("transcript_id", "gene_id"))
  }
  return(isodata)
}

find_chimera <- function(df){
  # df: tbl from GRanges obj.
  # return: 1-column tbl with the column "transcript_id".
  df %>%
    filter_('type == "exon"') %>%
    mutate_at(vars("transcript_id", "seqnames", "strand"), as.character) %>%
    count(transcript_id, seqnames, strand) %>%
    select("transcript_id") -> tmp
  anti_join(tmp, tmp %>% distinct_("transcript_id"), by = "transcript_id")
}

construct_tx <- function(df){
  # construct tx level data from exon level data.
  # df: exonic tbl from GRanges class.
}


#' import variant file into an isodata class object.
#'
#' @inheritParams iso_append
isodata_append.variant <- function(obj, isodata, append){
  if(append == F){
    isodata@variant <- obj
  }else {
    dupli <- inner_join(isodata@variant, obj)
    dupli <- dupli %>% mutate(tx_sample = str_c(seq_id, sample))
    if(nrow(dupli) > 0) message(str_c("Variants are duplicated between newly-added data and old data. These are combined into unique data set: ", str_c(dupli_exon, collapse = ", ")))
    # when one has additional info and others does not, these can be merged. However, conflicts can be occured. It should call warnings.
    stop("append == T has not been supported yet.")
  }
  return(isodata)
}
