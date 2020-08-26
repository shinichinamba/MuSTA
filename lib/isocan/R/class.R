check_isodata_valid <- function(object){
  if(ncol(object@exp) - 1 != length(object@sample_group)){
    paste0("Unequal number of sample, group: ", ncol(object@exp) - 1, length(object@sample_group))
  }else if(distinguish_slot(object@exon_ranges) != "exon_ranges"){
    "The exon_ranges slot does not satisfy its requirements."
  }else if(distinguish_slot(object@variant) != "variant"){
    "The variant slot does not satisfy its requirements."
  }else TRUE
  #
}


#' S4 class for data storage
#'
#' @description You should NOT use this function for inputting your data because data stored in different slots may be inconsistent with each other.
#' Instead, you can use 'store_isodata' function to input and extend your data,
#' and next you can use 'construct_isodata' for consistent position data storage in exon, transcript, and gene levels.
#' For checking if your isodata object has consistent structure, use 'check_isodata_consistent'.
#'
#' @slot exon_ranges A tibble for exon-level information, including exons' positions and additional meta data.
#' Any positional data should be stored first in this slot, and higher-level positional data should be constructed by the function 'construct_isodata'.
#' @slot gene_attr A tibble for gene-level information, including transcripts' positions and additional meta data.
#' @slot tx_attr A tibble for transcript-level information, possibly including transcripts' positions and additional meta data.
#' @slot variant A tibble for variant information.
#' @slot tx_seq_correspondence A tibble with 'tx_id' and 'seq_id' columns.
#' @slot gene_exp A tibble for gene-level expression.
#' @slot tx_exp A tibble for transcript-level expression.
#' @slot exp_type 1-length character vector. For the time being, only "TPM" is allowed.
#' @slot sample_group A character vector. It must have the smaller length than colnames of tx_exp by one.
#' @export
isodata <- setClass("isodata",
                        slots = list(
                          exon_ranges = "tbl",
                          gene_attr = "tbl",
                          tx_attr = "tbl",
                          variant = "tbl",
                          tx_seq_correspondence = "tbl",
                          gene_exp = "tbl",
                          tx_exp = "tbl",
                          exp_type = "character",
                          sample_group = "character"
                        ),
                        prototype = list(
                          exp_type = "TPM"
                          #
                        ),
                        validity = check_isodata_valid)

#' append data to an isodata class object.
#'
#' @description "isodata" class object consists of many slots of tibbles.
#' It is aimed to be able to store almost all information about full-length transcripts.
#' "import_iso" is a S4 generic for importing data into an isodata class object.
#' @param isodata An isodata class object.
#' @param obj A data to be imported.
#' @param append Logical. If false, old data in the slot concerned will be deleted.
iso_append <- function(isodata, obj, append = T) {
  if(!is(isodata, "isodata")) stop("The first argument must be an isodata class object.")
  isodata_append(append_class(obj), isodata, append) #change arrangement for method dispatch.
}

isodata_append <- function(obj, isodata, append) {
  UseMethod("isodata_append", obj)
}

distinguish_slot <- function(obj){
  col_obj <- colnames(obj)
  attr <- c()
  if(is(obj, "tbl")){
    if(all(c("seqnames", "start", "end", "strand", "transcript_id") %in% col_obj))           attr <- c(attr, "exon_ranges")
    if(all(c("seq_id", "chr", "type", "start", "end", "ref", "alt", "sample") %in% col_obj)) attr <- c(attr, "variant")
    if(all(c("chr", "type", "start", "end", "ref", "alt", "sample") %in% col_obj))           attr <- c(attr, "variant_without_seq_id")
    if(all(c("seq_id", "transcript_id") %in% col_obj))                                       attr <- c(attr, "tx_seq_correspondence")

  }
  return(attr)
}

append_class <- function(obj){
  c_obj <- class(obj)
  if(!"GRanges" %in% c_obj) try(class(obj) <- c(c_obj, str_c("iso_", distinguish_slot(obj))), silent = T)
  return(obj)
}

is_isodata <- function(isodata_obj) is(isodata_obj, "isodata")

construct_isodata <- function(isodata_obj){
  if(!is_isodata(isodata_obj)) stop("The argument must be 'isodata' class object.")
  #
}

check_isodata_consistent <- function(isodata_obj){
  dplyr::if_else(isodata_obj == construct_isodata(isodata_obj), TRUE, FALSE, missing = FALSE)
}

# A dammy S3 class for methodDispatch.
#
setOldClass("iso_exon_ranges")
setOldClass("iso_variant")
setOldClass("iso_variant_without_seq_id")
