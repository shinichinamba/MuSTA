## import a GRanges object into an isodata class object.
##
## @inheritParams import_iso
#isodata_append.vcfR <- funciton(isodata, obj, append){
#
#}

read_MANOmutect <- function(path, sample_name){
  readr::read_tsv(path, comment = "##") %>%
    rename(
      chr = contig,
      end = position,
      ref = ref_allele,
      alt = alt_allele) %>%
    mutate(start = end - 1, #point mutation only.
           type = "*",
           sample = sample_name,
           somatic = NA,
           ref = str_to_lower(ref),
           alt = str_to_lower(alt)) #for the time being.
}

read_MANOsid <- function(path, sample_name){
  readr::read_tsv(path, comment = "##") %>%
    rename(
      chr = `#CHROM`,
      start = POS,
      ref = REF,
      alt = ALT,
      somatic = INFO) %>%
    mutate(ref = str_sub(ref, start = 2L) %>% str_to_lower(),
           alt = str_sub(alt, start = 2L) %>% str_to_lower(),
           end = start + str_length(ref),
           somatic = (somatic == "SOMATIC"),
           type = sapply(str_length(alt) - str_length(ref),
                         function(x) if(x > 0){"+"}else if(x < 0){"-"}else{stop("The argument was supposed to contain only Ins/Del, but does not.")}),
           sample = sample_name)
}

#'@export
read_genomon_variant <- function(path){
  readr::read_tsv(path, col_types = readr::cols_only(Chr = "c", Start = "i", End = "i", Ref = "c", Alt = "c"), comment = "#") %>%
    `colnames<-`(c("chr", "start", "end", "ref", "alt")) %>%
    mutate(start = if_else(ref == "-", start, start - 1L),
           #End = End,
           ref = str_to_lower(if_else(ref == "-", NA_character_, ref)),
           alt = str_to_lower(if_else(alt == "-", NA_character_, alt)),
           type = if_else(is.na(ref), "+", if_else(is.na(alt), "-", "*")))
}
