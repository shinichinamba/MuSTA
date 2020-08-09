
args <- commandArgs(TRUE)
IO_summary_file <- args[1]
quality_threshold <- as.integer(args[2]) #50L

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
})
#load
input_files <- read_csv(IO_summary_file, col_types = cols(.default = col_character()))

suppressWarnings({
  paf_hq <- input_files[["minimap2_hq_paf"]]
  paf_lq <- input_files[["minimap2_lq_paf"]]
  
  sam_hq <- input_files[["minimap2_hq_sam"]]
  sam_lq <- input_files[["minimap2_lq_sam"]]
  
  paf_out_hq <- input_files[["qfilt_hq_paf"]]
  paf_out_lq <- input_files[["qfilt_lq_paf"]]
  
  sam_out_hq <- input_files[["qfilt_hq_sam"]]
  sam_out_lq <- input_files[["qfilt_lq_sam"]]
})

#functions
read_paf_quality <- function(file){
  suppressWarnings(
    readr::read_tsv(file,
                    col_names = c("seq_id", "map_quality"),
                    col_types = "c__________i____________"))
} #query_start, query_end, start, end is 0-based


qfilt <- function(paf, sam, out_paf, out_sam, thres = quality_threshold){
  paf_quality <- read_paf_quality(paf)
  pass_id <- 
    paf_quality %>% 
    group_by(seq_id) %>% 
    summarise(quality_pass = all(map_quality > quality_threshold)) %>% 
    filter(quality_pass) %>% 
    `[[`("seq_id")
  
  readLines(paf)[paf_quality$seq_id %in% pass_id] %>% 
    write(out_paf)
  
  sam <- readLines(sam)
  is_sam_head <- str_detect(sam, "^@")
  sam_head <- sam[is_sam_head]
  sam_body <- sam[!is_sam_head]
  sam_id <- 
    sam_body %>% 
    str_split_fixed("\t", 2) %>% 
    `[`(, 1L)
  c(
    sam_head,
    sam_body[sam_id %in% pass_id]
  ) %>% 
    write(out_sam)
}

#main
invisible(pwalk(list(paf_hq, sam_hq, paf_out_hq, sam_out_hq), qfilt))
if (!(is.null(paf_lq) || is.null(sam_lq) || is.null(paf_out_lq) || is.null(sam_out_lq))) invisible(pwalk(list(paf_lq, sam_lq, paf_out_lq, sam_out_lq), qfilt))
