args <- commandArgs(TRUE)

IO_summary_file <- args[1]
inter_merge_gtf <- args[2]
sqanti_gtf <- args[3]

number_report_samples_out <- args[4]
number_report_merged_out <- args[5]

workers <- as.integer(args[6])

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(furrr)
})
plan(multicore, workers = workers)

#fastq length stats:
#seqkit stats */long*.fastq -t dna -j 1 -a -T > length_report.txt

#load
input_files <- read_csv(IO_summary_file, col_types = cols(.default = col_character()))

#functions
line_count <- function(path) system(paste0("wc -l ", path), intern = T)

line_nonvacant_count <- function(path) system(paste0("grep -c '.' ", path), intern = T) 
count_txt <- function(paths) future_map_chr(paths, line_nonvacant_count) %>% as.integer()

grep_count <- function(path, search_chr) system(paste0("grep -c '", search_chr, "' ", path), intern = T)
count_fasta <- function(paths) future_map_chr(paths, grep_count, search_chr = "^>") %>% as.integer()

grep_exclude_count <- function(path, search_chr) system(paste0("grep -c -v '", search_chr, "' ", path), intern = T)
count_fastq <- function(paths) future_map_chr(paths, grep_exclude_count, search_chr = "^@") %>% as.integer()
# count_sam <- count_fastq

get_paf_unique <- function(path){
  suppressWarnings(
    readr::read_tsv(path,
                    col_names = c("seq_id"),
                    col_types = "c_______________________")
  ) %>% 
    `[[`("seq_id") %>% 
    n_distinct()
}
count_paf_unique <- function(paths) future_map_int(paths, get_paf_unique)

get_gtf_unique <- function(path){
  suppressWarnings(
    readr::read_tsv(path,
                    col_names = c("attribute"),
                    col_types = "________c")
  ) %>% 
    `[[`("attribute") %>% 
    n_distinct()
}
count_gtf_unique <- function(paths) future_map_int(paths, get_gtf_unique)
null_if_null <- function(x, fun) if (is.null(x) || any(is.na(x))) NULL else fun(x)
#main 1
message("number_report_samples...")

number_report_samples <- 
  tibble(
    samples = input_files[["sample"]],
    .name_repair = "minimal"
  ) %>% 
  mutate(
    polished_hq = input_files[["long_read_hq"]] %>% count_fastq(),
    polished_lq = input_files[["long_read_lq"]] %>% null_if_null(count_fastq),
    lordec_hq = input_files[["lordec_hq"]] %>% count_fasta(),
    lordec_lq = input_files[["lordec_lq"]] %>% null_if_null(count_fasta),
    minimap2_hq_paf = input_files[["minimap2_hq_paf"]] %>% count_paf_unique(),
    minimap2_lq_paf = input_files[["minimap2_lq_paf"]] %>% null_if_null(count_paf_unique),
    qfilt_hq_paf = input_files[["qfilt_hq_paf"]] %>% count_paf_unique(),
    qfilt_lq_paf = input_files[["qfilt_lq_paf"]] %>% null_if_null(count_paf_unique),
    intra_merge = input_files[["intra_merge_gtf"]] %>% count_gtf_unique()
    )

#main 2
message("number_report_merged...")

number_report_merged <- 
  tibble(
    inter_merge = inter_merge_gtf %>% count_gtf_unique(),
    sqanti_filtered = sqanti_gtf %>% count_gtf_unique(),
    .name_repair = "minimal"
  )

message("write out...")
write_tsv(number_report_samples, number_report_samples_out)
write_tsv(number_report_merged, number_report_merged_out)