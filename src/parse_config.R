####general####
mut_sv_mode <- "genomon" #TODO VCF or Genomon

####software####
softwares <- 
  c("sqanti_filter.py", "sqanti_qc.py", "lordec-build-SR-graph", "lordec-correct", 
    "minimap2", "samtools", "salmon")
software_path_argnames <- 
  c("sqanti_filter", "sqanti_qc", "lordec_build_SR_graph", "lordec_correct", 
    "minimap2", "samtools", "salmon")
soft_ver_comm <- 
  c("-v", "-v", "", "", "--version", "--version", "--version")
soft_ver_parse_fun <- 
  list(
    (function(x) str_remove(x, "^SQANTI ")),
    (function(x) str_remove(x, "^SQANTI ")),
    (function(x) str_remove(x, "^LoRDEC v")),
    (function(x) str_remove(x, "^LoRDEC v")),
    (function(x) x),
    (function(x) str_remove(x, "^samtools ")),
    (function(x) str_remove(x, "^salmon "))
  )
soft_names <- softwares %>% str_to_lower() %>% str_replace_all("-", "_")
stopifnot(length(softwares) == length(soft_ver_comm))
stopifnot(length(softwares) == length(soft_ver_parse_fun))

soft <- args[software_path_argnames] %>% setNames(soft_names)
soft_ver <- check_soft_path_and_get_ver(softwares, soft, soft_ver_comm, soft_ver_parse_fun)

####que_requirement####
process_names <- 
  c("fasta2BSgenome", "lordec_build", "lordec", "minimap2", "qual_filter", "samtools", "intra_merge", "inter_merge",
    "gtf2fasta", "salmon_index", "salmon_quant", "merge_salmon", "sqanti", 
    "filter_variant", "link_original_range", "fusion_bind", "fusion_parse", "fusion_sqanti", "fusion_summary", "report_number")
que_req <- 
  map(process_names, ~ config::get(config = ..1, file = src$config.yml)) %>% 
  setNames(process_names)

####input####
output_dir <- args$output %>% str_remove("/$") %>% fs::path_expand() %>% fs::path_abs()
ref_gtf <- args$gtf %>% fs::path_expand() %>% fs::path_abs()
genome_fa <- args$genome %>% fs::path_expand() %>% fs::path_abs()

result_dir <- path(output_dir, "result")
merge_dir <- path(output_dir, "merge")
script_dir <- path(output_dir, "script")

####input files####
input_files <- read_csv(
  args$input, 
  col_types = cols(.default = col_character())
)

possible_cols <- c("sample", "long_read_hq", "long_read_lq", "short_read_1", "short_read_2", "cluster_report", "SJ", "mutation", "sv")
selected_cols <- intersect(colnames(input_files), possible_cols)
allna_cols <- setdiff(possible_cols, selected_cols)

input_files <- input_files[selected_cols]

for (col in allna_cols) {
  input_files[[col]] <- rep(NA_character_, nrow(input_files))
}

samples <- input_files[["sample"]]

# check
mandatory_cols <- c("sample", "long_read_hq")
if (args$no_lordec == FALSE) mandatory_cols <- c(mandatory_cols, "short_read_1", "short_read_2")
if (args$use_lq == TRUE) mandatory_cols <- c(mandatory_cols, "long_read_lq")

missing_cols <- setdiff(mandatory_cols, colnames(input_files))
if (length(missing_cols) > 0L) {
  abort(paste0("Some mandatory columns are missing in the input csv file: ", str_c(missing_cols, collapse = ", ")), "requirement error")
}

dup_sample_id <- samples[duplicated(samples)]
if (length(dup_sample_id) > 0L) {
  abort(paste0("All sample identifiers in the input csv file must be unique: ", str_c(dup_sample_id, collapse = ", ")), "requirement error")
}

if (!args$no_lordec) {
  walk2(c("long_read_hq", "short_read_1", "short_read_2", "cluster_report", "SJ"), c(FALSE, FALSE, FALSE, TRUE, TRUE), check_no_na)
} else {
  walk2(c("long_read_hq", "cluster_report", "SJ"), c(FALSE, TRUE, TRUE), check_no_na)
}
if (args$use_lq) check_no_na("long_read_lq")
