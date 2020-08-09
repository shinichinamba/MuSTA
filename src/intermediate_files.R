salmon_dirname <- list()
salmon_dirname$pre_sqanti <- "salmon_init"
salmon_dirname$post_sqanti <- paste0(args$prefix, ".salmon")
salmon_dirname$ref <- "ref.salmon"

sqanti_prefix <- "initial_sqanti"
re_sqanti_prefix <- "sqanti"
fusion_sqanti_prefix <- paste0("fusion_", sqanti_prefix)

genome_name <- fs::path_ext_remove(fs::path_file(genome_fa))

file <- list()
file$genome_mmi <- path(dir$genome, genome_name, ext = "mmi")
file$genome_dir <- path(dir$genome, genome_name)
file$genome_seed <- path(dir$genome, genome_name, ext = "seed")
file$genome_pkg <- path(dir$genome, str_glue("BSgenome.{fa_name}.{user}.{pipeline_name}", fa_name = genome_name, user = "user", pipeline_name = .name)) #coupled with forgeBSgenome
#merge
file$inter_merge_gtf <- path(dir$merge, "inter_merge.gtf") #inter_merge.R
file$inter_merge_corr <- path(dir$merge, "inter_merge.correlation.txt") #inter_merge.R
file$inter_merge_fa <- path(dir$merge, "inter_merge.fa") #gtf_to_fasta
#salmon_pre_sqanti
file$salmon_index_pre_sqanti <- path(dir$merge, "salmon", "pre_sqanti_index") #salmon_index_pre_sqanti
file$salmon_tpm_pre_sqanti <- path(dir$merge, paste0(salmon_dirname$pre_sqanti, ".tpm.txt")) #salmon_quant_pre_sqanti
file$salmon_count_pre_sqanti <- path(dir$merge, paste0(salmon_dirname$pre_sqanti, ".count.txt")) #salmon_quant_pre_sqanti
#sqanti
file$sqanti_class <- paste0(path(dir$merge, sqanti_prefix), "_classification.txt") #sqanti_qc
file$inter_merge_corrected <- path(dir$merge, paste0("inter_merge_corrected.", c("gtf", "fasta", "faa"))) #gtf_to_fasta
file$sqanti_filter_supplement_folder <- path(dir$merge, "sqanti_filter_supplement") #sqanti_filter
file$sqanti_filter_result <- paste0(file$sqanti_class, "_filterResults.txt") #sqanti_filter
#sqanti filter result
file$sqanti_gtf <- path(dir$result, paste0(args$prefix, ".gtf")) #re_inter_merge
file$sqanti_corr <- path(dir$result, paste0(args$prefix, ".correlation.txt")) #re_inter_merge
file$sqanti_fa <- path(dir$result, paste0(args$prefix, ".fa")) #gtf_to_fasta_distinct
file$tx_count_unique <- path(dir$result, paste0(args$prefix, ".pbcount.txt"))
file$tx_existence <- path(dir$result, paste0(args$prefix, ".existence.txt"))
file$re_sqanti_corrected <- path(dir$merge, paste0(args$prefix, "_corrected.", c("gtf", "fasta", "faa"))) #gtf_to_fasta
#link original
file$variant_filtered <- path(dir$result, "variant.coincident.txt")
file$variant_filtered_somatic <- path(dir$result, "variant.coincident.somatic.txt")
file$original_range_multiex <- path(dir$result, "original_range.multiex.txt")
#salmon_post_sqanti
file$salmon_index_post_sqanti <- path(dir$merge, "salmon", "post_sqanti_index") #salmon_index_post_sqanti
file$salmon_duplicate_clusters_post_sqanti <- path(dir$merge, "salmon", "post_sqanti_index", "duplicate_clusters.tsv") #salmon_index_post_sqanti
file$salmon_tpm_post_sqanti <- path(dir$result, paste0(salmon_dirname$post_sqanti, ".tpm.txt")) #salmon_quant_post_sqanti
file$salmon_count_post_sqanti <- path(dir$result, paste0(salmon_dirname$post_sqanti, ".count.txt")) #salmon_quant_post_sqanti
#re_sqanti
file$re_sqanti_junc <- paste0(path(dir$result, re_sqanti_prefix), "_junctions.txt")

#ref (e.g. gencode)
file$ref_fa <- path(dir$result, "ref.fa")
file$salmon_index_ref <- path(dir$merge, "salmon", "ref_index") #salmon_index_ref
file$salmon_tpm_ref <- path(dir$result, paste0(salmon_dirname$ref, ".tpm.txt")) #salmon_quant_ref
file$salmon_count_ref <- path(dir$result, paste0(salmon_dirname$ref, ".count.txt")) #salmon_quant_ref
#fusion
file$chimera_binded <- path(dir$fusion, "fusion.paf.bind.txt") #bind_chimera
file$fusion_rds <- path(dir$fusion, "fusion.RDS")
file$fusion_gene_fragment_gtf <- path(dir$fusion, "fusion.gene_fragment.gtf")
file$fusion_sqanti_class <- paste0(path(dir$fusion, fusion_sqanti_prefix), "_classification.txt")
file$fusion_summary <- path(dir$result, "fusion_summary.txt")
file$fusion_gtf <- path(dir$result, "fusion.gtf")
file$fusion_sv <- path(dir$result, "fusion_sv.txt")
#report
file$number_report_samples <- path(dir$report, "number_report_samples.txt")
file$number_report_merged <- path(dir$report, "number_report_merged.txt")

file <- map(file, path_abs)