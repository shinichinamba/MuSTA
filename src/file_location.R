#source files
src <- list()

src$utils.R <- "src/utils.R"
src$parse_config.R <- "src/parse_config.R"
src$update_input_files.R <- "src/update_input_files.R"
src$intermediate_files.R <- "src/intermediate_files.R"
src$make_qscript.R <- "src/make_qscripts.R"
src$plan.R <- "src/plan.R"
src$config.yml <- "src/config.yml"

lib <- list()
#self-contained packages
##self-made R packages
lib$isocan <- "lib/isocan"
lib$jobwatcher <- "lib/jobwatcher"
##imported
lib$interleave_fastq <- "lib/interleave-fastq/interleave-fastq"

#sub-processes
##minimap2
lib$quality_filter_paf_sam.R <- "lib/quality_filter_paf_sam.R"
##merge
lib$intra_merge.R <- "lib/merge_iso_intrasample.R"
lib$inter_merge.R <- "lib/merge_iso_intersample.R"
lib$merge_iso_commonpart.R <- "lib/merge_iso_commonpart.R"
##salmon
lib$merge_salmon_iso.R <- "lib/merge_salmon_iso.R"
##sqanti
lib$update_post_sqanti.R <- "lib/update_post_sqanti.R"
##others
lib$filter_variants.R <- "lib/filter_variants.R"
lib$report_number.R <- "lib/report_number.R"
lib$keep_major_isoform.R <- "lib/keep_major_isoform.R"
lib$link_original_range.R <- "lib/link_original_ranges.R"
##make BSgenome package from genome.fa
lib$multifa2singlefa.sh <- "lib/multifa2singlefa.sh"
lib$forgeBSgenome.R <- "lib/forgeBSgenome.R"
##gtf -> fa
lib$gtf2fa.R <- "lib/gtf2fa.R"
##fusion
lib$bind_chimera.R <- "lib/bind_chimera.R"
lib$parse_chimera.R <- "lib/parse_chimera.R"
lib$summarise_chimera.R <- "lib/summarise_chimera.R"