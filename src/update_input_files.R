input_files <- input_files %>%
  mutate(sample_dir = sample_dirs)

if (args$copy_input == TRUE) {
  input_files <- input_files %>%
    mutate(
      long_copy_hq = path(sample_dir, "long_hq.fastq"),
      long_copy_lq = path(sample_dir, "long_lq.fastq"),
      short_copy_1 = path(sample_dir, "short_1.fastq"),
      short_copy_2 = path(sample_dir, "short_2.fastq")
    )
} else {
  input_files <- input_files %>%
    mutate(
      long_copy_hq = long_read_hq,
      long_copy_lq = long_read_lq,
      short_copy_1 = short_read_1,
      short_copy_2 = short_read_2
    )
}

if (args$no_lordec == FALSE) {
  input_files <- input_files %>%
    mutate(interleave = path(sample_dir, "interleave.fastq"),
           interleave_lordectmp = paste0(interleave, "_k19_s3.h5"),
           lordec_hq = path(sample_dir, "hq.lordec.fa"),
           lordec_lq = path(sample_dir, "lq.lordec.fa")
    )
} else {
  input_files <- input_files %>%
    mutate(lordec_hq = path(sample_dir, "hq.nolordec.fa"), # TODO fastq2fasta
           lordec_lq = path(sample_dir, "lq.nolordec.fa")  # TODO fastq2fasta
    )
}

input_files <- input_files %>%
  mutate(merged_fa = path(sample_dir, "merged.lordec.fa"),
         minimap2_hq_sam = path(sample_dir, "hq.lordec.sam"),
         minimap2_lq_sam = path(sample_dir, "lq.lordec.sam"),
         minimap2_hq_paf = path(sample_dir, "hq.lordec.paf"),
         minimap2_lq_paf = path(sample_dir, "lq.lordec.paf"),
         qfilt_hq_paf = path(sample_dir, "hq.lordec.qfilt.paf"),
         qfilt_lq_paf = path(sample_dir, "lq.lordec.qfilt.paf"),
         qfilt_hq_sam = path(sample_dir, "hq.lordec.qfilt.sam"),
         qfilt_lq_sam = path(sample_dir, "lq.lordec.qfilt.sam"),
         samtools_hq_bam = path(sample_dir, "hq.lordec.qfilt.sort.bam"),
         samtools_lq_bam = path(sample_dir, "lq.lordec.qfilt.sort.bam"),
         merged_sam = path(sample_dir, "merged.lordec.qfilt.sam"),
         samtools_bam = path(sample_dir, "merged.lordec.qfilt.sort.bam"),
         chimera_paf = path(sample_dir, "lordec.qfilt.fusion.txt"),
         intra_merge_gtf = path(sample_dir, "intra_merge.gtf"),
         intra_merge_corr = path(sample_dir, "intra_merge.correlation.txt"),
         variant_local = path(sample_dir, "variant.txt"),
         fl_count = path(sample_dir, "fl_count.txt"),
         salmon_exp_pre_sqanti = path(sample_dir, salmon_dirname$pre_sqanti, "quant.sf"),
         salmon_exp_post_sqanti = path(sample_dir, salmon_dirname$post_sqanti, "quant.sf"),
         salmon_exp_ref = path(sample_dir, salmon_dirname$ref, "quant.sf")
  )

if (!args$use_lq) {
  input_files <- input_files %>%
    select(-c(long_copy_lq, lordec_lq, minimap2_lq_sam, minimap2_lq_paf,
              qfilt_lq_paf, qfilt_lq_sam, samtools_lq_bam)) %>% 
    mutate(
      merged_fa = lordec_hq,
      merged_sam = qfilt_hq_sam,
      samtools_bam = samtools_hq_bam
    )
}

path_abs2 <- function(x) if_else(x == "", "", as.character(path_abs(x)))
input_files <- 
  bind_cols(
    input_files %>% select(sample),
    input_files %>% select(-sample) %>% mutate_all(path_abs2)
  )
