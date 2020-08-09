#tx.gtf -> tx.fa
args <- commandArgs(TRUE)

gtf <- args[1]
genome_pkg <- args[2]
out_fa <- args[3]

suppressPackageStartupMessages({
  pkgname <- fs::path_file(genome_pkg)
  pkgload::load_all(genome_pkg, export_all = FALSE)
  genome <- eval(parse(text = paste0(pkgname, "::", pkgname)))
  library(GenomicFeatures)
  library(GenomicRanges)
})

gtf_gr <- rtracklayer::import.gff(gtf)
types <- unique(gtf_gr$type)
if (length(types) > 1 && "exon" %in% types) gtf_gr <- gtf_gr[gtf_gr$type == "exon"]
gtf_gr <- split(gtf_gr, gtf_gr$transcript_id)
strands <- as.character(runValue(strand(gtf_gr)))
gtf_gr <- 
  c(
    sort(gtf_gr[strands == "+"]),
    sort(gtf_gr[strands == "-"], decreasing = TRUE)
  )
tx_seq <- extractTranscriptSeqs(genome, gtf_gr)
Biostrings::writeXStringSet(tx_seq, out_fa)