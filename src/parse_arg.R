parser <- ArgumentParser(description = paste0(.name, " v", .version, ': Multi-Sample Transcriptome Assembly'))

parser$add_argument('-v', '--version', action = 'version', version = .version)

#### INPUT ####
group_input <- parser$add_argument_group('INPUT')

group_input$add_argument(
  '-i', '--input', required = TRUE,
  help = 'An comma-separated file describing paths to input files. This file should contain at least four headers: "sample", "long_read_hq", "short_read_1", "short_read_2". Four headers can be additionally specified: "cluster_report", "long_read_lq", "SJ", "mutation", "sv".'
)

group_input$add_argument(
  '-t', '--gtf', required = TRUE,
  help = 'A gtf file of a reference transcriptome'
)

group_input$add_argument(
  '-g', '--genome', required = TRUE,
  help = 'A fasta file of a reference genome'
)

#### OUTPUT ####
group_output <- parser$add_argument_group('OUTPUT')

group_output$add_argument(
  '-o', '--output', default = 'output',
  help = 'A directory where all outputs will be placed (default: output)'
)

group_output$add_argument(
  '-p', '--prefix', default = 'musta',
  help = 'Prefix of result files (default: musta)'
)

#### CONFIGULATION ####
group_configulation <- parser$add_argument_group('CONFIGULATION')

group_configulation$add_argument(
  '--thread', type = "integer", default = 1L, 
  help = 'Number of threads (default: 1)'
)

group_configulation$add_argument(
  '--copy-input', action = "store_true",
  help = 'Copy input files to a temporary directory rather than using them as is. This option is recommended when input files are stored in external storage systems, and In/Out load is troublesome'
)

group_configulation$add_argument(
  '--no-lordec', action = "store_true",
  help = 'Do NOT run LoRDEC in order to correct long-read sequences with RNA-seq reads'
)

group_configulation$add_argument(
  '--no-short-read', action = "store_true",
  help = 'Use only long-read reads and do NOT use short-read RNA-seq reads. LoRDEC and salmon will be skipped.'
)


group_configulation$add_argument(
  '--use-lq', action = "store_true",
  help = 'Activate the {long_read_lq} field of the input config file (for using polished_lq.fastq in addition to polished_hq.fastq)'
)

group_configulation$add_argument(
  '--map-qual', type = 'integer', default = 50L,
  help = 'A threshold of mapping quality. reads with quality larger than this threshold are selected (default: 50)'
)

group_configulation$add_argument(
  '--salmon-ref', action = "store_true",
  help = 'Run salmon against a reference gtf as well as against a MuSTA-derived gtf'
)

group_configulation$add_argument(
  '--keep', action = "store_true", 
  help = 'Keep intermediate files'
)

#### EXTERNAL SOFTWARES ####
group_ext <- parser$add_argument_group('EXTERNAL SOFTWARES', 
                                      'This option group is need not to be specified unless you do not export external softwares in the Bash environment')
group_ext$add_argument(
  '--sqanti-filter', default = 'sqanti_filter.py',
  help = 'Path to {sqanti_filter.py}'
)

group_ext$add_argument(
  '--sqanti-qc', default = 'sqanti_qc.py',
  help = 'Path to {sqanti_qc.py}'
)

group_ext$add_argument(
  '--lordec-build-SR-graph', default = 'lordec-build-SR-graph',
  help = 'Path to {lordec-build-SR-graph}'
)

group_ext$add_argument(
  '--lordec-correct', default = 'lordec-correct',
  help = 'Path to {lordec-correct}'
)

group_ext$add_argument(
  '--minimap2', default = 'minimap2',
  help = 'Path to {minimap2}'
)

group_ext$add_argument(
  '--samtools', default = 'samtools',
  help = 'Path to {samtools}'
)

group_ext$add_argument(
  '--salmon', default = 'salmon',
  help = 'Path to {salmon}'
)
group_ext$add_argument(
  '--seqkit', default = 'seqkit',
  help = '[Deprecated] MuSTA does not require seqkit any more.'
)

#### Miscellaneous ####
group_miscellaneous <- parser$add_argument_group('MISCELLANEOUS')

group_miscellaneous$add_argument(
  '--os', dest = 'os_version', choices = c("6", "7"), default = '7',
  help = 'Redhat OS version for running qsub jobs. This option is effective only on SHIROKANE super computer in Human Genome Center, Japan (default: 7)'
)

group_miscellaneous$add_argument(
  '--force', action = "store_true",
  help = 'Re-run the entire pipeline, rather than resuming from the suspended point of the last run.'
)

group_miscellaneous$add_argument(
  '--verbose', action = "store_true", 
  help = 'Verbose mode'
)

group_miscellaneous$add_argument(
  '--dry-run', action = "store_true",
  help = 'Check the provided arguments but not run subprocesses. All reports will be created in the current directory'
)

args <- parser$parse_args()
if (args$no_short_read) {
  args$no_lordec <- TRUE
  args$salmon_ref <- FALSE
}
