Usage
=====

``` bash
./MuSTA.R -h
```

    usage: ./MuSTA.R [-h] [-v] -i INPUT -t GTF -g GENOME [-o OUTPUT] [-p PREFIX] [--thread THREAD] [--copy-input] [--no-lordec] [--no-short-read] [--use-lq]
                     [--map-qual MAP_QUAL] [--salmon-ref] [--keep] [--sqanti-filter SQANTI_FILTER] [--sqanti-qc SQANTI_QC]
                     [--lordec-build-SR-graph LORDEC_BUILD_SR_GRAPH] [--lordec-correct LORDEC_CORRECT] [--minimap2 MINIMAP2]
                     [--samtools SAMTOOLS] [--salmon SALMON] [--seqkit SEQKIT] [--os {6,7}] [--force] [--verbose] [--dry-run]
    
    MuSTA v0.0.3: Multi-Sample Transcriptome Assembly
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --version         show program's version number and exit
    
    INPUT:
      -i INPUT, --input INPUT
                            An comma-separated file describing paths to input files.
                            This file should contain at least four headers: "sample", "long_read_hq", "short_read_1", "short_read_2".
                            Four headers can be additionally specified: "cluster_report", "long_read_lq", "SJ", "mutation", "sv".
      -t GTF, --gtf GTF     A gtf file of a reference transcriptome
      -g GENOME, --genome GENOME
                            A fasta file of a reference genome
    
    OUTPUT:
      -o OUTPUT, --output OUTPUT
                            A directory where all outputs will be placed (default: output)
      -p PREFIX, --prefix PREFIX
                            Prefix of result files (default: musta)
    
    CONFIGULATION:
      --thread THREAD       Number of threads (default: 1)
      --copy-input          Copy input files to a temporary directory rather than using them as is.
                            This option is recommended when input files are stored in external
                            storage systems, and In/Out load is troublesome
      --no-lordec           Do NOT run LoRDEC in order to correct long-read sequences with RNA-seq reads
      --no-short-read       Use only long-read reads and do NOT use short-read RNA-seq reads. LoRDEC and salmon will be skipped.
      --use-lq              Activate the {long_read_lq} field of the input config file (for using polished_lq.fastq in addition to polished_hq.fastq)
      --map-qual MAP_QUAL   A threshold of mapping quality. reads with quality larger than this threshold are selected (default: 50)
      --salmon-ref          Run salmon against a reference gtf as well as against a MuSTA-derived gtf
      --keep                Keep intermediate files
    
    EXTERNAL SOFTWARES:
      This option group is need not to be specified unless you do not export external softwares in the Bash environment
    
      --sqanti-filter SQANTI_FILTER
                            Path to {sqanti_filter.py}
      --sqanti-qc SQANTI_QC
                            Path to {sqanti_qc.py}
      --lordec-build-SR-graph LORDEC_BUILD_SR_GRAPH
                            Path to {lordec-build-SR-graph}
      --lordec-correct LORDEC_CORRECT
                            Path to {lordec-correct}
      --minimap2 MINIMAP2   Path to {minimap2}
      --samtools SAMTOOLS   Path to {samtools}
      --salmon SALMON       Path to {salmon}
      --seqkit SEQKIT       [Deprecated] MuSTA does not require seqkit any more.
    
    MISCELLANEOUS:
      --os {6,7}            Redhat OS version for running qsub jobs.
                            This option is effective only on SHIROKANE super computer in Human Genome Center, Japan (default: 7)
      --force               Re-run the entire pipeline, rather than resuming from the suspended point of the last run.
      --verbose             Verbose mode
      --dry-run             Check the provided arguments but not run subprocesses. All reports will be created in the current directory
