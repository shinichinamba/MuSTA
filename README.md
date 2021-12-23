# MuSTA: Multi-Sample Transcript Assembly

[![Documentation Status](https://readthedocs.org/projects/musta/badge/?version=latest)](https://musta.readthedocs.io/en/latest/?badge=latest)

This is a pipeline and associated R codes which take long-read isoform sequencing results as input and genarate an assembled transcriptome.

You can see the documents and the quick-start guides at https://musta.readthedocs.io.

A long-read transcriptome simlator, simlady, can be found at https://github.com/shinichinamba/simlady.


## Usage

`./MuSTA.R --help`

## Citation

If you feel this pipeline (or related codes) is useful, please cite:

S Namba *et al*. Multi-sample Full-length Transcriptome Analysis of 22 Breast Cancer Clinical Specimens with Long-Read Sequencing. ***BioRxiv*** (2020) https://doi.org/10.1101/2020.07.15.199851


## Update

### v0.0.3

* Allowing gzipped inputs
* Allowing both fasta and fastq files for short-read RNA-seq reads
* New --no-short-read option
* Bug fixed
