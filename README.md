# MuSTA: Multi-Sample Transcript Assembly

[![Documentation Status](https://readthedocs.org/projects/musta/badge/?version=latest)](https://musta.readthedocs.io/en/latest/?badge=latest)

This is a pipeline and associated R codes which take long-read isoform sequencing results as input and genarate an assembled transcriptome.

You can see the documents and the quick-start guides at https://musta.readthedocs.io.

A long-read transcriptome simlator, simlady, can be found at https://github.com/shinichinamba/simlady.


## Usage

`./MuSTA.R --help`

## Citation

If you feel this pipeline (or related codes) is useful, please cite:

S Namba *et al*. Transcript-targeted analysis reveals isoform alterations and double-hop fusions in breast cancer. ***communications biology*** (2021) https://www.nature.com/articles/s42003-021-02833-4


## Update

### v0.0.4

* Bug fixed


### v0.0.3

* Allowing gzipped inputs
* Allowing both fasta and fastq files for short-read RNA-seq reads
* New --no-short-read option
* Bug fixed
