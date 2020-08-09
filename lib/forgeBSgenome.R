args <- commandArgs(TRUE)

pipeline_name <- args[1] #pipe
srcdir <- args[2] # abs_path to fafile_dir
dcf_path <- args[3] # path to write DCF (seed) file
pkgdir <- args[4] # abs_path to write package files

suppressPackageStartupMessages({
  library(stringr)
  library(Biostrings) # mandatory for BSgenome::forgeBSgenomeDataPkg
})

date <- Sys.Date()
user <- "user"
fa_name <- fs::path_file(srcdir) #hg38
suffix <- ".fa"
chroms <- fs::path_ext_remove(fs::path_file(fs::dir_ls(srcdir, type = "file")))
chroms_Rexp <- paste0("c('", str_c(chroms, collapse = "', '"), "')")
pkgname <- str_glue("BSgenome.{fa_name}.{user}.{pipeline_name}")
pkgpath <- fs::path(pkgdir, pkgname)

dcf <- str_glue("Package: {pkgname}
Title: Full genome sequences for {fa_name}
Description: Full genome sequences for {fa_name} as provided by {user} and stored in Biostrings objects.
Version: 0.0.1
Author: {user}
Maintainer: {user} <your_email_address>
License: Depend on your sequence data file
organism: {fa_name}
organism_biocview: {fa_name}
common_name: unknown
provider: {user}
provider_version: {pipeline_name}
release_date: {date}
release_name: {pipeline_name}
BSgenomeObjname: {fa_name}
seqs_srcdir: {srcdir}
seqnames: {chroms_Rexp}
seqfiles_suffix: {suffix}")

fileConn <- file(dcf_path)
writeLines(dcf, fileConn)
close(fileConn)

if (fs::dir_exists(pkgpath)) fs::dir_delete(pkgpath)
BSgenome::forgeBSgenomeDataPkg(dcf_path, destdir = pkgdir)
