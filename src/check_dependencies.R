required_packages_cran <- 
  c("rlang", "dplyr", "tidyr", "stringi", "stringr", "purrr", "future", "furrr", "ggplot2",
    "remotes", "pkgload", "config", "fs", "drake", "argparse", "clisymbols")
required_packages_bioc <- 
  c("BSgenome", "Biostrings", "rtracklayer", "GenomicRanges", "GenomicFeatures")

not_installed_flag <- FALSE
not_installed_packages_cran <- c()
not_installed_packages_bioc <- c()

message_header <- "There are some packages to be installed. Please install them with the code below:"
message_cran <- ""
message_bioc <- ""

####CRAN####
for (p in required_packages_cran) {
  if (!requireNamespace(p, quietly = TRUE)) not_installed_packages_cran <- c(not_installed_packages_cran, p)
}

if (length(not_installed_packages_cran) > 0L) {
  not_installed_flag <- TRUE
  message_cran <- paste0('install.packages("', not_installed_packages_cran, '", dependencies = TRUE)', collapse = "\n")
}

####Bioconductor####
for (p in required_packages_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) not_installed_packages_bioc <- c(not_installed_packages_bioc, p)
}

if (length(not_installed_packages_bioc) > 0L) {
  not_installed_flag <- TRUE
  message_bioc_header <- 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")\n'
  message_bioc_main <- paste0('BiocManager::install("', not_installed_packages_bioc, '")', collapse = "\n")
  message_bioc <- paste0(message_bioc_header, message_bioc_main)
}

####Summary####
if (not_installed_flag) {
  stop_message <- paste(message_header, message_cran, message_bioc, sep = "\n")
  stop(stop_message)
}
