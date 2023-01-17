read_collapsed_group_txt <- function(file){
  readr::read_tsv(file,
                  col_names = c("PB_id", "seq_id")) %>%
    mutate(seq_id = map(seq_id, str_split, ",")) %>%
    unnest() %>% unnest()
}

#' @export
read_gtf <- function(file) {
  col_name <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame")
  rtracklayer::import(con = file, format = "gtf") %>%
    as_tibble() %>%
    mutate_if(is.factor, as.character) %>%
    select(-width) %>%
    rename(feature = "type", frame = "phase", seqname = "seqnames") %>%
    mutate(score = as.character(score), frame = as.character(frame)) %>%
    tidyr::replace_na(list(score = ".", frame = ".")) -> gtf
  gtf %>% select(!!!syms(c(col_name, setdiff(colnames(gtf), col_name))))
}

read_gtf_old <- function(file){
  col_name <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "group")
  readr::read_tsv(file,
                  col_names = col_name,
                  col_types = "ccciicccc",
                  comment = "#"
  ) %>%
    mutate(index = row_number()) -> tbl_file
  tbl_file %>%
    mutate(group = group %>% str_remove(';$') %>% str_split("; ")) %>%
    tidyr::unnest() %>%
    tidyr::separate(group, c("key", "value"), sep = '[:space:]') %>%
    mutate(value = value %>% str_remove('^"') %>% str_remove('"$')) %>%
    group_by(!!!syms(c(col_name[1:8], "key", "index"))) %>%
    summarise(value = paste(value, collapse = "; ")) %>%
    ungroup() %>%
    tidyr::spread(key, value) %>%
    arrange(index) %>%
    select(-index) -> tbl_file
  mandatory_cols <- c(col_name[1:8], "gene_id", "transcript_id")
  tbl_file %>%
    select(!!!syms(c(mandatory_cols, setdiff(colnames(tbl_file), mandatory_cols))))
}

#' @export
write_gtf <- function(x, path, rename_list = NULL, default = ".") {
  res <- try(write_gtf_default(x, path, rename_list, default))
  if (is(res, "try-error")) {
    warning("write_gtf() falled back to rtracklayer::export()")
    x <- dplyr::rename(x, !!!rename_list)
    y <- dplyr::rename(x, type = feature)
    y <- GenomicRanges::makeGRangesFromDataFrame(y)
    rtracklayer::export(y, path, "gtf")
  } # for unknown error in stringi...
  invisible(x)
}

write_gtf_default <- function(x, path, rename_list = NULL, default = "."){
  stopifnot(!(is.null(default) || is.na(default)))
  # opt_scipen <- getOption("scipen")
  # options(scipen = 999L)
  required_cols <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "gene_id", "transcript_id")
  if (nrow(x) > 0L) x <- tidyr::replace_na(x, list(strand = "*"))
  x <-
    x %>%
    dplyr::mutate_if(is.character, dplyr::recode, .missing = default) %>%
    dplyr::rename(!!!rename_list)
  deficient_cols <- dplyr::setdiff(required_cols, colnames(x))
  x <- purrr::reduce(deficient_cols, ~ mutate(.x, !!sym(.y) := default), .init = x)
  meta_cols <-
    c(
      "gene_id",
      "transcript_id",
      dplyr::setdiff(colnames(x), required_cols)
    )
  add_colname_core <- function(df, colname){
    if (nrow(df) > 0L) df[, colname] <- str_c(colname, ' "', df[[colname]], '";') else df$colname <- character()
    df
  }
  add_colname <- function(df, col_vec) purrr::reduce(col_vec, add_colname_core, .init = df)
  x %>%
    add_colname(meta_cols) %>%
    tidyr::unite(!!!as.list(meta_cols), col = ".meta", sep = " ") %>%
    select(!!!syms(c(required_cols[-(-1:0 + length(required_cols))], ".meta"))) %>%
    readr::format_tsv(col_names = FALSE, quote_escape = "double") %>%
    stringr::str_replace('^\"\"', '\"') %>%
    stringr::str_replace_all('([^\"])\"', "\\1") %>%
    stringr::str_remove("\n$") %>%
    stringr::str_split("\n") %>% `[[`(1L) %>%
    write(file = path)
  # options(scipen = opt_scipen)
  invisible(x)
}


