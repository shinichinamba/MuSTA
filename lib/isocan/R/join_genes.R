#' @export
filter_by_list <- function(x, y_list){
  purrr::reduce2(y_list, names(y_list), .init = x, .f = ~ filter(..1, (!! sym(..3) == ..2)))
}

#' @export
count_exon <- function(gtf_tbl, cols_exon_depending = c("start", "end"), col_feature = NULL, sort = FALSE){
  if (!is.null(col_feature)) filter_by_list(gtf_tbl, col_feature) -> gtf_tbl
  dplyr::setdiff(colnames(gtf_tbl), cols_exon_depending) -> cols_rest
  dplyr::count(gtf_tbl, !!!syms(cols_rest), sort = sort)
}

#' @export
gather_exon <- function(gtf_tbl, cols_start_end = c("start", "end"), into = "exon_structure",
                        drop = TRUE, original_into = "exon_nest",
                        col_feature = NULL, sep_exon = "-", sep_intron = " "){ #TODO expand: gather_intron
  if (!is.null(col_feature)) filter_by_list(gtf_tbl, col_feature) -> gtf_tbl
  dplyr::setdiff(colnames(gtf_tbl), cols_start_end) -> cols_rest
  gtf_tbl %>%
    dplyr::arrange(!!!syms(c(cols_rest, cols_start_end))) %>%
    dplyr::mutate(.dammy = into,
           .start_end = stringr::str_c(!!! dplyr::syms(cols_start_end), sep = sep_exon)) %>%
    group_by(!!!dplyr::syms(cols_rest)) -> gtf_tbl
  if (drop) {
    gtf_tbl %>%
      dplyr::summarise(!!into := stringr::str_c(!! dplyr::sym(".start_end"), collapse = sep_intron)) %>%
      dplyr::ungroup()
  }else{
    gtf_tbl %>%
      dplyr::summarise(!!into := stringr::str_c(!! sym(".start_end"), collapse = sep_intron),
                !!original_into := list(dplyr::tibble(!!!dplyr::syms(cols_start_end)))) %>%
      dplyr::ungroup()
  }
}

#' @export
remove_ends <- function(gathered_gtf, col_structure = "exon_structure"){
  gathered_gtf %>%
    dplyr::mutate(!! col_structure := (!! dplyr::sym(col_structure)) %>% stringr::str_remove("^[:digit:]+") %>% stringr::str_remove("[:digit:]+$"))
}

#' @export
separate_ends <- function(gathered_gtf, col_structure = "exon_structure", into = c("tx_start", "tx_end")){
  gathered_gtf %>%
    dplyr::mutate(!! into[1] := (!! dplyr::sym(col_structure)) %>% stringr::str_extract("^[:digit:]+"),
           !! into[2] := (!! dplyr::sym(col_structure)) %>% stringr::str_extract("[:digit:]+$"),
           !! col_structure := (!! dplyr::sym(col_structure)) %>% stringr::str_remove("^[:digit:]+") %>% stringr::str_remove("[:digit:]+$")) %>%
    dplyr::mutate_at(dplyr::vars(!!!into), as.integer)
}

#' @export
tx_join.intronic_include <- function(x, y, by = c("seqname", "strand"),
                                     cols_start_end.x = c("start", "end"), cols_start_end.y = c("start", "end"),
                                     col_feature.x = NULL, col_feature.y = NULL){ #TODO inner_join's arguments
  x %>%
    gather_exon(cols_start_end = cols_start_end.x, col_feature = col_feature.x, into = ".exon_structure") %>%
    separate_ends(col_structure = ".exon_structure", into = c(".tx_start.x", ".tx_end.x")) -> x
  y %>%
    gather_exon(cols_start_end = cols_start_end.y, col_feature = col_feature.y, into = ".exon_structure") %>%
    separate_ends(col_structure = ".exon_structure", into = c(".tx_start.y", ".tx_end.y")) -> y
  left_join(x, y, by = c(by, ".exon_structure")) %>%
    dplyr::filter(!!dplyr::sym(".tx_start.y") <= !!dplyr::sym(".tx_start.x") & !!dplyr::sym(".tx_end.x") <= !!dplyr::sym(".tx_end.y")) %>%
    dplyr::select(-c(!!!c(".exon_structure", ".tx_start.x", ".tx_end.x", ".tx_start.y", ".tx_end.y")))
}

#' @export
tx_join.edge_tolerate_match <- function(x, y, by = c("seqname", "strand"),
                                        cols_start_end.x = c("start", "end"), cols_start_end.y = c("start", "end"),
                                        col_feature.x = NULL, col_feature.y = NULL) {
  dplyr::bind_rows(
    tx_join.intronic_match(x, y, by, cols_start_end.x, cols_start_end.y, col_feature.x, col_feature.y),
    tx_join.singleexon_overlap(x, y, by, cols_start_end.x, cols_start_end.y, col_feature.x, col_feature.y)
  )
}

tx_join.intronic_match <- function(x, y, by = c("seqname", "strand"),
                                   cols_start_end.x = c("start", "end"), cols_start_end.y = c("start", "end"),
                                   col_feature.x = NULL, col_feature.y = NULL){ #TODO inner_join's arguments
  #warning("This function drops single-exon transcripts. Please use carefully.") #FIXME warning only if it is needed.
  x %>%
    gather_exon(cols_start_end = cols_start_end.x, col_feature = col_feature.x, into = ".exon_structure") %>%
    remove_ends(col_structure = ".exon_structure") %>%
    dplyr::filter(!!dplyr::sym(".exon_structure") != "-") -> x
  y %>%
    gather_exon(cols_start_end = cols_start_end.y, col_feature = col_feature.y, into = ".exon_structure") %>%
    remove_ends(col_structure = ".exon_structure") %>%
    dplyr::filter(!!dplyr::sym(".exon_structure") != "-") -> y
  dplyr::inner_join(x, y, by = c(by, ".exon_structure")) %>%
    dplyr::select(-c(!!dplyr::sym(".exon_structure")))
}

tx_join.singleexon_overlap <- function(x, y, by = c("seqname", "strand"),
                                       cols_start_end.x = c("start", "end"), cols_start_end.y = c("start", "end"),
                                       col_feature.x = NULL, col_feature.y = NULL){ #TODO inner_join's arguments
  #warning("This function drops multi-exon transcripts. Please use carefully.") #FIXME warning only if it is needed.
  x %>%
    gather_exon(cols_start_end = cols_start_end.x, col_feature = col_feature.x, into = ".exon_structure") %>%
    separate_ends(col_structure = ".exon_structure", into = c(".tx_start.x", ".tx_end.x")) %>%
    dplyr::filter(!!dplyr::sym(".exon_structure") == "-") %>% dplyr::select(-.exon_structure) -> x
  y %>%
    gather_exon(cols_start_end = cols_start_end.y, col_feature = col_feature.y, into = ".exon_structure") %>%
    separate_ends(col_structure = ".exon_structure", into = c(".tx_start.y", ".tx_end.y")) %>%
    dplyr::filter(!!sym(".exon_structure") == "-") %>% dplyr::select(-.exon_structure) -> y
  dplyr::left_join(x, y, by = by) %>%
    dplyr::filter(.tx_start.x <= .tx_end.y & .tx_start.y <= .tx_end.x) %>%
    dplyr::select(-c(.tx_start.x, .tx_end.x, .tx_start.y, .tx_end.y))
}

#' @export
tx_join.include <- function(x, y, by = c("seqname", "strand"), #TODO "by" only takes a character vector unlike *_join
                            cols_start_end.x = c("start", "end"), cols_start_end.y = c("start", "end"),
                            col_feature.x = NULL, col_feature.y = NULL, #TODO not included tx of x should not be dropped implicitly, but matched to NA explicitly. i.e. inner_join to left_join.
                            parallel = FALSE, verbose = FALSE, big_size = FALSE){
  if (length(cols_start_end.x) != 2 || length(cols_start_end.y) != 2) stop("cols_start_end must have length of 2.")
  if (!is.null(col_feature.x)) filter_by_list(x, col_feature.x) -> x
  if (!is.null(col_feature.y)) filter_by_list(y, col_feature.y) -> y
  if (parallel && !pacman::p_loaded("furrr")) {
    warning("Package 'furrr' is not loaded. 'parallel' option was changed to FALSE.")
    parallel <- FALSE
  }
  .map2_lgl <- if (parallel) {furrr::future_map2_lgl} else {purrr::map2_lgl}
  .pmap_lgl <- if (parallel) {furrr::future_pmap_lgl} else {purrr::pmap_lgl}

  #identifier
  dplyr::setdiff(colnames(x), c(by, cols_start_end.x)) -> cols_rest.x
  dplyr::setdiff(colnames(y), c(by, cols_start_end.y)) -> cols_rest.y
  #add index & make distinct
  if (verbose) message("add index & make distinct...")
  x %>% gather_exon(cols_start_end = cols_start_end.x, into = ".exon_structure.x",
                    drop = TRUE) -> x_structure
  x_structure %>%
    dplyr::mutate(.group.x = dplyr::group_indices(x_structure,
                                    !!!dplyr::syms(c(by, ".exon_structure.x")))
    ) -> x_structure
  x_structure %>%
    dplyr::select(!!!dplyr::syms(c(".group.x", ".exon_structure.x", by))) %>%
    dplyr::distinct(!!dplyr::sym(".group.x"), .keep_all = TRUE) %>%
    separate_ends(col_structure = ".exon_structure.x", into = c(".tx_start.x", ".tx_end.x")) -> x_distinct
  x_structure %>%
    dplyr::distinct(!!!dplyr::syms(c(cols_rest.x, by, ".group.x"))) -> x_correspondence #TODO as a stand-alone function getting distinct subset and a corresponding tbl.
  rm(x);rm(x_structure)
  .cols_start_end.y <- list(.start.y = cols_start_end.y[1], .end.y = cols_start_end.y[2])
  y %>% dplyr::rename(!!!.cols_start_end.y) %>% #rename for escaping name conflict
    gather_exon(cols_start_end = c(".start.y", ".end.y"), into = ".exon_structure.y",
                drop = FALSE, original_into = ".exon_nest.y") -> y_structure
  y_structure %>%
    dplyr::mutate(.group.y = dplyr::group_indices(y_structure %>% dplyr::select(-c(!!dplyr::sym(".exon_nest.y"))),
                                    !!!dplyr::syms(c(by, ".exon_structure.y")))
    ) -> y_structure
  y_structure %>%
    dplyr::select(!!!dplyr::syms(c(".group.y", ".exon_structure.y", ".exon_nest.y", by))) %>%
    dplyr::distinct(!!dplyr::sym(".group.y"), .keep_all = TRUE) %>%
    separate_ends(col_structure = ".exon_structure.y", into = c(".tx_start.y", ".tx_end.y")) -> y_distinct
  y_structure %>%
    dplyr::distinct(!!!dplyr::syms(c(cols_rest.y, ".group.y"))) -> y_correspondence
  rm(y);rm(y_structure)
  #multi-exon
  if (verbose) message("processing multi-exon transcripts...")
  if (!big_size) {
    #try 1time join
    x_distinct %>%
      dplyr::filter(!.exon_structure.x == "-") %>% #multi-exon
      dplyr::left_join(y_distinct, by = by) %>%
      dplyr::filter(!!dplyr::sym(".tx_start.y") <= !!dplyr::sym(".tx_start.x") & !!dplyr::sym(".tx_end.x") <= !!dplyr::sym(".tx_end.y")) %>% #necessary condition
      dplyr::filter(.map2_lgl(!!dplyr::sym(".exon_structure.y"), !!dplyr::sym(".exon_structure.x"), stringr::str_detect)) %>%
      dplyr::filter(.pmap_lgl(list(!!!dplyr::syms(c(".tx_start.x", ".tx_end.x", ".exon_nest.y"))),
                       ~ any(..3[[".start.y"]] <= ..1 & ..1 <= ..3[[".end.y"]]) &&
                         any(..3[[".start.y"]] <= ..2 & ..2 <= ..3[[".end.y"]]))) %>%
      dplyr::select(!!!dplyr::syms(c(".group.x", ".group.y"))) -> res_multi
  }else{
    if (verbose) message("big_size option is TRUE. Splitting data...")
    x_distinct %>%
      dplyr::filter(!.exon_structure.x == "-") %>% #multi-exon
      dplyr::group_split(!!!dplyr::syms(by)) %>%
      purrr::map_dfr(~ {
        .x %>% dplyr::left_join(y_distinct, by = by) %>%
          dplyr::filter(!!dplyr::sym(".tx_start.y") <= !!dplyr::sym(".tx_start.x") & !!dplyr::sym(".tx_end.x") <= !!dplyr::sym(".tx_end.y")) %>% #necessary condition
          dplyr::filter(.map2_lgl(!!dplyr::sym(".exon_structure.y"), !!dplyr::sym(".exon_structure.x"), stringr::str_detect)) %>%
          dplyr::filter(.pmap_lgl(list(!!!dplyr::syms(c(".tx_start.x", ".tx_end.x", ".exon_nest.y"))),
                           ~ any(..3[[".start.y"]] <= ..1 & ..1 <= ..3[[".end.y"]]) &&
                             any(..3[[".start.y"]] <= ..2 & ..2 <= ..3[[".end.y"]]))) %>%
          dplyr::select(!!!dplyr::syms(c(".group.x", ".group.y")))
      }) -> res_multi
  }
  #single-exon
  if (verbose) message("processing single-exon transcripts...")
  if (!big_size) {
    #try 1time join
    x_distinct %>%
      dplyr::filter(.exon_structure.x == "-") %>% #single-exon
      dplyr::left_join(y_distinct %>% dplyr::select(-(!!dplyr::sym(".exon_structure.y"))), by = by) %>%
      dplyr::filter(!!dplyr::sym(".tx_start.y") <= !!dplyr::sym(".tx_start.x") & !!dplyr::sym(".tx_end.x") <= !!dplyr::sym(".tx_end.y")) %>%  #necessary condition
      dplyr::filter(.pmap_lgl(list(!!!dplyr::syms(c(".tx_start.x", ".tx_end.x", ".exon_nest.y"))),
                       ~ any(..3[[".start.y"]] <= ..1 & ..1 <= ..3[[".end.y"]]) &&
                         any(..3[[".start.y"]] <= ..2 & ..2 <= ..3[[".end.y"]]))) %>%
      dplyr::select(!!!dplyr::syms(c(".group.x", ".group.y"))) -> res_single
  }else{
    if (verbose) message("big_size option is TRUE. Splitting data...")
    x_distinct %>%
      dplyr::filter(.exon_structure.x == "-") %>% #single-exon
      dplyr::group_split(!!!dplyr::syms(by)) %>%
      purrr::map_dfr(~ {
        .x %>%
          dplyr::left_join(y_distinct %>% dplyr::select(-(!!dplyr::sym(".exon_structure.y"))), by = by) %>%
          dplyr::filter(!!dplyr::sym(".tx_start.y") <= !!dplyr::sym(".tx_start.x") & !!dplyr::sym(".tx_end.x") <= !!dplyr::sym(".tx_end.y")) %>%  #necessary condition
          dplyr::filter(.pmap_lgl(list(!!!dplyr::syms(c(".tx_start.x", ".tx_end.x", ".exon_nest.y"))),
                           ~ any(..3[[".start.y"]] <= ..1 & ..1 <= ..3[[".end.y"]]) &&
                             any(..3[[".start.y"]] <= ..2 & ..2 <= ..3[[".end.y"]]))) %>%
          dplyr::select(!!!dplyr::syms(c(".group.x", ".group.y")))
      }) -> res_single
  }
  rm(x_distinct);rm(y_distinct)
  #bind
  dplyr::bind_rows(res_multi, res_single) %>%
    dplyr::left_join(x_correspondence, by = ".group.x") %>%
    dplyr::left_join(y_correspondence, by = ".group.y") %>%
    dplyr::select(-c(!!!dplyr::syms(c(".group.x", ".group.y"))))
}
