####pipeline####
print_list <- function(x) {
  paste0(
    substitute(x), ": \n",
    str_c("  ", names(x), ": ", as.character(x), collapse = "\n")
  )
}

####perse config####
check_soft_ver <- function(soft_name, soft_path, ver_comm, ver_parse_fun) {
  res <- try(suppressWarnings(system2(soft_path, ver_comm, stdout = TRUE, stderr = TRUE)), silent = TRUE)
  if (class(res) == "try-error") {
    list(FALSE, str_glue("{soft_name} is not found in {soft_path}."))
  } else {
    return(list(TRUE, ver_parse_fun(res[1])))
  }
}

check_soft_path_and_get_ver <- function(softwares, soft, soft_ver_comm, soft_ver_parse_fun) {
  res <- pmap(list(softwares, soft, soft_ver_comm, soft_ver_parse_fun), check_soft_ver)
  is_valid <- map_lgl(res, `[[`, 1L)
  if (!all(is_valid)) stop(map_chr(res[!is_valid], `[[`, 2L) %>% str_c(collapse = "\n"))
  else map_chr(res, `[[`, 2L) %>% setNames(softwares)
}

# stop_for_warnings <- function(x, f, ...) tryCatch(f(x, ...), warning = function(e) stop(e))
check_no_na <- function(x, arrow_all_na = FALSE) {
  na_vec <- is.na(c(input_files[[x]]))
  if (any(na_vec)) {
    if (!arrow_all_na || !all(na_vec)) abort(str_glue("The '{y}' field of the config file must be specified for all samples listed in the 'general' field.", y = str_remove(x, "_[:digit:]+$")))
  }
} 

####plan####
plan_add <- function(condition, true, false = "", plan = "drake_txt") {
  plan_txt <- eval(parse(text = plan))
  # print(plan_txt)
  new_txt <- paste0(plan_txt, if (condition) true else false)
  assign(plan, new_txt, pos = 1L)
}

####make_qscripts####
