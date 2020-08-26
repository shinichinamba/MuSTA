verify_scalar <- function(...) purrr::walk(list(...), ~ assertthat::assert_that(assertthat::is.scalar(.x)))
verify_no_na <- function(x) assertthat::assert_that(assertthat::noNA(x))
verify_in <- function(x, y) assertthat::assert_that(all(x %in% y), msg = paste0(deparse(x), " is NOT included in ", deparse(y), "."))
verify_file_exists <- function(x) assertthat::assert_that(fs::file_exists(x))
  
verify_path <- function(path, recursive) {
  if (recursive) {
    assertthat::assert_that(is.character(path))
    fs::dir_create(fs::path_dir(path))
  }else{
    assertthat::assert_that(fs::dir_exists(fs::path_dir(path)))
  }
}

verify_hgc <- function() {
  assertthat::assert_that(get_jobwatcher_mode() == "hgc", msg = "This function works only in HGC super computer.")
}
verify_uge <- function() {
  assertthat::assert_that(get_jobwatcher_mode() %in% c("hgc", "uge"), msg = "This function works only in UGE environment.")
}

.as_character <- function(x) vctrs::vec_cast(x, character())
.map_as_character <- function(...) {
  x <- rlang::list2(...)
  purrr::map(x, .as_character)
}

is_number <- function(x) stringr::str_detect(x, "^[:digit:]+$")

assign_oneof_default_arg_chr <- function(arg_vec) {
  default_args <- formals(sys.function(sys.parent()))
  arg_vec <- arg_vec[arg_vec %in% names(default_args)]
  for (i in seq_along(arg_vec)) {
    true_arg <- .as_character(eval(as.name(arg_vec[i]), envir = parent.frame(1L)))[1L]
    verify_in(true_arg, eval(default_args[[arg_vec[i]]]))
    assign(arg_vec[i], true_arg, envir = parent.frame(1L))
  }
  invisible(arg_vec)
}

# .tester0 <- function() print(eval(as.name("a"), envir = parent.frame(1)))
# .tester1 <- function(a = "1", b = c("1","2","3")) {
#   assign_oneof_default_arg_chr(c("a", "b"))
#   # .tester0()
# }
# .tester2 <- function(a = 3) .tester1()
# .tester2_2 <- function(a = 3) .tester1(a = "4", b = "2")
# .tester3 <- function() .tester2()