###unexported object derived from usethis###
todo_bullet <- function() crayon::red(clisymbols::symbol$bullet)
done_bullet <- function() crayon::green(clisymbols::symbol$tick)
fail_bullet <- function() crayon::bgRed(clisymbols::symbol$cross)

bulletize <- function(line, bullet) paste0(bullet, " ", line)
todo <- function(..., .envir = parent.frame()) {
  out <- stringr::str_glue(..., .envir = .envir)
  cli::cat_line(bulletize(out, bullet = todo_bullet()))
}

done <- function(..., .envir = parent.frame()) {
  out <- stringr::str_glue(..., .envir = .envir)
  cli::cat_line(bulletize(out, bullet = done_bullet()))
}

fail <- function(..., .envir = parent.frame()) {
  out <- stringr::str_glue(..., .envir = .envir)
  cli::cat_line(bulletize(out, bullet = fail_bullet()))
}
