context("qsub")

test_that("dots_parser", {
  expect_equal(dots_parser("a", 4), "a\n4")
  expect_equal(dots_parser("a", c(2, 4), "z", sep_collapse = " "), "a 2 4 z")
})

#make_qsubfile: skip.
#verify_path: skip.

test_that("write_job", {
  testthat::expect_equal(write_job("temporary contents", fs::path(tempdir(), "job_test.txt"), recursive = TRUE, add_time = FALSE)[[1]],
                         fs::path(tempdir(), "job_test.txt"))#path name check only
  testthat::expect_error(write_job("temporary contents2", fs::path(tempdir(), "___not__exist_directory__", "job_test.txt"), recursive = FALSE, add_time = FALSE))
})


#only test in the HGC supercomputer environment.
#if (tryCatch(system("qstat", intern = T) %>% is.character(), error = function(e) FALSE)) {}

#qsub & qrecall ...