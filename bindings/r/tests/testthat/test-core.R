test_that("the core version is a version string", {
  expect_match(core_version(), "^\\d+\\.\\d+\\.\\d+")
})

test_that("stream seeds are deterministic", {
  path <- list("subject", 17, "ae")
  expect_identical(stream_seed(1, path), stream_seed(1, path))
  expect_false(identical(stream_seed(1, path), stream_seed(2, path)))
  expect_identical(stream_seed(1, list(17L)), stream_seed(1, list(17)))
})

test_that("bad inputs are rejected", {
  expect_error(stream_seed(-1, list()), "whole number")
  expect_error(stream_seed(1, list(1.5)), "whole number")
  expect_error(stream_seed(1, list(TRUE)), "string or a whole number")
})
