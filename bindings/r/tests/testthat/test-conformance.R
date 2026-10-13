# Replay the shared conformance cases in spec/conformance/ (also replayed by the Python suite).
# The spec directory is outside the package, so these tests skip when the package is checked
# from a built tarball.

spec_path <- function(name) test_path("..", "..", "..", "..", "spec", "conformance", name)

test_that("RNG stream seeds match the shared golden cases", {
  path <- spec_path("rng_streams.json")
  skip_if_not(file.exists(path), "spec/conformance is not available")
  cases <- jsonlite::read_json(path)$cases
  for (case in cases) {
    expect_identical(stream_seed(case$seed, case$path), case$expected)
  }
})
