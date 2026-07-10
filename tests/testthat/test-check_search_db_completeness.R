# Offline unit tests for the search-database completeness gate (issue #20).
# These exercise the pure manifest-validation logic without any network.

make_manifest <- function(...) {
  tibble::tribble(
    ~proteome_id, ~resp_status, ~source,       ~n_sequences, ~expected,
    ...
  )
}

test_that("a complete manifest passes silently", {
  m <- make_manifest(
    "UP000000625", 200L, "uniprotkb", 4403L, 4403L,
    "UP001610837", 200L, "uniparc",   5022L, 5022L
  )
  expect_true(check_search_db_completeness(m))
})

test_that("an all-empty manifest is a hard error", {
  m <- make_manifest(
    "UP000000625", 500L, "not_downloaded", 0L, 4403L,
    "UP001610837", 500L, "not_downloaded", 0L, 5022L
  )
  # partial-failure warning fires first; assert the terminal error.
  expect_error(
    suppressWarnings(check_search_db_completeness(m)),
    "EMPTY"
  )
})

test_that("a grossly truncated total is a hard error (5 of 6064)", {
  # The #20 silent-truncation case: a UniParc fallback delivered 5 of 6064.
  m <- make_manifest(
    "UP001702503", 200L, "uniparc", 5L, 6064L
  )
  expect_error(
    check_search_db_completeness(m),
    "grossly incomplete"
  )
})

test_that("a partial failure warns but does not error when the total is adequate", {
  # One large complete proteome + one small failed one: total is well above the
  # 90% floor, so it should warn (not error).
  m <- make_manifest(
    "UP000000625", 200L, "uniprotkb",     4403L, 4403L,
    "UP999999999", 500L, "not_downloaded",   0L,   50L
  )
  expect_warning(
    expect_true(check_search_db_completeness(m)),
    "did not download completely"
  )
})

test_that("missing expected counts skip the shortfall check (endpoint count trusted)", {
  # expected = NA everywhere -> total_expected 0 -> shortfall check skipped;
  # non-zero delivery passes.
  m <- make_manifest(
    "UP000000625", 200L, "uniprotkb", 4403L, NA_integer_
  )
  expect_true(check_search_db_completeness(m))
})

test_that("the error message includes a per-proteome manifest", {
  m <- make_manifest(
    "UP001702503", 200L, "uniparc", 5L, 6064L
  )
  err <- tryCatch(check_search_db_completeness(m), error = function(e) conditionMessage(e))
  expect_match(err, "UP001702503")
  expect_match(err, "expected=6064")
  expect_match(err, "delivered=5")
})
