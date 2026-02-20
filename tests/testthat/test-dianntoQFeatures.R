test_that("diannToQFeatures works", {
  expect_no_error(conduitR::diannToQFeatures("inst/extdata/v2/diann.parquet"))
})
