test_that("diann_to_qfeatures works", {
  expect_no_error(diann_to_qfeatures("inst/extdata/diann.parquet"))
})
