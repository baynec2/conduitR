test_that("Function works", {
  skip_if_offline()
  expect_no_error(get_eggnog_function("ENOG502S5P5"))
})
