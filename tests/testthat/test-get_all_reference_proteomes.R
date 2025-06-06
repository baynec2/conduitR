test_that("we can download all reference proteomes", {
  # Note that this has been down in the past due to what seems to be service
  # interuptions to uniprots ftp server.
  expect_no_error(get_all_reference_proteomes())
  })
