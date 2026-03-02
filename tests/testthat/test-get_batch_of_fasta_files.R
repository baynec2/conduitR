test_that("get_fasta_file works", {
  skip_if_offline()
  bt_proteome <- "UP000436858"
  expect_no_error(get_fasta_file(bt_proteome))
})

test_that("get_fasta_file errors", {
  skip_if_offline()
  wrong_id <- "notvalidid"
  expect_error(get_fasta_file(wrong_id))
})
