test_that("get_fasta_file works", {
  skip_if_offline()
  bt_proteome <- "UP000436858"
  expect_no_error(get_fasta_file(bt_proteome))
})

test_that("get_fasta_file reports an unusable id instead of throwing", {
  skip_if_offline()
  # Since the completeness gate (#20), get_fasta_file does not throw on a bad
  # or incomplete proteome: it reports the shortfall, writes no file, and lets
  # the caller's gate decide whether the overall database is usable.
  wrong_id <- "notvalidid"

  out <- expect_no_error(get_fasta_file(wrong_id))

  expect_equal(out$proteome_id, wrong_id)
  expect_equal(out$source, "not_downloaded")
  expect_equal(out$n_sequences, 0L)
  # A failed download must not leave a file behind.
  expect_false(file.exists(paste0(wrong_id, ".fasta")))
})
