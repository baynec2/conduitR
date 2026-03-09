test_that("this works with character and numeric ids", {
  skip_if_offline()
  tmp_dir <- withr::local_tempdir()
  temp_fp <- paste0(file.path(tmp_dir), "test.fasta")
  # Test both character and numeric input in one pass using the same small organism
  expect_no_error(download_fasta_from_organism_ids(c("818", 818)))
  unlink(temp_fp)
})


test_that("this doesn't work with a wrong id", {
  skip_if_offline()
  tmp_dir <- tempdir()
  temp_fp <- paste0(file.path(tmp_dir), "test.fasta")
  ids = c(818, "thisisawrongid")
  expect_error(download_fasta_from_organism_ids(ids))
  unlink(temp_fp)
})


test_that("this works with previously problematic ids", {
  skip_if_offline()
  tmp_dir <- tempdir()
  temp_fp <- paste0(file.path(tmp_dir), "test.fasta")
  problematic_id = c(818, 77133, 2320102)
  expect_no_error(download_fasta_from_organism_ids(problematic_id, fasta_destination_fp = temp_fp))
  unlink(temp_fp)
})
