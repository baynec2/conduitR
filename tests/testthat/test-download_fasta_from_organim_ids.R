test_that("this works with character ids", {
  tmp_dir <- tempdir()  # Automatically deleted after session
  temp_fp =paste0(file.path(tmp_dir),"test.fasta")
  ids = c("818","9606")
  expect_no_error(download_fasta_from_organism_ids(ids))
  unlink(temp_fp)
})


test_that("this works with numeric ids", {
  tmp_dir <- tempdir()  # Automatically deleted after session
  temp_fp =paste0(file.path(tmp_dir),"test.fasta")
  ids = c(818,9606)
  expect_no_error(download_fasta_from_organism_ids(ids))
  unlink(temp_fp)
})


test_that("this doesn't work with a wrong id", {
  tmp_dir <- tempdir()  # Automatically deleted after session
  temp_fp =paste0(file.path(tmp_dir),"test.fasta")
  ids = c(818,"thisisawrongid")
  expect_error(download_fasta_from_organism_ids(ids))
  unlink(temp_fp)
})


test_that("this works with a problematic id",{
  # Create a temporary directory
  tmp_dir <- tempdir()  # Automatically deleted after session
  problematic_id  = c(818,77133,2320102)
  # Define temporary file path
  temp_fp =paste0(file.path(tmp_dir),"test.fasta")
  expect_no_error(download_fasta_from_organism_ids(problematic_id,destination_fp = temp_fp))
  unlink(temp_fp)
})
