test_that("this works", {
  tmp_dir <- tempdir()  # Automatically deleted after session
  proteome_id = "UP000270929"
  expect_no_error(get_fasta_file(proteome_id))
  unlink(tmp_dir)
})

# UP000050895 has no sequences in UniProtKB; we fall back to UniParc (may or may not have sequences).
test_that("proteome with empty UniProtKB returns valid source (uniparc or not_downloaded)", {
  tmp_dir <- tempdir()
  proteome_id <- "UP000050895"
  result <- get_fasta_file(proteome_id, tmp_dir)
  expect_true(result$source %in% c("uniparc", "not_downloaded"))
  unlink(tmp_dir)
})

test_that("this works with a problematic id",{
  # Create a temporary directory
  tmp_dir <- tempdir()  # Automatically deleted after session
  proteome_id = "UP000041370"
  # Define temporary file path
  expect_no_error(get_fasta_file(proteome_id,file.path(tmp_dir)))
  unlink(tmp_dir)
})
