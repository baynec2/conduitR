test_that("this works", {
  tmp_dir <- tempdir()  # Automatically deleted after session
  proteome_id = "UP000270929"
  expect_no_error(get_fasta_file(proteome_id))
  unlink(tmp_dir)
})

# For whatever reason, this id downloads but has no body. Making sure a message is found for this case.
test_that("this gives an error with problematic id", {
  tmp_dir <- tempdir()  # Automatically deleted after session
  proteome_id = "UP000050895"
  expect_message(get_fasta_file(proteome_id),"Failed to download FASTA for Proteome: UP000050895")
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
