test_that("this works", {
  proteome_id = "UP000270929"
  expect_no_error(get_fasta_file(proteome_id))
})

# For whatever reason, this id downloads but has no body. Making sure a message is found for this case.
test_that("this gives an error with problematic id", {
  proteome_id = "UP000050895"
  expect_message(get_fasta_file(proteome_id),"Failed to download FASTA for Proteome: UP000050895")
})
