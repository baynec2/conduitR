test_that("this works with character ids", {
  expect_no_error(download_fasta_from_organism_ids(c("818","9606")))
})


test_that("this works with numeric ids", {
  expect_no_error(download_fasta_from_organism_ids(c(818,9606)))
})


test_that("this doesn't work with a wrong id", {
  expect_error(download_fasta_from_organism_ids(c(818,"9606dfsaasd")))
})
