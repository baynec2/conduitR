# Validation works for good id
test_that("validation works with insulin", {
  expect_equal(validate_uniprot_accession_ids("P01308"), "P01308")
})

# Validation works for bad ID
test_that("validation removes bad IDs", {
  expect_message(validate_uniprot_accession_ids("NotanID"), "1 invalid UniProt IDs removed.")
})

# Validation works for good and bad ids
test_that("validation removes bad IDs", {
  expect_message(
    validate_uniprot_accession_ids(c("P01308", rep("NotanID", 100), "P01308")),
    "100 invalid UniProt IDs removed."
  )
})
