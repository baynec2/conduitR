test_that("this works with simple test set.", {
  # Human and B.theta
  human_btheta = c(9606,818)
  # getting the taxonomy
  out = get_ncbi_taxonomy(human_btheta)
  # Should have Human and b.theta.
  expect_contains(out$species,c("Homo sapiens","Bacteroides thetaiotaomicron"))
})
test_that("there are no NAs",{
  # Human and B.theta
  human_btheta <- c(9606, 818)

  # Getting the taxonomy
  out <- get_ncbi_taxonomy(human_btheta)

  # Ensure no missing values
  expect_false(any(is.na(out$domain)))
  expect_false(any(is.na(out$kingdom)))
  expect_false(any(is.na(out$phylum)))
  expect_false(any(is.na(out$class)))
  expect_false(any(is.na(out$order)))
  expect_false(any(is.na(out$family)))
  expect_false(any(is.na(out$genus)))
  expect_false(any(is.na(out$species)))

})

test_that("problematic IDs work",{
  problematic_ids = c(77133,2320102)

  out = get_ncbi_taxonomy(problematic_ids)
})
