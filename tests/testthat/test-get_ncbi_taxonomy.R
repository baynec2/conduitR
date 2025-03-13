test_that("this works with simple test set.", {
  # Human and B.theta
  human_btheta = c(9606,818)
  # getting the taxonomy
  out = get_ncbi_taxonomy(human_btheta)
  # Should have Human and b.theta.
  expect_contains(out$species,c("Homo sapiens","Bacteroides thetaiotaomicron"))
})
