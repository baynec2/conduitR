test_that("can retreive proteome ids for  human and b.theta", {
  out <- get_batch_of_proteome_ids(c("9606", "818"))

  expect_equal(unique(out$Organism), c(
    "Homo sapiens (Human)",
    "Bacteroides thetaiotaomicron"
  ))
})
