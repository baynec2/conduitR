test_that("can retreive proteome ids for  human and b.theta", {
  skip_if_offline()
  out <- get_proteome_ids_from_organism_ids(c("9606", "818"))

  expect_true(grepl("Homo sapiens", unique(out$organism)[1]))
  expect_true(grepl("Bacteroides thetaiotaomicron", unique(out$organism)[2]))
})
