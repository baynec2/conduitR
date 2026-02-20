test_that("create_conduit_obj() works", {
expect_no_error({
  conduit_obj = create_conduit_obj(QFeatures.rds = "",
                                 combined_metrics.tsv = "tests/data/conduit_output/combined_metrics.tsv",
                                 database_protein_taxonomy.tsv ="tests/data/conduit_output/database_taxonomy.tsv",
                                 detected_protein_taxonomy.tsv ="tests/data/conduit_output/detected_protein_taxonomy.tsv")


saveRDS(conduit_obj,"tests/data/conduit.rds")

})
  })
