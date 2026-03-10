test_that("plot_biplot() works", {
  skip_if_not(file.exists(conduit_rds()), "conduit.rds fixture not available")
  qf = readRDS(conduit_rds())@QFeatures
  plot_biplot(qf,"species",color = "microbiome_treatment",
              facet_formula = ~day)
})
