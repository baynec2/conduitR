test_that("plot_biplot() works", {
  qf = readRDS(conduit_rds())@QFeatures
  plot_biplot(qf,"species",color = "microbiome_treatment",
              facet_formula = ~day)
})
