test_that("plot_biplot() works", {
  qf = readRDS("tests/data/conduit.rds")@QFeatures
  plot_biplot(qf,"species",color = "microbiome_treatment",
              facet_formula = ~day)
})
