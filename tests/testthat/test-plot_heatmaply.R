test_that("plot_heatmaply works", {
  qf = readRDS("tests/data/conduit.rds")@QFeatures

  plot_heatmaply(qf,"species","microbiome_treatment")
})
