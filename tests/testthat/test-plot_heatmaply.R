test_that("plot_heatmaply works", {
  qf = readRDS(conduit_rds())@QFeatures

  plot_heatmaply(qf,"species","microbiome_treatment")
})
