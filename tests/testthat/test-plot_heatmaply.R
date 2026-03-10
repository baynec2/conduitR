test_that("plot_heatmaply works", {
  skip_if_not(file.exists(conduit_rds()), "conduit.rds fixture not available")
  qf = readRDS(conduit_rds())@QFeatures

  plot_heatmaply(qf,"species","microbiome_treatment")
})
