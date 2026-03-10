test_that("plot_heatmap() works", {
  skip_if_not(file.exists(conduit_rds()), "conduit.rds fixture not available")
  qf = readRDS(conduit_rds())@QFeatures
  plot_heatmap(qf,
               "species",
               c("microbiome_treatment","day","tumor_classification"),
               "organism_type")
})
