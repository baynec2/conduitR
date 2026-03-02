test_that("plot_sample_cor_heatmap() works", {
  qf = readRDS(conduit_rds())@QFeatures
  qf = add_log_imputed_norm_assays(qf)

  plot

})
