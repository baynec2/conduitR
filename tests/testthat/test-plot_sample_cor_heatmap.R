test_that("plot_sample_cor_heatmap() works", {
  qf = readRDS("tests/data/conduit.rds")@QFeatures
  qf = add_log_imputed_norm_assays(qf)

  plot

})
