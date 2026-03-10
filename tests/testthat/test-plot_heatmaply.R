test_that("plot_heatmaply works", {
  qf = add_log_imputed_norm_assays(make_minimal_conduit()@QFeatures)
  plot_heatmaply(qf, "species_log2_imputed_norm", "microbiome_treatment")
})
