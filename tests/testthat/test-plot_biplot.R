test_that("plot_biplot() works", {
  qf = add_log_imputed_norm_assays(make_minimal_conduit()@QFeatures)
  plot_biplot(qf, "species_log2_imputed_norm", color = "microbiome_treatment",
              facet_formula = ~day)
})
