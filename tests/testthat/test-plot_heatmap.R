test_that("plot_heatmap() works", {
  qf = add_log_imputed_norm_assays(make_minimal_conduit()@QFeatures)
  expect_no_error(
    plot_heatmap(qf,
                 "species_log2_imputed_norm",
                 c("microbiome_treatment", "day", "tumor_classification"),
                 "species",
                 scale = FALSE)
  )
})
