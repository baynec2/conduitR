test_that("plot_sample_cor_heatmap() works", {
  qf = add_log_imputed_norm_assays(make_minimal_conduit()@QFeatures)

  expect_no_error(
    plot_sample_cor_heatmap(qf,
                            assay_name = "protein_group_log2_imputed_norm",
                            sample_annotation_variables = character(0))
  )
})
