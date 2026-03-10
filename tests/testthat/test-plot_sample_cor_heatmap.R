test_that("plot_sample_cor_heatmap() works", {
  skip_if_not(file.exists(conduit_rds()), "conduit.rds fixture not available")
  qf = readRDS(conduit_rds())@QFeatures
  qf = add_log_imputed_norm_assays(qf)

  expect_no_error(
    plot_sample_cor_heatmap(qf,
                            assay_name = "protein_group_log2_imputed_norm",
                            sample_annotation_variables = character(0))
  )
})
