test_that("predict classification lasso regression works", {
  skip_if_not(file.exists(conduit_rds()), "conduit.rds fixture not available")
  skip_if_not_installed("rsample")
  qf = conduitR::add_log_imputed_norm_assays(readRDS(conduit_rds())@QFeatures)
  assay_name = "protein_group_log2_imputed"
  outcome = "microbiome_treatment"
  train_percent = 70
  model_type = "lasso_regression"
  v=1
  grid_size = 20

  expect_no_error(
    predict_classification(qf, assay_name = "protein_group_log2_imputed", outcome = "microbiome_treatment", train_percent = 70, model_type = "lasso_regression", v = 1, grid_size = 20)
  )
})
