test_that("predict classification lasso regression works", {
  skip_if_not_installed("rsample")
  skip_if_not_installed("glmnet")
  qf = conduitR::add_log_imputed_norm_assays(make_minimal_conduit()@QFeatures)

  expect_no_error(
    predict_classification(qf, assay_name = "protein_group_log2_MinDet", outcome = "microbiome_treatment", train_percent = 70, model_type = "lasso_regression", v = 1, grid_size = 20)
  )
})
