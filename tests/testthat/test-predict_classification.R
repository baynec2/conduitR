qf = conduitR::add_log_imputed_norm_assays(readRDS("tests/data/conduit.rds")@QFeatures)


test_that("predict classification lasso regression works", {
  assay_name = "protein_group_log2_imputed"
  outcome = "microbiome_treatment"
  train_percent = 70
  model_type = "lasso_regression"
  v=1
  grid_size = 20

  predict_classification(qf,assay_name = "protein_group_log2_imputed",outcome = "microbiome_treatment",train_percent = 70,model_type = "lasso_regression", v=1, grid_size = 20)




})
