test_that("add_log_imputed_norm_assays() works", {
  qf = make_minimal_conduit()@QFeatures
  qf_final = add_log_imputed_norm_assays(qf)
})
