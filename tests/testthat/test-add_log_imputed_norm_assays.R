test_that("add_log_imputed_norm_assays() works", {
  qf = readRDS(conduit_rds())@QFeatures
  qf_final = add_log_imputed_norm_assays(qf)
    })
