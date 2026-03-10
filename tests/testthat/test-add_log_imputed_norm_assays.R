test_that("add_log_imputed_norm_assays() works", {
  skip_if_not(file.exists(conduit_rds()), "conduit.rds fixture not available")
  qf = readRDS(conduit_rds())@QFeatures
  qf_final = add_log_imputed_norm_assays(qf)
    })
