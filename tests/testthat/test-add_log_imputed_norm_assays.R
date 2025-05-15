test_that("add_log_imputed_norm_assays() works", {
  qf = readRDS("tests/data/conduit.rds")@QFeatures
  qf_final = add_log_imputed_norm_assays(qf)
    })
