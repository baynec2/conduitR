test_that("find_possible_contrasts() works", {
  skip_if_not(file.exists(conduit_rds()), "conduit.rds fixture not available")
  qf = add_log_imputed_norm_assays(readRDS(conduit_rds())@QFeatures)
  find_possible_contrast_terms(qf,"protein_group_log2_imputed",
                               ~microbiome_treatment + immunotherapy_treatment)
})
