test_that("find_possible_contrasts() works", {
  qf = add_log_imputed_norm_assays(readRDS("tests/data/conduit.rds")@QFeatures)
  find_possible_contrast_terms(qf,"protein_group_log2_imputed",
                               ~microbiome_treatment + immunotherapy_treatment)
})
