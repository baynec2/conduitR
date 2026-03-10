test_that("find_possible_contrasts() works", {
  qf = add_log_imputed_norm_assays(make_minimal_conduit()@QFeatures)
  find_possible_contrast_terms(qf, "protein_group_log2_imputed",
                               ~microbiome_treatment + immunotherapy_treatment)
})
