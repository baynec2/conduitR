test_that("perform_limma_analysis() works", {
  qf = add_log_imputed_norm_assays(make_minimal_conduit()@QFeatures)
  assay_name = "protein_group_log2_MinDet"
  formula = ~microbiome_treatment

  possible_contrasts = find_possible_contrast_terms(qf, assay_name, formula)

  out = perform_limma_analysis(qf, assay_name, formula, contrast = "microbiome_treatmenttreated")

  out
})
