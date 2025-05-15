test_that("perform_limma_analysis() works", {
  qf = add_log_imputed_norm_assays(readRDS("tests/data/conduit.rds")@QFeatures)
  assay_name = "protein_group_log2_imputed"
  formula = ~microbiome_treatment

  possible_contrasts = find_possible_contrast_terms(qf,assay_name,formula)


  out = perform_limma_analysis(qf,assay_name,formula,contrast = "microbiome_treatmentnone")

  out
})
