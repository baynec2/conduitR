test_that("add_relative_abundance_assays() works", {
  qf = make_minimal_conduit()@QFeatures
  qf_final = add_log_imputed_norm_assays(qf)
  rel_abundance = add_relative_abundance_assays(qf_final, "species")
})
