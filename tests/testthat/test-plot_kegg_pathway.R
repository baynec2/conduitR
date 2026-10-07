test_that("plot_kegg_pathway works", {
  skip_if_not_installed("ggkegg")
  skip_if_offline()
  # "ko00010" is a KO pathway, for which kegg_col defaults to "kegg_orthology"
  # -- and K-numbers are KO identifiers. The fixture previously supplied them
  # under xref_kegg, which is the column used for organism-specific pathways.
  stats <- tibble::tibble(kegg_orthology = list("K00001"), logFC = 1.5)
  expect_no_error(plot_kegg_pathway(stats_results = stats, kegg_pathway_id = "ko00010"))
})
