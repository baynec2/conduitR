test_that("plot_kegg_pathway works", {
  skip_if_not_installed("ggkegg")
  skip_if_offline()
  stats <- tibble::tibble(xref_kegg = list("K00001"), logFC = 1.5)
  expect_no_error(plot_kegg_pathway(stats_results = stats, kegg_pathway_id = "ko00010"))
})
