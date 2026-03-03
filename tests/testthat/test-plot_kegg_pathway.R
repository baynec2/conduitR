test_that("plot_kegg_pathway works", {
  skip_if_not_installed("ggkegg")
  skip_if_offline()
  expect_no_error(plot_kegg_pathway())
})
