test_that("plot_percent_detected_taxa_tree() works", {
  conduit_obj = make_minimal_conduit()
  expect_no_error(plot_percent_detected_taxa_tree(conduit_obj))
})
