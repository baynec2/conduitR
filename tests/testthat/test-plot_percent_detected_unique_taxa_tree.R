test_that("plot_percent_detected_taxa_tree() works", {
  skip_if_not(file.exists(conduit_rds()), "conduit.rds fixture not available")
  conduit_obj = readRDS(conduit_rds())
  has_metrics <- !inherits(
    tryCatch(slot(conduit_obj, "metrics"), error = identity),
    "error"
  )
  skip_if(!has_metrics, "conduit fixture predates the metrics slot — regenerate conduit.rds")
  expect_no_error(plot_percent_detected_taxa_tree(conduit_obj))
})
