test_that("plot_upset() works", {
  qf <- make_minimal_conduit()@QFeatures
  expect_no_error(
    plot_upset(qf, "protein_group", group_by = "microbiome_treatment")
  )
})
