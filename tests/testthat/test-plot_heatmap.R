test_that("plot_heatmap() works", {
  qf = readRDS("tests/data/conduit.rds")@QFeatures
  plot_heatmap(qf,
               "species",
               c("microbiome_treatment","day","tumor_classification"),
               "organism_type")
})
