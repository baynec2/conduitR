qf <- readRDS("inst/extdata/add_annotation/qf.rds")

conduit_annotations_wide <- readr::read_delim("inst/extdata/add_annotation/conduit_annotations.txt") |>
  dplyr::select(Protein.Group, annotation_type, term) |>  # keep columns of interest
  # pivot so each annotation_type becomes a column
  tidyr::pivot_wider(
    names_from = annotation_type,
    values_from = term,
    values_fn = \(x) paste(unique(x), collapse = ";")  # collapse multiple terms per protein
  )

test_that("it works with go", {
  expect_no_error({
  go <- conduitR::add_annotation_to_qf(qf,
                             id_column = Protein.Group,
                             column_name = go,
                             conduit_annotations = conduit_annotations_wide)
   })
})


test_that("it works with pfam", {
  expect_no_error({
    pfam <- conduitR::add_annotation_to_qf(qf,
                                         id_column = Protein.Group,
                                         column_name = pfam,
                                         conduit_annotations = conduit_annotations_wide)
  })
})

test_that("it works with eggnog", {
  expect_no_error({
    eggnog <- conduitR::add_annotation_to_qf(qf,
                                           id_column = Protein.Group,
                                           column_name = eggnog,
                                           conduit_annotations = conduit_annotations_wide)
  })
})

test_that("it works with eggnog code", {
  expect_no_error({
    eggnog_code <- conduitR::add_annotation_to_qf(qf,
                                           id_column = Protein.Group,
                                           column_name = eggnog_code,
                                           conduit_annotations = conduit_annotations_wide,
                                           regex = "[^;]+(?=;)")
  })
})


test_that("it works with kegg_pathway", {
  expect_no_error({
    test <- conduitR::add_annotation_to_qf(qf,
                                           id_column = Protein.Group,
                                           column_name = kegg_pathway,
                                           conduit_annotations = conduit_annotations_wide)
  })
})

test_that("it works with kegg_map_pathway", {
  expect_no_error({
    test <- conduitR::add_annotation_to_qf(qf,
                                           id_column = Protein.Group,
                                           column_name = kegg_map_pathway,
                                           conduit_annotations = conduit_annotations_wide
                                           )
  })
})

test_that("it works with kegg_orthology", {
  expect_no_error({
    test <- conduitR::add_annotation_to_qf(qf,
                                           id_column = Protein.Group,
                                           column_name = kegg_orthology,
                                           conduit_annotations = conduit_annotations_wide
                                           )
  })
})


test_that("it works with cazy class", {
  expect_no_error({
    test <- conduitR::add_annotation_to_qf(qf,
                                           id_column = Protein.Group,
                                           column_name = cazy_class,
                                           conduit_annotations = conduit_annotations_wide)
  })
})


test_that("it works with cazy family", {
  expect_no_error({
    test <- conduitR::add_annotation_to_qf(qf,
                                           id_column = Protein.Group,
                                           column_name = cazy_family,
                                           conduit_annotations = conduit_annotations_wide)
  })
})

