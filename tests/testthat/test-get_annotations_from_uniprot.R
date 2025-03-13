test_that("this works", {
  pg_groups = readr::read_tsv("inst/extdata/report.pg_matrix.tsv")

  uniprot_ids = pg_groups |>
    dplyr::pull("Protein.Group") |>
    strsplit(";") |>
    unlist() |>
    unique()

  out = get_annotations_from_uniprot(uniprot_ids)
})


