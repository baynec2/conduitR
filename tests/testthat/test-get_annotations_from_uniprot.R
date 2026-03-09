test_that("this works", {
  skip_if_offline()
  pg_groups = readr::read_tsv(system.file("extdata/report.pg_matrix.tsv", package = "conduitR"))

  uniprot_ids = pg_groups |>
    dplyr::pull("Protein.Group") |>
    strsplit(";") |>
    unlist() |>
    unique() |>
    head(20)

  out = get_annotations_from_uniprot(uniprot_ids)
})

test_that("this works with ids that were previously having col type problems",{
  skip_if_offline()
  expect_no_error({
  problem_id = "A0A0R4J083"
  conduitR::get_annotations_from_uniprot(problem_id)

  })

})

test_that("this works with ids that were problematic previously",{
  skip_if_offline()
  expect_no_error({
    problems = readr::read_delim(detected_proteins_txt()) |>
      dplyr::pull(protein_id) |>
      head(20)

    out = conduitR::get_annotations_from_uniprot(problems)

  })

})
