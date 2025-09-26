test_that("this works", {
  pg_groups = readr::read_tsv("inst/extdata/report.pg_matrix.tsv")

  uniprot_ids = pg_groups |>
    dplyr::pull("Protein.Group") |>
    strsplit(";") |>
    unlist() |>
    unique()

  out = get_annotations_from_uniprot(uniprot_ids)
})

test_that("this works with ids that were previously having col type problems",{
  expect_no_error({
  problem_id = "A0A0R4J083"
  conduitR::get_annotations_from_uniprot(problem_id)

  })

})

test_that("this works with ids that were problematic previously",{
  expect_no_error({
    problems = readr::read_delim("tests/data/detected_protein_info.txt") |>
      pull(protein_id)


    out = conduitR::get_annotations_from_uniprot(problems)

  })

})


