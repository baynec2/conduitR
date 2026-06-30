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

test_that("all-invalid IDs (e.g. UniParc) return an empty schema, not an error", {
  # Regression: an all-UniParc proteome strips every ID in validation, leaving
  # zero valid accessions. Previously this produced a 0-column tibble and the
  # downstream `mutate(go = ...)` failed with "object 'go' not found". No
  # network is hit (all IDs are stripped before any request), so this runs
  # offline. The function must return a well-formed empty tibble instead.
  out <- expect_no_error(
    suppressMessages(
      get_annotations_from_uniprot(c("UPI0000000053", "UPI0000000054", "not_an_id"),
                                   workers = 1)
    )
  )
  expect_s3_class(out, "tbl_df")
  expect_equal(nrow(out), 0L)
  expect_true(all(c("go", "xref_kegg", "xref_eggnog") %in% names(out)))
})
