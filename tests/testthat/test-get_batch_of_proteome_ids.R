test_that("this works", {
  skip_if_offline()
  expect_no_error(get_proteome_id_from_organism_id("818"))
})
