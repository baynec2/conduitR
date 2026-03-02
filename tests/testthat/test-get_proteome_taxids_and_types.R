test_that("this works with UP001369780, a uniprot id that is an edge case and
          not actually in the database", {
  skip_if_offline()
  expect_no_error(get_proteome_taxids_and_types("UP001369780"))
})

test_that("this works with UP000005640, the human proteome id", {
  skip_if_offline()
            expect_no_error(get_proteome_taxids_and_types("UP000005640"))
          })

test_that("this works with bad id", {
  skip_if_offline()
  expect_no_error(get_proteome_taxids_and_types("badid"))
})


test_that("this works when bad ids are provided in combination with good ids", {
  skip_if_offline()
  expect_no_error(get_proteome_taxids_and_types(c("badid","UP000005640")))
})
