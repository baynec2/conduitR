test_that("this works with UP001369780, a uniprot id that is an edge case and
          not actually in the database", {
  expect_no_error(get_proteome_taxids_and_types("UP001369780"))
})

test_that("this works with UP000005640, the human proteome id", {
            expect_no_error(get_proteome_taxids_and_types("UP000005640"))
          })

test_that("this works with bad id", {
  expect_no_error(get_proteome_taxids_and_types("badid"))
})


test_that("this works when bad ids are provided in combination with good ids", {
  expect_no_error(get_proteome_taxids_and_types(c("badid","UP000005640")))
})
