test_that("uniprot api works", {
  skip_if_offline()
  expect_no_error({
    check_api_service("https://rest.uniprot.org/proteomes/UP000005640")})

})

test_that("uniprot api works", {
  skip_if_offline()
  expect_error({
    check_api_service("bad")})
})

test_that("ncbi taxonomy api works ", {
  skip_if_offline()
  expect_no_error({
    check_api_service("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=taxonomy")
  })
})
