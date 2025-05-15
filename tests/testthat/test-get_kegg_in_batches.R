test_that("works with no trailing ;", {
  kegg_ids = "ype:YPO2670"
  expect_no_error(get_kegg_in_batches(kegg_ids))
})


test_that("works with trailing ;", {
  kegg_ids = "ype:YPO2670;"
  expect_no_error(get_kegg_in_batches(kegg_ids))

})


test_that("works with multiple ;", {
  kegg_ids = "ype:YPO2670;ypk:y1242;ypm:YP_2471;"
  expect_no_error(get_kegg_in_batches(kegg_ids))

  })


test_that("works with problem ;", {
  kegg_ids = "rbc:BN938_0910;"
  expect_no_error(get_kegg_in_batches(kegg_ids))

})

