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
  # Note this only extracts the first kegg term.
  expect_no_error(get_kegg_in_batches(kegg_ids))
  })


test_that("works with problem ;", {
  kegg_ids = "rbc:BN938_0910;"
  expect_no_error(get_kegg_in_batches(kegg_ids))
})


# Having trouble getting it to run with multiple batches.
test_that("works with multiple batches",{
  kegg_ids <- c(
    "hsa:1636;", "hsa:999; hsa:4625", "hsa:6714", "hsa:5184", "hsa:1327", "hsa:4624", "hsa:1938",
    "hsa:3858", "hsa:3860", "hsa:3852", "hsa:9601", "hsa:54", "hsa:210", "hsa:3936", "hsa:5358",
    "hsa:327", "hsa:2108", "hsa:5576", "hsa:6523", "hsa:653145", "hsa:728113", "hsa:2027", "hsa:2224",
    "hsa:1368", "hsa:6476", "hsa:10327", "hsa:5315", "hsa:7184", "hsa:6628", "hsa:1340", "hsa:3728",
    "hsa:7381", "hsa:1357", "hsa:1360", "hsa:290", "hsa:5224", "hsa:410", "hsa:2683", "hsa:7430",
    "hsa:4830", "hsa:2799", "hsa:6187", "hsa:1832", "hsa:1350", "hsa:3014", "hsa:873", "hsa:35",
    "hsa:5406", "hsa:2720", "hsa:3009"
  )

expect_no_error(get_kegg_in_batches(kegg_ids))

})

# Having a problem with a specific term

test_that("works with problematic term",{
          problem = "fbs:N510_000502"
          get_kegg_in_batches(problem)



# Having trouble getting it to run with multiple batches.
test_that("pathways and pathway ids match",{
  kegg_ids <- c(
    "hsa:1636;"
  )

t3 = get_kegg_in_batches(kegg_ids)

})
