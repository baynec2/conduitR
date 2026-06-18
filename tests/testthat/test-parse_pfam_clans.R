test_that("maps Pfam accessions to their descriptions", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c(
    "PF00001\tCL0192\tGPCR_A\t7tm_1\t7 transmembrane receptor (rhodopsin family)",
    "PF00002\tCL0192\tGPCR_A\t7tm_2\t7 transmembrane receptor (Secretin family)"
  ), tmp)
  on.exit(unlink(tmp))

  dict <- parse_pfam_clans(tmp)

  expect_equal(dict$description[dict$term == "PF00001"],
               "7 transmembrane receptor (rhodopsin family)")
  expect_equal(dict$description[dict$term == "PF00002"],
               "7 transmembrane receptor (Secretin family)")
  expect_equal(nrow(dict), 2L)
})
