test_that("maps EC numbers to enzyme names and trims trailing period", {
  dat <- c(
    "CC   -----------------------------------------------------------------------",
    "ID   1.1.1.1",
    "DE   Alcohol dehydrogenase.",
    "AN   Aldehyde reductase.",
    "//",
    "ID   1.1.1.2",
    "DE   Alcohol dehydrogenase (NADP(+)).",
    "//"
  )
  tmp <- tempfile(fileext = ".dat")
  writeLines(dat, tmp)
  on.exit(unlink(tmp))

  dict <- parse_enzyme_dat(tmp)

  expect_equal(dict$description[dict$term == "1.1.1.1"], "Alcohol dehydrogenase")
  expect_equal(dict$description[dict$term == "1.1.1.2"],
               "Alcohol dehydrogenase (NADP(+))")
  expect_equal(nrow(dict), 2L)
})
