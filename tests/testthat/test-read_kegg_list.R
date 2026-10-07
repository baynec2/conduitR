test_that("strips db prefixes for ko, pathway, module and brite lists", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c(
    "ko:K00001\tE1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]",
    "path:map00010\tGlycolysis / Gluconeogenesis",
    "md:M00001\tGlycolysis (Embden-Meyerhof pathway)",
    "br:ko00000\tKEGG Orthology (KO)"
  ), tmp)
  on.exit(unlink(tmp))

  dict <- read_kegg_list(tmp)

  expect_equal(dict$description[dict$term == "K00001"],
               "E1.1.1.1, adh; alcohol dehydrogenase [EC:1.1.1.1]")
  expect_equal(dict$description[dict$term == "map00010"],
               "Glycolysis / Gluconeogenesis")
  expect_equal(dict$description[dict$term == "M00001"],
               "Glycolysis (Embden-Meyerhof pathway)")
  expect_equal(dict$description[dict$term == "ko00000"], "KEGG Orthology (KO)")
})
