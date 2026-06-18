test_that("parses primary ids, names and alt_ids from an OBO file", {
  obo <- c(
    "format-version: 1.2",
    "",
    "[Term]",
    "id: GO:0008150",
    "name: biological_process",
    "namespace: biological_process",
    "alt_id: GO:0000004",
    "alt_id: GO:0007582",
    "",
    "[Term]",
    "id: GO:0003674",
    "name: molecular_function",
    "",
    "[Typedef]",
    "id: part_of",
    "name: part of"
  )
  tmp <- tempfile(fileext = ".obo")
  writeLines(obo, tmp)
  on.exit(unlink(tmp))

  dict <- parse_go_obo(tmp)

  # Primary terms present with their names.
  expect_equal(dict$description[dict$term == "GO:0008150"], "biological_process")
  expect_equal(dict$description[dict$term == "GO:0003674"], "molecular_function")

  # alt_ids resolve to the primary term's name.
  expect_equal(dict$description[dict$term == "GO:0000004"], "biological_process")
  expect_equal(dict$description[dict$term == "GO:0007582"], "biological_process")

  # Typedef stanzas are ignored (part_of must not appear as a GO term).
  expect_false("part_of" %in% dict$term)
})
