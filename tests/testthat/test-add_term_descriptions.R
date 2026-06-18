test_that("fills blank descriptions from the dictionary by (annotation_type, term)", {
  annotations <- tibble::tibble(
    protein_id      = c("P1", "P2", "P3"),
    annotation_type = c("go", "go", "kegg_orthology"),
    term            = c("GO:0008150", "GO:0003674", "K00001"),
    description     = c(NA_character_, "", NA_character_)
  )
  dictionary <- tibble::tibble(
    annotation_type = c("go", "go", "kegg_orthology"),
    term            = c("GO:0008150", "GO:0003674", "K00001"),
    description     = c("biological_process", "molecular_function", "alcohol dehydrogenase")
  )

  out <- add_term_descriptions(annotations, dictionary)

  expect_equal(out$description,
               c("biological_process", "molecular_function", "alcohol dehydrogenase"))
  # Row count and column set preserved.
  expect_equal(nrow(out), nrow(annotations))
  expect_setequal(names(out), names(annotations))
})

test_that("never overwrites a description that is already present", {
  annotations <- tibble::tibble(
    annotation_type = "go",
    term            = "GO:0008150",
    description     = "existing source description"
  )
  dictionary <- tibble::tibble(
    annotation_type = "go",
    term            = "GO:0008150",
    description     = "dictionary description"
  )

  out <- add_term_descriptions(annotations, dictionary)
  expect_equal(out$description, "existing source description")
})

test_that("does not borrow descriptions across annotation_type (no in-object fallback)", {
  # uniprot_go has a name, eggNOG 'go' does not; they must stay independent.
  annotations <- tibble::tibble(
    annotation_type = c("uniprot_go", "go"),
    term            = c("GO:0008150", "GO:0008150"),
    description     = c("biological_process", NA_character_)
  )
  # Dictionary only describes 'go'; uniprot_go keeps its own.
  dictionary <- tibble::tibble(
    annotation_type = "go",
    term            = "GO:0008150",
    description     = "biological_process"
  )

  out <- add_term_descriptions(annotations, dictionary)
  expect_equal(out$description, c("biological_process", "biological_process"))

  # With no dictionary entry for 'go', it must remain NA rather than borrowing
  # the uniprot_go name.
  out2 <- add_term_descriptions(
    annotations,
    tibble::tibble(annotation_type = character(),
                   term = character(),
                   description = character())
  )
  expect_equal(out2$description, c("biological_process", NA_character_))
})

test_that("adds a description column when annotations lacks one", {
  annotations <- tibble::tibble(
    annotation_type = "go",
    term            = "GO:0008150"
  )
  dictionary <- tibble::tibble(
    annotation_type = "go",
    term            = "GO:0008150",
    description     = "biological_process"
  )
  out <- add_term_descriptions(annotations, dictionary)
  expect_true("description" %in% names(out))
  expect_equal(out$description, "biological_process")
})

test_that("unmatched terms keep NA and rows are not duplicated", {
  annotations <- tibble::tibble(
    annotation_type = c("go", "go"),
    term            = c("GO:0000001", "GO:0008150"),
    description     = NA_character_
  )
  # Duplicate dictionary rows for one term should not multiply output rows.
  dictionary <- tibble::tibble(
    annotation_type = c("go", "go"),
    term            = c("GO:0008150", "GO:0008150"),
    description     = c("biological_process", "biological_process (dup)")
  )
  out <- add_term_descriptions(annotations, dictionary)
  expect_equal(nrow(out), 2L)
  expect_equal(out$description, c(NA_character_, "biological_process"))
})

test_that("errors when required columns are missing", {
  expect_error(
    add_term_descriptions(tibble::tibble(term = "GO:1"),
                          tibble::tibble(annotation_type = "go", term = "GO:1",
                                         description = "x")),
    "annotation_type"
  )
  expect_error(
    add_term_descriptions(tibble::tibble(annotation_type = "go", term = "GO:1"),
                          tibble::tibble(annotation_type = "go", term = "GO:1")),
    "description"
  )
})
