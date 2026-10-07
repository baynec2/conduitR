qf <- readRDS(test_path("fixtures/add_annotation/qf.rds"))

conduit_annotations_wide <- readr::read_delim(test_path("fixtures/add_annotation/conduit_annotations.txt")) |>
  dplyr::select(Protein.Group, annotation_type, term) |>
  tidyr::pivot_wider(
    names_from = annotation_type,
    values_from = term,
    values_fn = \(x) paste(unique(x), collapse = ";")
  )

n_assays_before <- length(qf)

annotation_columns <- list(
  list(name = "go",                regex = "[^;]+(?=;|$)"),
  list(name = "pfam",              regex = "[^;]+(?=;|$)"),
  list(name = "eggnog",            regex = "[^;]+(?=;|$)"),
  list(name = "eggnog_code",       regex = "[^;]+(?=;)"),
  list(name = "kegg_pathway",      regex = "[^;]+(?=;|$)"),
  list(name = "kegg_map_pathway",  regex = "[^;]+(?=;|$)"),
  list(name = "kegg_orthology",    regex = "[^;]+(?=;|$)"),
  list(name = "cazy_class",        regex = "[^;]+(?=;|$)"),
  list(name = "cazy_family",       regex = "[^;]+(?=;|$)")
)

for (spec in annotation_columns) {
  local({
    nm    <- spec$name
    regex <- spec$regex
    test_that(paste("add_annotation_to_qf adds", nm, "as a rowData list-column and registers it"), {
      sym <- rlang::sym(nm)
      out <- rlang::eval_tidy(rlang::call2(
        "add_annotation_to_qf",
        qf                 = quote(qf),
        id_column          = quote(Protein.Group),
        column_name        = sym,
        conduit_annotations = quote(conduit_annotations_wide),
        regex              = regex,
        .ns                = "conduitR"
      ))

      expect_equal(length(out), n_assays_before)
      rd <- SummarizedExperiment::rowData(out[["protein_groups"]])
      expect_true(nm %in% colnames(rd))
      expect_true(is.list(rd[[nm]]))

      tgts <- conduitR::aggregation_targets(out)
      expect_true(nm %in% names(tgts))
      expect_equal(tgts[[nm]]$kind, "functional")
      expect_equal(tgts[[nm]]$from, "protein_groups")
    })
  })
}
