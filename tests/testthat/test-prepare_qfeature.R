test_that("prepare_qfeature works as intended", {
  sample_annotation_fp <- test_path("fixtures/prepare_qfeature/sample_annotation.txt")
  skip_if_not(file.exists(sample_annotation_fp), "prepare_qfeature fixtures missing — generate them first")

  vector_of_matrix_fps <- c(
    test_path("fixtures/prepare_qfeature/superkingdom_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/subcellular_locations_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/kingdom_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/phylum_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/order_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/class_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/family_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/genus_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/species_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/precursor_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/peptide_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/protein_group_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/go_matrix.tsv"),
    test_path("fixtures/prepare_qfeature/go_taxa_matrix.tsv")
  )


  expect_no_error(prepare_qfeature(
    sample_annotation_fp,
    vector_of_matrix_fps
  ))
})
