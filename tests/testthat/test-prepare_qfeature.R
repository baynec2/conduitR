test_that("prepare_qfeature works as intended", {
  sample_annotation_fp <- "tests/data/prepare_qfeature/sample_annotation.txt"

  vector_of_matrix_fps <- c(
    "tests/data/prepare_qfeature/superkingdom_matrix.tsv",
    "tests/data/prepare_qfeature/subcellular_locations_matrix.tsv",
    "tests/data/prepare_qfeature/kingdom_matrix.tsv",
    "tests/data/prepare_qfeature/phylum_matrix.tsv",
    "tests/data/prepare_qfeature/order_matrix.tsv",
    "tests/data/prepare_qfeature/class_matrix.tsv",
    "tests/data/prepare_qfeature/family_matrix.tsv",
    "tests/data/prepare_qfeature/genus_matrix.tsv",
    "tests/data/prepare_qfeature/species_matrix.tsv",
    "tests/data/prepare_qfeature/precursor_matrix.tsv",
    "tests/data/prepare_qfeature/peptide_matrix.tsv",
    "tests/data/prepare_qfeature/protein_group_matrix.tsv",
    "tests/data/prepare_qfeature/go_matrix.tsv",
    "tests/data/prepare_qfeature/go_taxa_matrix.tsv"
  )


  expect_no_error(prepare_qfeature(
    sample_annotation_fp,
    vector_of_matrix_fps
  ))
})
