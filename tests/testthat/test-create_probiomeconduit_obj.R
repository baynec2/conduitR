test_that("can create a conduit_obj", {
 expect_no_error(create_conduitt_obj("tests/data/prepare_qfeature/qf.rds",
                           database_taxonomy.tsv = "tests/data/prepare_qfeature/database_taxonomy.tsv",
                           database_metrics.tsv = "tests/data/prepare_qfeature/database_metrics.tsv",
                           detected_protein_taxonomy.tsv = "tests/data/prepare_qfeature/detected_protein_taxonomy.tsv",
                           detected_protein_metrics.tsv = "tests/data/prepare_qfeature/detected_protein_metrics.tsv"))
})
