# Path helpers for test fixtures. Use these instead of hardcoding paths.
# test_path() resolves relative to tests/testthat/, which works under both
# devtools::test() and R CMD check.

fixture_path  <- function(...) test_path("fixtures", ...)
conduit_rds   <- function() test_path("fixtures/conduit.rds")
taxonomy_txt  <- function() test_path("fixtures/taxonomy.txt")
organism_ids_txt <- function() test_path("fixtures/organism_ids.txt")
detected_proteins_txt <- function() test_path("fixtures/detected_proteins.txt")
uniprot_ids_5000_rds  <- function() test_path("fixtures/uniprot_ids_5000.rds")
diann_parquet <- function() system.file("extdata/diann.parquet", package = "conduitR")
