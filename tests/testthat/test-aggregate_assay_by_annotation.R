make_listcol_test_qf <- function() {
  feats <- paste0("PG", 1:5)
  samps <- c("S1", "S2", "S3")
  mat <- matrix(
    c( 10, 100, 1000,
       20, 200, 2000,
       30, 300, 3000,
       40, 400, 4000,
       50, 500, 5000),
    nrow = 5, byrow = TRUE,
    dimnames = list(feats, samps)
  )
  go <- list(
    PG1 = c("GO:1", "GO:2"),
    PG2 = c("GO:1"),
    PG3 = c("GO:3"),
    PG4 = NA_character_,
    PG5 = character(0)
  )
  pg_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = mat),
    rowData = S4Vectors::DataFrame(row.names = feats, go = I(go))
  )
  QFeatures::QFeatures(list(protein_groups = pg_se))
}

make_atomic_test_qf <- function() {
  feats <- paste0("pep", 1:6)
  samps <- c("S1", "S2", "S3")
  mat <- matrix(
    c( 1,  10,  100,
       2,  20,  200,
       3,  30,  300,
       4,  40,  400,
       5,  50,  500,
       6,  60,  600),
    nrow = 6, byrow = TRUE,
    dimnames = list(feats, samps)
  )
  species <- c("E_coli", "E_coli", "S_aureus", "S_aureus", NA, NA)
  pep_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(intensity = mat),
    rowData = S4Vectors::DataFrame(row.names = feats, species = species)
  )
  QFeatures::QFeatures(list(peptides = pep_se))
}

# list-col, include_na = "drop" --------------------------------------------

test_that("listcol drop: features without terms contribute to no term row", {
  qf  <- make_listcol_test_qf()
  out <- conduitR::aggregate_assay_by_annotation(
    qf, i = "protein_groups", fcol = "go", include_na = "drop"
  )

  expect_true("go" %in% names(out))
  ag <- SummarizedExperiment::assay(out[["go"]])
  expect_setequal(rownames(ag), c("GO:1", "GO:2", "GO:3"))

  # GO:1 = PG1+PG2 ; GO:2 = PG1 ; GO:3 = PG3
  expect_equal(unname(ag["GO:1", ]), c(10+20, 100+200, 1000+2000))
  expect_equal(unname(ag["GO:2", ]), c(10, 100, 1000))
  expect_equal(unname(ag["GO:3", ]), c(30, 300, 3000))
})

# list-col, include_na = "group" -------------------------------------------

test_that("listcol group: features without terms aggregate into Unassigned bucket", {
  qf  <- make_listcol_test_qf()
  out <- conduitR::aggregate_assay_by_annotation(
    qf, i = "protein_groups", fcol = "go", include_na = "group"
  )

  expect_true("Unassigned" %in% rownames(SummarizedExperiment::assay(out[["go"]])))
  ag <- SummarizedExperiment::assay(out[["go"]])

  base <- SummarizedExperiment::assay(qf[["protein_groups"]])
  unmapped <- c("PG4", "PG5")
  expect_equal(unname(ag["Unassigned", ]),
               unname(colSums(base[unmapped, , drop = FALSE])))
})

# atomic, include_na = "drop" ----------------------------------------------

test_that("atomic drop: NA-annotated features dropped, column sums shrink", {
  qf  <- make_atomic_test_qf()
  out <- conduitR::aggregate_assay_by_annotation(
    qf, i = "peptides", fcol = "species", include_na = "drop"
  )

  expect_setequal(rownames(SummarizedExperiment::assay(out[["species"]])),
                  c("E_coli", "S_aureus"))
  ag   <- SummarizedExperiment::assay(out[["species"]])
  base <- SummarizedExperiment::assay(qf[["peptides"]])
  # NA-species peptides (pep5, pep6) dropped
  expect_lt(sum(ag[, "S1"]), sum(base[, "S1"]))
  expect_equal(unname(colSums(ag)),
               unname(colSums(base[c("pep1","pep2","pep3","pep4"), ])))
})

# atomic, include_na = "group" ---------------------------------------------

test_that("atomic group: NA-annotated features bucketed; column sums preserved", {
  qf  <- make_atomic_test_qf()
  out <- conduitR::aggregate_assay_by_annotation(
    qf, i = "peptides", fcol = "species", include_na = "group"
  )

  ag <- SummarizedExperiment::assay(out[["species"]])
  expect_setequal(rownames(ag), c("E_coli", "S_aureus", "Unassigned"))

  base <- SummarizedExperiment::assay(qf[["peptides"]])
  expect_equal(unname(colSums(ag)), unname(colSums(base)))
  expect_equal(unname(ag["Unassigned", ]),
               unname(colSums(base[c("pep5","pep6"), ])))

  # Base assay rowData should NOT carry the temp ".species_with_unassigned"
  # column after aggregation completes
  base_rd <- SummarizedExperiment::rowData(out[["peptides"]])
  expect_false(".species_with_unassigned" %in% colnames(base_rd))
})

test_that("missing assay or column produces a clear error", {
  qf <- make_atomic_test_qf()
  expect_error(
    conduitR::aggregate_assay_by_annotation(qf, i = "nope", fcol = "species"),
    "not found"
  )
  expect_error(
    conduitR::aggregate_assay_by_annotation(qf, i = "peptides", fcol = "nope"),
    "not found"
  )
})
