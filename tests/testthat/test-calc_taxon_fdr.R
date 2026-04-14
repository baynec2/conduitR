# Shared test vectors: two families, each with one target and one decoy PSM.
# This mirrors the real DIA-NN data structure where decoy PSMs share the same
# family label as their source peptides.
pep   <- c(0.01,  0.05,  0.10,  0.50)
taxon <- c("A",   "A",   "B",   "B")
decoy <- c(FALSE, TRUE,  FALSE, TRUE)

# Hand-worked expected values:
#   A/target: score = -log(0.01) = 4.6052
#   A/decoy:  score = -log(0.05) = 2.9957
#   B/target: score = -log(0.10) = 2.3026
#   B/decoy:  score = -log(0.50) = 0.6931
#
# Ranked descending:
#   rank1: A/target  n_t=1 n_d=0 fdr=0.00  qvalue=0.00
#   rank2: A/decoy   n_t=1 n_d=1 fdr=1.00  qvalue=0.50
#   rank3: B/target  n_t=2 n_d=1 fdr=0.50  qvalue=0.50
#   rank4: B/decoy   n_t=2 n_d=2 fdr=1.00  qvalue=1.00

# --- input validation ---

test_that("errors on wrong pep type", {
  expect_error(calc_taxon_fdr("a", taxon, decoy), "`pep` must be a numeric vector")
})

test_that("errors on wrong taxon type", {
  expect_error(calc_taxon_fdr(pep, 1:4, decoy), "`taxon` must be a character vector")
})

test_that("errors on wrong decoy type", {
  expect_error(
    calc_taxon_fdr(pep, taxon, c(1L, 0L, 1L, 0L)),
    "`decoy` must be a logical vector"
  )
})

test_that("errors on fdr_threshold out of range or wrong type", {
  expect_error(calc_taxon_fdr(pep, taxon, decoy, fdr_threshold = 0),   "`fdr_threshold`")
  expect_error(calc_taxon_fdr(pep, taxon, decoy, fdr_threshold = 1),   "`fdr_threshold`")
  expect_error(calc_taxon_fdr(pep, taxon, decoy, fdr_threshold = "a"), "`fdr_threshold`")
})

test_that("errors on mismatched vector lengths", {
  expect_error(calc_taxon_fdr(c(0.1, 0.2), taxon, decoy), "same length")
})

test_that("errors on PEP outside [0, 1]", {
  expect_error(calc_taxon_fdr(c(-0.1, 0.05, 0.1, 0.5), taxon, decoy), "`pep` values")
  expect_error(calc_taxon_fdr(c(0.01, 0.05, 0.1, 1.1),  taxon, decoy), "`pep` values")
})

test_that("errors on NA in any input vector", {
  expect_error(calc_taxon_fdr(c(NA, 0.05, 0.1, 0.5), taxon, decoy),    "NA")
  expect_error(calc_taxon_fdr(pep, c(NA, "A", "B", "B"), decoy),        "NA")
  expect_error(calc_taxon_fdr(pep, taxon, c(NA, TRUE, FALSE, TRUE)),     "NA")
})

# --- correctness ---

test_that("aggregates target and decoy PSMs separately per taxon", {
  res <- calc_taxon_fdr(pep, taxon, decoy)$results
  expect_equal(res$score[res$taxon == "A" & !res$decoy], -log(0.01), tolerance = 1e-6)
  expect_equal(res$score[res$taxon == "A" &  res$decoy], -log(0.05), tolerance = 1e-6)
})

test_that("computes correct q-values from hand-worked example", {
  res <- calc_taxon_fdr(pep, taxon, decoy)$results
  expect_equal(res$qvalue[res$taxon == "A" & !res$decoy], 0.0, tolerance = 1e-6)
  expect_equal(res$qvalue[res$taxon == "A" &  res$decoy], 0.5, tolerance = 1e-6)
  expect_equal(res$qvalue[res$taxon == "B" & !res$decoy], 0.5, tolerance = 1e-6)
  expect_equal(res$qvalue[res$taxon == "B" &  res$decoy], 1.0, tolerance = 1e-6)
})

test_that("detected contains correct targets at fdr_threshold = 0.01", {
  res <- calc_taxon_fdr(pep, taxon, decoy, fdr_threshold = 0.01)
  expect_equal(nrow(res$detected), 1L)
  expect_equal(res$detected$taxon, "A")
})

test_that("detected contains correct targets at fdr_threshold = 0.6", {
  res <- calc_taxon_fdr(pep, taxon, decoy, fdr_threshold = 0.6)
  expect_equal(nrow(res$detected), 2L)
  expect_setequal(res$detected$taxon, c("A", "B"))
})

test_that("reports correct n_targets and n_decoys", {
  res <- calc_taxon_fdr(pep, taxon, decoy)
  expect_equal(res$n_targets, 2L)
  expect_equal(res$n_decoys,  2L)
})

# --- return structure ---

test_that("returns list with correct names and types", {
  res <- calc_taxon_fdr(pep, taxon, decoy)
  expect_named(res, c("results", "detected", "n_targets", "n_decoys", "fdr_threshold"))
  expect_s3_class(res$results, "tbl_df")
  expect_named(res$results, c("taxon", "score", "decoy", "fdr", "qvalue"))
  expect_equal(res$fdr_threshold, 0.01)
})

test_that("results are sorted descending by score", {
  expect_true(all(diff(calc_taxon_fdr(pep, taxon, decoy)$results$score) <= 0))
})

# --- edge cases ---

test_that("PEP of 0 is clamped to machine epsilon without error", {
  expect_no_error(
    res <- calc_taxon_fdr(c(0, 0.05, 0.1, 0.5), taxon, decoy)
  )
  expect_true(all(is.finite(res$results$score)))
})

test_that("peptide argument collapses to best PSM per (peptide, taxon, decoy)", {
  # Two target PSMs for peptide p1 / taxon A: PEP 0.1 and 0.5 → best is 0.1
  # One decoy PSM for peptide p1 / taxon A: PEP 0.2
  # One target PSM for peptide p2 / taxon B: PEP 0.05
  res <- calc_taxon_fdr(
    pep     = c(0.1,   0.5,   0.2,   0.05),
    taxon   = c("A",   "A",   "A",   "B"),
    decoy   = c(FALSE, FALSE, TRUE,  FALSE),
    peptide = c("p1",  "p1",  "p1",  "p2")
  )
  expect_equal(
    res$results$score[res$results$taxon == "A" & !res$results$decoy],
    -log(0.1), tolerance = 1e-6
  )
  expect_equal(
    res$results$score[res$results$taxon == "A" &  res$results$decoy],
    -log(0.2), tolerance = 1e-6
  )
})
