# =============================================================================
# Input validation — calc_taxon_fdr (PSM-level entry point)
# =============================================================================
pep     <- c(0.01,  0.80,  0.50,  0.02)
taxon   <- c("A",   "A",   "B",   "B")
decoy   <- c(FALSE, TRUE,  FALSE, TRUE)
peptide <- c("PEPA1", "PEPA2", "PEPB1", "PEPB2")

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

test_that("errors on qvalue_threshold / min_peptides out of range", {
  expect_error(calc_taxon_fdr(pep, taxon, decoy, qvalue_threshold = 0),  "qvalue_threshold")
  expect_error(calc_taxon_fdr(pep, taxon, decoy, qvalue_threshold = 1),  "qvalue_threshold")
  expect_error(calc_taxon_fdr(pep, taxon, decoy, qvalue_threshold = "a"), "qvalue_threshold")
  expect_error(calc_taxon_fdr(pep, taxon, decoy, min_peptides = -1),     "min_peptides")
})

test_that("errors on mismatched lengths and NA inputs", {
  expect_error(calc_taxon_fdr(pep[-1], taxon, decoy), "same length")
  expect_error(calc_taxon_fdr(c(NA, 0.1, 0.1, 0.1), taxon, decoy), "must not contain NA")
})

test_that("errors on out-of-range pep", {
  expect_error(calc_taxon_fdr(c(1.5, 0.1, 0.1, 0.1), taxon, decoy), "between 0 and 1")
})

# =============================================================================
# calc_taxon_fdr — picked competition end to end (PSM level)
# =============================================================================
# A: target PEP 0.01 (strong) vs decoy 0.80 (weak) -> target wins the pick
# B: target PEP 0.50 (weak)   vs decoy 0.02 (strong) -> decoy wins the pick
test_that("picks the higher-scoring side of each taxon and aggregates score", {
  out  <- calc_taxon_fdr(pep, taxon, decoy, peptide)
  reps <- out$results[out$results$picked_winner, ]
  expect_equal(nrow(reps), 2L)
  expect_false(reps$decoy[reps$taxon == "A"])  # A target won
  expect_true(reps$decoy[reps$taxon == "B"])   # B decoy won
  # score is sum(-log(PEP)); A/target = -log(0.01)
  expect_equal(reps$score[reps$taxon == "A"], -log(0.01))
})

test_that("discarded losers carry NA fdr/qvalue and pass = FALSE", {
  out    <- calc_taxon_fdr(pep, taxon, decoy, peptide)
  losers <- out$results[!out$results$picked_winner, ]
  expect_equal(nrow(losers), 2L)
  expect_true(all(is.na(losers$fdr)))
  expect_true(all(is.na(losers$qvalue)))
  expect_true(all(!losers$pass))
})

test_that("min_peptides gate is skipped when peptide is not supplied", {
  out <- calc_taxon_fdr(pep, taxon, decoy)        # no peptide
  expect_true(all(is.na(out$results$n_unique_peptides_all)))
  expect_type(out$results$pass, "logical")
})

test_that("reports rep counts, first_decoy_rank and pairing", {
  out <- calc_taxon_fdr(pep, taxon, decoy, peptide)
  expect_equal(out$n_targets, 1L)         # only A/target is a target rep
  expect_equal(out$n_decoys, 1L)          # only B/decoy is a decoy rep
  expect_equal(out$n_missing_pair, 0L)    # both taxa have target and decoy
  expect_true(is.numeric(out$first_decoy_rank))
})

# =============================================================================
# Internal picked competition — pick_taxon_fdr_compete (aggregated table)
# =============================================================================
# Hand-worked: two taxa already aggregated.
#   A: target score 100 (40 pep), decoy 3   -> target wins
#   B: target score 2   (1 pep),  decoy 5   -> decoy wins
# Representatives, score-descending:
#   rank1 A/target  cum_d=0 cum_t=1 fdr=(0+1)/1=1.0
#   rank2 B/decoy   cum_d=1 cum_t=1 fdr=(1+1)/1=2.0
agg_taxon <- c("A", "A", "B", "B")
agg_score <- c(100,  3,   2,   5)
agg_decoy <- c(FALSE, TRUE, FALSE, TRUE)
agg_npep  <- c(40,   2,   1,   3)

test_that("compete uses the +1 guard and a monotone q-value", {
  out  <- conduitR:::pick_taxon_fdr_compete(agg_taxon, agg_score, agg_decoy, agg_npep)
  reps <- out$results[out$results$picked_winner, ]
  reps <- reps[order(-reps$score), ]
  expect_equal(reps$fdr, c(1.0, 2.0))
  # q-value is non-increasing as score increases -> non-decreasing here.
  expect_false(is.unsorted(reps$qvalue))
})

test_that("compete pass requires target winner AND qvalue AND min_peptides", {
  tx <- c("A", "A", "noise", "noise")
  sc <- c(100, 1, 0.5, 0.4)
  dc <- c(FALSE, TRUE, FALSE, TRUE)
  np <- c(40, 1, 1, 1)
  a  <- conduitR:::pick_taxon_fdr_compete(tx, sc, dc, np, qvalue_threshold = 0.5, min_peptides = 2)
  expect_true(a$results$pass[a$results$taxon == "A" & a$results$picked_winner])
  a2 <- conduitR:::pick_taxon_fdr_compete(tx, sc, dc, np, qvalue_threshold = 0.5, min_peptides = 99)
  expect_false(a2$results$pass[a2$results$taxon == "A" & a2$results$picked_winner])
})

test_that("compete flags taxa present on only one side", {
  out <- conduitR:::pick_taxon_fdr_compete(
    taxon = c("A", "A", "B"), score = c(10, 2, 5),
    decoy = c(FALSE, TRUE, FALSE), n_unique_peptides = c(3, 1, 2)
  )
  expect_equal(out$n_missing_pair, 1L)  # B has no decoy
})

# =============================================================================
# Enrichment presence rule (method = "enrichment")
# =============================================================================
test_that("calc_taxon_fdr enrichment requires peptide and a positive margin", {
  expect_error(
    calc_taxon_fdr(pep, taxon, decoy, method = "enrichment"),
    "requires `peptide`"
  )
  expect_error(
    calc_taxon_fdr(pep, taxon, decoy, peptide, method = "enrichment", margin = 0),
    "`margin` must be a positive"
  )
})

test_that("enrichment passes taxa whose per-peptide score beats margin x decoy rate", {
  # Decoy reps set the noise rate; median per-peptide decoy score = median(2, 3) = 2.5.
  #   HIGH: score 100 over 20 pep -> 5/pep = 2.0x rate  -> passes at margin 2
  #   LOW : score 6   over 3  pep -> 2/pep = 0.8x rate  -> fails  at margin 2
  tx <- c("HIGH","HIGH","LOW","LOW","dA","dB")
  sc <- c(100,   4,     6,    5,    100, 150)   # dA/dB are decoy-only noise reps
  dc <- c(FALSE, TRUE,  FALSE,TRUE, TRUE, TRUE)
  np <- c(20,    2,     3,    2,    40,  50)
  out <- conduitR:::pick_taxon_fdr_compete(tx, sc, dc, np,
           method = "enrichment", margin = 2, min_peptides = 2)
  reps <- out$results[out$results$picked_winner, ]
  expect_equal(out$method, "enrichment")
  expect_true(is.finite(out$decoy_rate))
  expect_true(reps$pass[reps$taxon == "HIGH"])
  expect_false(reps$pass[reps$taxon == "LOW"])
})

test_that("enrichment still enforces the min_peptides floor", {
  # A taxon with a huge per-peptide score but too few peptides must not pass.
  tx <- c("X","X","dA","dA")
  sc <- c(50, 1, 2, 3)
  dc <- c(FALSE, TRUE, TRUE, TRUE)
  np <- c(3,  1, 10, 10)                    # X has 3 peptides
  out <- conduitR:::pick_taxon_fdr_compete(tx, sc, dc, np,
           method = "enrichment", margin = 2, min_peptides = 10)
  expect_false(out$results$pass[out$results$taxon == "X" & out$results$picked_winner])
})

# =============================================================================
# Regression test — acceptance criteria on the real first-pass table
# =============================================================================
# Fixture is the first_pass_fdr_results.tsv produced for the murine_abx_treatment
# experiment (827 target + 663 decoy rows over 894 distinct lineages). It is the
# correctness spec for the picked-FDR fix (see picked_taxon_fdr.py).
test_that("picked FDR reproduces the documented acceptance criteria", {
  fixture <- testthat::test_path("fixtures", "first_pass_fdr_results.tsv")
  skip_if_not(file.exists(fixture), "regression fixture not available")

  d <- readr::read_tsv(fixture, show_col_types = FALSE)
  out <- conduitR:::pick_taxon_fdr_compete(
    taxon             = as.character(d$taxon),
    score             = d$score,
    decoy             = as.logical(d$decoy),
    n_unique_peptides = d$n_unique_peptides_all,
    qvalue_threshold  = 0.05,
    min_peptides      = 2
  )

  names <- dplyr::distinct(
    dplyr::mutate(d, taxon = as.character(taxon)),
    taxon, taxon_name, taxon_name_rank
  )
  reps <- dplyr::left_join(
    dplyr::filter(out$results, picked_winner), names, by = "taxon"
  )

  # 894 distinct lineages -> 894 representatives.
  expect_equal(nrow(reps), 894L)

  # Poaceae (a true low-abundance dietary family) is rescued: q ~= 0.0075, pass.
  poaceae <- reps[reps$taxon_name == "Poaceae", ]
  expect_false(poaceae$decoy)
  expect_true(poaceae$pass)
  expect_equal(round(poaceae$qvalue, 4), 0.0075)

  # ~133 units pass before any decoy enters: first decoy-winner at rank 134.
  expect_equal(out$first_decoy_rank, 134L)

  # Picking sinks the null to the bottom: decoy-winners score <= ~20,
  # target-winners run up to ~13,570.
  decoy_reps  <- reps[reps$decoy, ]
  target_reps <- reps[!reps$decoy, ]
  expect_lte(max(decoy_reps$score), 20)
  expect_gt(max(target_reps$score), 13000)
  expect_lt(max(decoy_reps$score), max(target_reps$score) / 100)

  # The large-family decoys that broke the old running FDR no longer compete —
  # their targets win the pick and the decoys are discarded.
  for (fam in c("Lachnospiraceae", "Enterobacteriaceae")) {
    expect_false(reps$decoy[reps$taxon_name == fam])
  }

  # ~413 units pass at q<=0.05, min_peptides=2 (187 fam / 141 sp / 64 gen / 21 strain).
  passers <- reps[reps$pass, ]
  expect_equal(nrow(passers), 413L)
  rank_counts <- table(passers$taxon_name_rank)
  expect_equal(as.integer(rank_counts[["family"]]),  187L)
  expect_equal(as.integer(rank_counts[["species"]]), 141L)
  expect_equal(as.integer(rank_counts[["genus"]]),   64L)
  expect_equal(as.integer(rank_counts[["strain"]]),  21L)

  # ~88 families pass at q<=0.01.
  out01 <- conduitR:::pick_taxon_fdr_compete(
    taxon             = as.character(d$taxon),
    score             = d$score,
    decoy             = as.logical(d$decoy),
    n_unique_peptides = d$n_unique_peptides_all,
    qvalue_threshold  = 0.01,
    min_peptides      = 2
  )
  reps01 <- dplyr::left_join(
    dplyr::filter(out01$results, picked_winner, pass), names, by = "taxon"
  )
  expect_equal(nrow(reps01), 133L)
  expect_equal(sum(reps01$taxon_name_rank == "family"), 88L)
})
