#' Calculate taxonomic FDR using picked target-decoy competition
#'
#' Aggregates PSM scores to the taxon level and applies the \strong{picked}
#' target-decoy strategy (Savitski et al. 2015) to call which taxa are present.
#' Each taxon accumulates a score (sum of \code{-log(PEP)}) separately for its
#' target PSMs and its decoy PSMs; the higher-scoring of the two is kept as that
#' taxon's single representative and the loser is discarded, then all
#' representatives compete in one descending-score ranking.
#'
#' @details
#' \strong{Why picked.} A naive running FDR that ranks every taxon's target row
#' \emph{and} its decoy row in one list and walks a running \code{decoys /
#' targets} down it is biased when score scales with abundance: a high-abundance
#' lineage's reversed decoy can outscore a genuinely present low-abundance
#' lineage's forward target, inflating the running FDR before the low-abundance
#' true taxon is reached and wrongly rejecting it. Picking removes this — a
#' present lineage's target wins the pick, so its decoy is discarded and can no
#' longer outrank smaller true lineages. The decoys that survive are genuine
#' noise units that score near zero regardless of rank, so taxa resolved at
#' mixed ranks (family / genus / species / strain) can compete in one list
#' without per-rank stratification (picking sinks the null to the bottom).
#'
#' \strong{Pooling across runs:} This function is designed to operate on PSMs
#' pooled from all runs in an experiment to produce experiment-level taxon
#' detection calls. DIA-NN reports one row per (precursor, run), so a peptide
#' detected in N runs contributes N scores; taxa consistently detected across
#' many samples naturally accumulate higher scores. For per-sample detection,
#' filter the input vectors to a single run before calling.
#'
#' @param pep numeric vector of PEP scores (0-1) for each PSM. Must come from
#'   \strong{unfiltered} DIA-NN output (\code{--qvalue 1}); pre-filtered output
#'   strips most decoys and produces an unrepresentative null distribution.
#' @param taxon character vector of taxon assignment for each PSM.
#' @param decoy logical vector indicating decoy status for each PSM. DIA-NN
#'   reports this as an integer (0/1); coerce with \code{as.logical()} before
#'   passing. Mixing of target and decoy PSMs within the same taxon is expected
#'   (DIA-NN decoys share the taxon label of their source peptide) and is
#'   handled correctly: scores are aggregated separately per
#'   \code{(taxon, decoy)} combination before the pick.
#' @param peptide optional character vector of peptide sequences for collapsing
#'   to best PSM per unique \code{(peptide, taxon, decoy)} combination before
#'   aggregation. Required for the \code{min_peptides} gate to apply (without it
#'   the per-taxon peptide count is unknown and the gate is skipped).
#' @param qvalue_threshold numeric q-value cutoff for a presence call. Default
#'   0.05 — permissive on purpose, because a taxon dropped here is unrecoverable
#'   in downstream database-reduction stages.
#' @param min_peptides numeric minimum distinct peptides
#'   (\code{n_unique_peptides_all}) for a presence call. Default 2 — a single
#'   chance peptide must not make a call. Only enforced when \code{peptide} is
#'   supplied.
#'
#' @returns a list with elements:
#'   \item{results}{tibble with one row per \code{(taxon, decoy)} combination,
#'     containing \code{score}, \code{n_unique_peptides_all} (NA when
#'     \code{peptide} is not provided), \code{decoy}, \code{picked_winner}
#'     (logical — this row won its taxon's pick), the picked \code{fdr} and
#'     \code{qvalue} (NA for discarded losers), and \code{pass} (target winner
#'     AND \code{qvalue <= qvalue_threshold} AND
#'     \code{n_unique_peptides_all >= min_peptides}). Winners sorted by
#'     descending score, then the discarded losers.}
#'   \item{detected}{tibble subset of \code{results} with \code{pass == TRUE}.}
#'   \item{n_targets}{count of target representatives (target winners).}
#'   \item{n_decoys}{count of decoy representatives (decoy winners).}
#'   \item{first_decoy_rank}{1-based rank of the first decoy winner in the
#'     score-sorted representative list (\code{NA} if no decoy wins).}
#'   \item{n_missing_pair}{count of taxa present on only one side (target-only
#'     or decoy-only); these cannot be truly competed and are flagged for the
#'     caller.}
#'   \item{qvalue_threshold}{the q-value threshold used.}
#'   \item{min_peptides}{the min-peptides floor used.}
#' @export
#'
#' @examples
#' # Two taxa, each with a target and a decoy PSM. Taxon A's target wins its
#' # pick; taxon B's decoy outscores its target so B's decoy becomes the
#' # representative and sinks to the bottom of the competition.
#' calc_taxon_fdr(
#'   pep     = c(0.01, 0.80, 0.50, 0.02),
#'   taxon   = c("A",  "A",  "B",  "B"),
#'   decoy   = c(FALSE, TRUE, FALSE, TRUE),
#'   peptide = c("PEPA1", "PEPA2", "PEPB1", "PEPB2")
#' )
calc_taxon_fdr <- function(pep, taxon, decoy, peptide = NULL,
                           qvalue_threshold = 0.05, min_peptides = 2) {

  # --- input validation ---
  if (!is.numeric(pep))     stop("`pep` must be a numeric vector")
  if (!is.character(taxon)) stop("`taxon` must be a character vector")
  if (!is.logical(decoy))   stop("`decoy` must be a logical vector")
  if (!is.numeric(qvalue_threshold) || qvalue_threshold <= 0 ||
      qvalue_threshold >= 1) {
    stop("`qvalue_threshold` must be a numeric value between 0 and 1")
  }
  if (!is.numeric(min_peptides) || min_peptides < 0) {
    stop("`min_peptides` must be a non-negative numeric value")
  }

  lengths <- c(length(pep), length(taxon), length(decoy))
  if (!is.null(peptide)) lengths <- c(lengths, length(peptide))
  if (length(unique(lengths)) != 1) {
    stop("all input vectors must be the same length")
  }

  if (any(pep < 0 | pep > 1, na.rm = TRUE)) {
    stop("`pep` values must be between 0 and 1")
  }

  if (any(is.na(pep) | is.na(taxon) | is.na(decoy))) {
    stop("input vectors must not contain NA values")
  }

  # --- clamp PEP to avoid Inf scores from log(0) ---
  pep <- pmax(pep, .Machine$double.eps)

  # --- build working tibble ---
  df <- tibble::tibble(pep = pep, taxon = taxon, decoy = decoy)

  if (!is.null(peptide)) {
    df <- df |>
      dplyr::mutate(peptide = peptide) |>
      # collapse to best PSM per unique (peptide, taxon, decoy) combination
      dplyr::summarise(
        pep = min(pep),
        .by = c(taxon, peptide, decoy)
      )
  }

  # --- aggregate to taxon level (separately per decoy status) ---
  # n_unique_peptides_all is meaningful only when peptide info is supplied — it
  # counts every distinct peptide sequence contributing to the (taxon, decoy)
  # group's score, regardless of per-PSM evidence strength (weak PSMs already
  # contribute little to score via -log(PEP) weighting). When peptide is NULL,
  # set it to NA so the min_peptides gate is skipped (it cannot be evaluated).
  agg_specs <- list(score = quote(sum(-log(pep))))
  if (!is.null(peptide)) {
    agg_specs$n_unique_peptides_all <- quote(dplyr::n_distinct(peptide))
  }

  agg <- df |>
    dplyr::summarise(!!!agg_specs, .by = c(taxon, decoy))

  if (is.null(peptide)) {
    agg$n_unique_peptides_all <- NA_integer_
  }

  # --- picked target-decoy competition over the per-taxon representatives ---
  pick_taxon_fdr_compete(
    taxon             = agg$taxon,
    score             = agg$score,
    decoy             = agg$decoy,
    n_unique_peptides = agg$n_unique_peptides_all,
    qvalue_threshold  = qvalue_threshold,
    min_peptides      = min_peptides
  )
}

# Picked target-decoy competition over an already-aggregated per-(taxon, decoy)
# score table. Internal: callers reach it through calc_taxon_fdr (PSM-level).
# Kept as a standalone function so the competition can be unit-tested directly
# against an aggregated reference table (the picked-FDR correctness spec).
#
# taxon/score/decoy/n_unique_peptides are parallel vectors, one per aggregated
# (taxon, decoy) row. Returns the same list shape as calc_taxon_fdr.
pick_taxon_fdr_compete <- function(taxon, score, decoy, n_unique_peptides,
                                   qvalue_threshold = 0.05, min_peptides = 2) {

  df <- tibble::tibble(
    taxon                 = as.character(taxon),
    score                 = as.numeric(score),
    decoy                 = as.logical(decoy),
    n_unique_peptides_all = n_unique_peptides
  )

  # --- pick one representative per taxon: the higher-scoring of its rows ---
  # Sort by descending score; tie-break toward the target (decoy = FALSE) so a
  # target never loses a pick to its own decoy on an exact-score tie. The first
  # row per taxon in this order is the winner; the rest are discarded losers.
  df <- df |>
    dplyr::arrange(dplyr::desc(score), decoy, taxon)
  df$picked_winner <- !duplicated(df$taxon)

  # Flag taxa present on only one side (no true target-vs-decoy contest).
  n_missing_pair <- sum(table(df$taxon) < 2L)

  # --- combined competition over the representatives, score-descending ---
  rep <- df[df$picked_winner, , drop = FALSE]  # already score-descending
  is_decoy <- rep$decoy
  cum_d <- cumsum(is_decoy)
  cum_t <- pmax(cumsum(!is_decoy), 1)
  fdr   <- (cum_d + 1) / cum_t              # +1: conservative small-count guard
  qval  <- rev(cummin(rev(fdr)))            # monotone q-value (cummin from tail)
  # n_unique_peptides_all is NA when no peptide info was supplied -> the
  # min_peptides gate cannot be evaluated, so it is skipped (treated as passed).
  npep_ok <- is.na(rep$n_unique_peptides_all) |
    rep$n_unique_peptides_all >= min_peptides
  pass <- (!is_decoy) & (qval <= qvalue_threshold) & npep_ok

  rep$fdr    <- fdr
  rep$qvalue <- qval
  rep$pass   <- pass

  # --- losers: carry through with NA fdr/qvalue and pass = FALSE ---
  losers <- df[!df$picked_winner, , drop = FALSE]
  if (nrow(losers) > 0) {
    losers$fdr    <- NA_real_
    losers$qvalue <- NA_real_
    losers$pass   <- FALSE
  }

  results <- dplyr::bind_rows(rep, losers) |>
    dplyr::select(taxon, score, n_unique_peptides_all, decoy,
                  picked_winner, fdr, qvalue, pass)

  list(
    results          = results,
    detected         = dplyr::filter(results, pass),
    n_targets        = sum(rep$decoy == FALSE),
    n_decoys         = sum(rep$decoy == TRUE),
    first_decoy_rank = if (any(is_decoy)) which(is_decoy)[1] else NA_integer_,
    n_missing_pair   = as.integer(n_missing_pair),
    qvalue_threshold = qvalue_threshold,
    min_peptides     = min_peptides
  )
}
