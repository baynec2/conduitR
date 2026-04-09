#' Calculate taxonomic FDR using target-decoy competition
#'
#' Applies a target-decoy competition strategy at the taxon level. Each taxon
#' accumulates a score (sum of \code{-log(PEP)}) separately for its target PSMs
#' and its decoy PSMs, producing one target row and one decoy row per taxon.
#' All rows compete in the same descending-score ranking to compute a running
#' FDR and q-value.
#'
#' @param pep numeric vector of PEP scores (0-1) for each PSM. Must come from
#'   \strong{unfiltered} DIA-NN output (\code{--qvalue 1}); pre-filtered output
#'   strips most decoys and produces an unrepresentative null distribution.
#' @param taxon character vector of taxonomic family assignment for each PSM.
#' @param decoy logical vector indicating decoy status for each PSM. DIA-NN
#'   reports this as an integer (0/1); coerce with \code{as.logical()} before
#'   passing. Mixing of target and decoy PSMs within the same taxon is expected
#'   (DIA-NN decoys share the family label of their source peptide) and is
#'   handled correctly: scores are aggregated separately per
#'   \code{(taxon, decoy)} combination.
#' @param peptide optional character vector of peptide sequences for collapsing
#'   to best PSM per unique \code{(peptide, taxon, decoy)} combination before
#'   aggregation. If NULL, PSMs are aggregated directly.
#' @param fdr_threshold numeric FDR threshold for calling taxa as detected.
#'   Defaults to 0.01.
#'
#' @details
#' \strong{Pooling across runs:} This function is designed to operate on PSMs
#' pooled from all runs in an experiment to produce experiment-level taxon
#' detection calls. DIA-NN reports one row per (precursor, run), so a peptide
#' detected in N runs contributes N scores; families consistently detected
#' across many samples naturally accumulate higher scores than those seen in
#' only one sample. For per-sample detection, filter the input vectors to a
#' single run before calling.
#'
#' The \code{peptide} argument controls within-experiment redundancy: passing
#' \code{Stripped.Sequence} collapses each unique peptide sequence to its best
#' PEP across all runs (one contribution per peptide regardless of how many
#' runs detected it). Omitting \code{peptide} (default) lets every
#' (precursor, run) detection contribute independently.
#'
#' @returns a list with elements:
#'   \item{results}{tibble with one row per \code{(taxon, decoy)} combination,
#'     containing score, decoy status, running FDR, and q-value, sorted by
#'     descending score}
#'   \item{detected}{tibble subset of results containing only target taxa
#'     passing fdr_threshold}
#'   \item{n_targets}{integer count of target taxa in results}
#'   \item{n_decoys}{integer count of decoy taxa in results}
#'   \item{fdr_threshold}{the FDR threshold used}
#' @export
#'
#' @examples
calc_taxon_fdr <- function(pep, taxon, decoy, peptide = NULL, fdr_threshold = 0.01) {

  # --- input validation ---=
  if (!is.numeric(pep))   stop("`pep` must be a numeric vector")
  if (!is.character(taxon)) stop("`taxon` must be a character vector")
  if (!is.logical(decoy)) stop("`decoy` must be a logical vector")
  if (!is.numeric(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold >= 1) {
    stop("`fdr_threshold` must be a numeric value between 0 and 1")
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

  results <- df |>
    dplyr::summarise(
      score = sum(-log(pep)),
      .by   = c(taxon, decoy)
    ) |>
    dplyr::arrange(dplyr::desc(score)) |>
    dplyr::mutate(
      n_targets = cumsum(!decoy),
      n_decoys  = cumsum(decoy),
      fdr = dplyr::case_when(
        n_targets == 0 ~ 1,
        .default = n_decoys / n_targets
      ),
      qvalue = rev(cummin(rev(fdr)))
    ) |>
    dplyr::select(taxon, score, decoy, fdr, qvalue)

  # --- assemble output ---

  list(
    results       = results,
    detected      = dplyr::filter(results, !decoy, qvalue <= fdr_threshold),
    n_targets     = sum(!results$decoy),
    n_decoys      = sum(results$decoy),
    fdr_threshold = fdr_threshold
  )
}
