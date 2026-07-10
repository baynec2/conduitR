#' Validate a downloaded search-database manifest for completeness
#'
#' Guards against a silently empty or grossly truncated search database (see
#' issue #20). Given the per-proteome download manifest produced by
#' [get_fasta_file()], this errors when no sequences were delivered at all, or
#' when the total delivered count is grossly short of the summed expected
#' `proteinCount`, and warns when some (but not all) proteomes failed to
#' download completely.
#'
#' @param download_results A tibble with (at least) columns `proteome_id`,
#'   `source`, `n_sequences`, and `expected`, as returned (row-bound) by
#'   [get_fasta_file()].
#' @param tol Numeric in (0, 1]. The database is rejected if the total delivered
#'   sequence count is below `tol` times the total expected `proteinCount`.
#'   Default 0.9. The expected-count check is skipped when no expected counts
#'   are available.
#'
#' @return Invisibly `TRUE` if the database passes. Throws an error (empty or
#'   grossly incomplete) or emits a warning (partial failures) otherwise.
#'
#' @keywords internal
#' @noRd
check_search_db_completeness <- function(download_results, tol = 0.9) {
  n_proteomes <- nrow(download_results)
  delivered <- suppressWarnings(as.numeric(download_results$n_sequences))
  delivered[is.na(delivered)] <- 0
  expected <- suppressWarnings(as.numeric(download_results$expected))

  total_delivered <- sum(delivered)
  total_expected <- sum(expected, na.rm = TRUE)

  # Per-proteome manifest for actionable error/warning messages.
  manifest <- paste(
    sprintf(
      "  %s: expected=%s delivered=%d source=%s",
      download_results$proteome_id,
      ifelse(is.na(expected), "?", as.character(expected)),
      as.integer(delivered),
      download_results$source
    ),
    collapse = "\n"
  )

  failed <- download_results$source == "not_downloaded" | delivered == 0
  if (any(failed)) {
    warning(sprintf(
      "%d of %d proteomes did not download completely:\n%s",
      sum(failed), n_proteomes, manifest
    ))
  }

  if (total_delivered == 0) {
    stop(sprintf(
      "Search database is EMPTY: 0 sequences delivered across %d proteome(s). Refusing to write an empty database.\n%s",
      n_proteomes, manifest
    ))
  }

  if (total_expected > 0 && total_delivered < tol * total_expected) {
    stop(sprintf(
      "Search database is grossly incomplete: delivered %d of ~%d expected proteins (%d%%, below the %d%% floor).\n%s",
      as.integer(total_delivered), as.integer(total_expected),
      round(100 * total_delivered / total_expected), round(100 * tol),
      manifest
    ))
  }

  invisible(TRUE)
}
