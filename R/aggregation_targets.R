#' Get registered aggregation targets
#'
#' Returns the registry of annotation columns that have been added to a
#' QFeatures object and are eligible for on-demand aggregation via
#' [aggregate_assay_by_annotation()]. Stale entries (those whose base assay
#' or column no longer exists in the object) are silently filtered out.
#'
#' Targets are registered by [add_annotation_to_qf()] (with `kind =
#' "functional"`) and [add_taxonomy_to_qf()] (with `kind = "taxonomic"`).
#'
#' @param qf A QFeatures object.
#'
#' @return A named list. Each element corresponds to a rowData column name
#'   that can be aggregated against, and contains:
#'   \itemize{
#'     \item `kind`: one of `"functional"`, `"taxonomic"` (more values may
#'           be added in the future).
#'     \item `from`: the name of the base assay this target is aggregated
#'           from.
#'   }
#'   Returns an empty list if nothing has been registered.
#'
#' @export
aggregation_targets <- function(qf) {
  tgts <- S4Vectors::metadata(qf)[["aggregation_targets"]]
  if (is.null(tgts) || length(tgts) == 0) {
    return(list())
  }
  valid <- vapply(names(tgts), function(nm) {
    spec <- tgts[[nm]]
    if (!is.list(spec) || is.null(spec$from)) return(FALSE)
    if (!spec$from %in% names(qf)) return(FALSE)
    nm %in% colnames(SummarizedExperiment::rowData(qf[[spec$from]]))
  }, logical(1))
  tgts[valid]
}

# Register a single aggregation target on a QFeatures object. Idempotent:
# re-registering the same target overwrites the prior entry. Not exported.
register_aggregation_target <- function(qf, target, kind, from) {
  md <- S4Vectors::metadata(qf)
  if (is.null(md[["aggregation_targets"]])) {
    md[["aggregation_targets"]] <- list()
  }
  md[["aggregation_targets"]][[target]] <- list(kind = kind, from = from)
  S4Vectors::metadata(qf) <- md
  qf
}
