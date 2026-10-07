#' Aggregate an Assay to Annotation Terms
#'
#' Aggregates features of a base assay into a new assay whose rows are the
#' levels of an annotation column. The annotation column may be either a
#' list-column (e.g. GO term sets, where one feature can map to many terms)
#' or an atomic column (e.g. species, where each feature maps to exactly one
#' term).
#'
#' @param qf A QFeatures object.
#' @param i Name of the base assay to aggregate from.
#' @param fcol Name of a rowData column on assay `i`. Either a list-column
#'   (each entry a character vector of terms) or an atomic character/factor
#'   column.
#' @param name Name of the new aggregated assay (defaults to `fcol`).
#' @param include_na One of `"drop"` (features with no annotation are
#'   excluded from the aggregation; column totals shrink) or `"group"`
#'   (unannotated features are collected into an `"Unassigned"` term so
#'   intensity is conserved).
#' @param na.rm Passed to `fun` for handling NA intensity values during
#'   summation. Independent of `include_na`.
#' @param fun Summation function passed to [QFeatures::aggregateFeatures()].
#'   If `NULL`, defaults to [MsCoreUtils::colSumsMat] for list-columns and
#'   [base::colSums] for atomic columns.
#'
#' @return A QFeatures object with one new assay named `name`.
#' @export
aggregate_assay_by_annotation <- function(qf,
                                          i,
                                          fcol,
                                          name = fcol,
                                          include_na = c("drop", "group"),
                                          na.rm = TRUE,
                                          fun = NULL) {
  include_na <- match.arg(include_na)

  if (!i %in% names(qf)) {
    stop(sprintf("Assay '%s' not found in qf", i))
  }
  rd <- SummarizedExperiment::rowData(qf[[i]])
  if (!fcol %in% colnames(rd)) {
    stop(sprintf("Column '%s' not found in rowData of assay '%s'", fcol, i))
  }

  if (is.list(rd[[fcol]])) {
    .aggregate_by_listcol(qf, i = i, fcol = fcol, name = name,
                          include_na = include_na, na.rm = na.rm,
                          fun = if (is.null(fun)) MsCoreUtils::colSumsMat else fun)
  } else {
    .aggregate_by_atomic(qf, i = i, fcol = fcol, name = name,
                         include_na = include_na, na.rm = na.rm,
                         fun = if (is.null(fun)) colSums else fun)
  }
}

.aggregate_by_listcol <- function(qf, i, fcol, name, include_na, na.rm, fun) {
  rd <- SummarizedExperiment::rowData(qf[[i]])
  feat_terms <- rd[[fcol]]
  feat_ids <- rownames(qf[[i]])

  no_terms <- vapply(feat_terms, function(x) {
    is.null(x) || length(x) == 0 || all(is.na(x))
  }, logical(1))

  real_terms <- unique(stats::na.omit(unlist(feat_terms[!no_terms],
                                             use.names = FALSE)))
  if (length(real_terms) == 0 && !any(no_terms)) {
    stop(sprintf("No annotation terms found in column '%s'", fcol))
  }

  # QFeatures requires every feature to map to at least one term, so unmapped
  # features get a sentinel. In "group" mode the sentinel is the visible
  # "Unassigned" row; in "drop" mode we strip it from the result.
  unmapped_term <- if (include_na == "group") "Unassigned" else ".__drop__"
  if (any(no_terms)) {
    feat_terms[no_terms] <- unmapped_term
    all_terms <- c(real_terms, unmapped_term)
  } else {
    all_terms <- real_terms
  }

  lens <- lengths(feat_terms)
  feat_long <- rep(feat_ids, lens)
  term_long <- unlist(feat_terms, use.names = FALSE)
  keep <- !is.na(term_long)
  feat_long <- feat_long[keep]
  term_long <- term_long[keep]

  adj <- Matrix::sparseMatrix(
    i = match(feat_long, feat_ids),
    j = match(term_long, all_terms),
    x = 1,
    dims = c(length(feat_ids), length(all_terms)),
    dimnames = list(feat_ids, all_terms)
  )

  # Multiple adjacency matrices in rowData confuse downstream dispatch; strip
  # any pre-existing matrix-typed columns before setting the new one.
  is_adj <- vapply(rd, function(x) is.matrix(x) || methods::is(x, "Matrix"),
                   logical(1))
  if (any(is_adj)) {
    SummarizedExperiment::rowData(qf[[i]]) <- rd[, !is_adj, drop = FALSE]
  }

  QFeatures::adjacencyMatrix(qf[[i]]) <- adj
  adj_col_name <- paste0(fcol, "_adjacency_matrix")
  rd_now <- SummarizedExperiment::rowData(qf[[i]])
  hit <- colnames(rd_now) == "adjacencyMatrix"
  if (any(hit)) colnames(rd_now)[hit] <- adj_col_name
  SummarizedExperiment::rowData(qf[[i]]) <- rd_now

  qf <- QFeatures::aggregateFeatures(
    qf,
    i = i,
    name = name,
    fcol = adj_col_name,
    fun = fun,
    na.rm = na.rm
  )

  new_rd <- SummarizedExperiment::rowData(qf[[name]])
  adj_cols <- grepl("_adjacency_matrix$", colnames(new_rd))
  if (any(adj_cols)) {
    SummarizedExperiment::rowData(qf[[name]]) <- new_rd[, !adj_cols, drop = FALSE]
  }

  if (include_na == "drop" && ".__drop__" %in% rownames(qf[[name]])) {
    se <- qf[[name]]
    se <- se[rownames(se) != ".__drop__", ]
    # Replacing the assay drops the auto-generated AssayLink to the base
    # assay; aggregation correctness is unaffected and downstream
    # consumers in conduit-summit use the aggregated assay directly.
    suppressWarnings(qf[[name]] <- se)
  }
  qf
}

.aggregate_by_atomic <- function(qf, i, fcol, name, include_na, na.rm, fun) {
  eff_fcol <- fcol
  added_tmp <- FALSE

  if (include_na == "group") {
    rd <- SummarizedExperiment::rowData(qf[[i]])
    vals <- rd[[fcol]]
    if (anyNA(vals)) {
      tmp_col <- paste0(".", fcol, "_with_unassigned")
      if (is.factor(vals)) {
        vals <- as.character(vals)
      }
      vals[is.na(vals)] <- "Unassigned"
      rd[[tmp_col]] <- vals
      SummarizedExperiment::rowData(qf[[i]]) <- rd
      eff_fcol <- tmp_col
      added_tmp <- TRUE
    }
  }

  qf <- QFeatures::aggregateFeatures(
    qf,
    i = i,
    fcol = eff_fcol,
    name = name,
    fun = fun,
    na.rm = na.rm
  )

  if (added_tmp) {
    rd_base <- SummarizedExperiment::rowData(qf[[i]])
    rd_base[[eff_fcol]] <- NULL
    SummarizedExperiment::rowData(qf[[i]]) <- rd_base
    rd_new <- SummarizedExperiment::rowData(qf[[name]])
    if (eff_fcol %in% colnames(rd_new)) {
      colnames(rd_new)[colnames(rd_new) == eff_fcol] <- fcol
      SummarizedExperiment::rowData(qf[[name]]) <- rd_new
    }
  }
  qf
}
