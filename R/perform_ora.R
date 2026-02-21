#' Over-Representation Analysis (ORA) on Limma Results
#'
#' Takes differential expression results from limma, filters by significance and
#' direction (up/down), and runs gene set enrichment (ORA) using
#' `clusterProfiler::enricher` with term-to-gene and term-to-name mappings
#' derived from the limma stats and optional conduit annotations.
#'
#' @param limma_stats Tibble or data frame of limma results; must include
#'   columns for adjusted p-value, logFC, protein/feature ID, and the
#'   annotation type column (e.g. species).
#' @param direction Character. `"up"` or `"down"` to test over-representation of
#'   upregulated or downregulated features (default: `"up"`).
#' @param conduit A conduit object whose `@annotations` slot can supply
#'   TERM2NAME when `annotation_type` matches.
#' @param annotation_type Character. Name of the column in `limma_stats` (and
#'   in conduit annotations) used as the term (e.g. `"species"`) (default:
#'   `"species"`).
#' @param adj_pval_threshold Numeric. Maximum adjusted p-value for significant
#'   features (default: 0.05).
#' @param logFC_threshold Numeric. Minimum absolute log fold change for
#'   inclusion; must be >= 0 for `direction = "up"` and <= 0 for `"down"`
#'   (default: 1).
#' @param ... Additional arguments passed to `clusterProfiler::enricher`.
#'
#' @return Result of `clusterProfiler::enricher` (enrichResult object).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After perform_limma_analysis and with a conduit object
#' limma_res <- perform_limma_analysis(qf, "protein_groups", ~ group, "B - A")
#' ora_up <- perform_ora(
#'   limma_res$top_table,
#'   direction = "up",
#'   conduit = conduit_obj,
#'   annotation_type = "species"
#' )
#' }
perform_ora <- function(limma_stats,
                        direction = "up",
                        conduit,
                        annotation_type = "species",
                        adj_pval_threshold = 0.05,
                        logFC_threshold = 1,
                        ...) {

  # Check direction
  if (!direction %in% c("up", "down")) {
    stop("`direction` must be either 'up' or 'down'")
  }

  # Check that logFC_threshold is consistent with direction
  if (direction == "up" && logFC_threshold < 0) {
    stop("For direction = 'up', logFC_threshold must be >= 0")
  }
  if (direction == "down" && logFC_threshold > 0) {
    stop("For direction = 'down', logFC_threshold must be <= 0")
  }

  # Filter by significance and direction
  sig_proteins <- limma_stats |>
    dplyr::filter(!is.na(adj.P.Val)) |>
    dplyr::filter(adj.P.Val <= adj_pval_threshold) |>
    dplyr::filter(!is.na(logFC)) |>
    dplyr::filter(
      (direction == "up" & logFC >= logFC_threshold) |
      (direction == "down" & logFC <= -logFC_threshold)
    ) |>
    dplyr::pull(Protein.Group)  # or whatever column holds protein IDs

  # 2. Expand TERM2GENE
  all_terms <- limma_stats |>
    dplyr::select(id, !!rlang::sym(annotation_type)) |>
    tidyr::separate_rows(!!rlang::sym(annotation_type), sep = ";") |>
    dplyr::filter(!is.na(!!rlang::sym(annotation_type)), !!rlang::sym(annotation_type) != "") |>
    dplyr::distinct() |>
    dplyr::mutate(term = !!rlang::sym(annotation_type))

  # Getting TERM2GENE
  TERM2GENE = all_terms |>
    dplyr::select(term = term,gene = id) |>
    as.data.frame()

  if(annotation_type %in% unique(conduit@annotations$annotation_type)) {

    TERM2NAME = conduit@annotations |>
      dplyr::filter(annotation_type == annotation_type) |>
      dplyr::select(term,name = description)

  } else {
    TERM2NAME = all_terms |>
      dplyr::select(term = term,name = term) |>
      as.data.frame()
  }

  ora_result <- clusterProfiler::enricher(gene = sig_proteins,
                                        TERM2GENE = TERM2GENE,
                                        TERM2NAME = TERM2NAME,
                                        ...)

  return(ora_result)
}
