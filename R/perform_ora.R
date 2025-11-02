#' perform_ora
#'
#' @param limma_stats
#' @param direction
#' @param conduit
#' @param annotation_type
#' @param adj_pval_threshold
#' @param logFC_threshold
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
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
