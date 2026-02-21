#' Gene Set Enrichment Analysis (GSEA) on Ranked Limma Results
#'
#' Runs GSEA using `clusterProfiler::GSEA` with a ranked gene list (e.g. by
#' logFC) and custom TERM2GENE/TERM2NAME derived from limma stats and optional
#' conduit annotations.
#'
#' @param limma_stats Tibble of limma results with an ID column and a numeric
#'   ranking column (e.g. `logFC`). Must also contain the column used for
#'   `annotation_type` (e.g. GO terms).
#' @param conduit A conduit object; its `@annotations` slot is used for
#'   TERM2NAME when `annotation_type` matches.
#' @param annotation_type Character. Column name in `limma_stats` and in
#'   conduit annotations for term IDs (default: `"go"`).
#' @param ranking_column Character. Column name in `limma_stats` used to rank
#'   genes (default: `"logFC"`). Values must be numeric, non-NA, with unique
#'   row names.
#'
#' @return Result of `clusterProfiler::GSEA` (enrichResult object).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' limma_res <- perform_limma_analysis(qf, "protein_groups", ~ group, "B - A")
#' gsea_res <- perform_gsea(limma_res$top_table, conduit = conduit_obj,
#'   annotation_type = "go", ranking_column = "logFC")
#' }
perform_gsea = function(limma_stats,
                        conduit,
                        annotation_type = "go",
                        ranking_column = "logFC"){

  # 1. get rid of NAs
  no_na <- limma_stats |>
    dplyr::filter(!is.na(!!rlang::sym(annotation_type))) |>
    dplyr::select(id, !!rlang::sym(ranking_column)) |>
    dplyr::filter(!is.na(!!rlang::sym(ranking_column))) |>
    dplyr::arrange(dplyr::desc(!!rlang::sym(ranking_column)))  # ensure sorted

  # Pull Gene list
  geneList <- no_na |>
    dplyr::pull(!!rlang::sym(ranking_column))
  # Name the gene list vector
  names(geneList) <- dplyr::pull(no_na,id)

  # sanity checks
  stopifnot(is.numeric(geneList))
  stopifnot(all(!is.na(geneList)))
  stopifnot(length(unique(names(geneList))) == length(geneList))

  # 2. Expand TERM2GENE
  all_terms <- limma_stats |>
    dplyr::select(id, !!rlang::sym(annotation_type)) |>
    tidyr::separate_rows(!!rlang::sym(annotation_type), sep = ";") |>
    dplyr::filter(!is.na(!!rlang::sym(annotation_type)), !!rlang::sym(annotation_type) != "") |>
    dplyr::distinct() |>
    dplyr::mutate(term = !!rlang::sym(annotation_type))

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

  # 3. Run clusterProfiler GSEA with custom TERM2GENE
  gsea_res <- clusterProfiler::GSEA(
    geneList  = geneList,
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE)

  return(gsea_res)
}
