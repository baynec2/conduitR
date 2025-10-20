#' Title
#'
#' @param limma_stats
#' @param conduit
#' @param annotation_type
#' @param ranking_column
#'
#' @returns
#' @export
#'
#' @examples
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
