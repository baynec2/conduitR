#' add_rowdata_to_limma_results
#'
#' Add rowdata to the statistics of limma data.
#'
#' @param limma_stats
#' @param qf
#' @param assay_name
#'
#' @returns
#' @export
#'
#' @examples
add_rowdata_to_limma_results = function(limma_stats,qf,assay_name){

  # Extracting rowData from the assay
  rowData = SummarizedExperiment::rowData(qf[[assay_name]]) |>
    as.data.frame() |>
    tibble::rownames_to_column("id")

  # Adding id column to limma_stats
  limma_stats = limma_stats |>
    tibble::rownames_to_column("id")

  # Combining
  out = dplyr::right_join(rowData,limma_stats,by = "id")

  return(out)
}
