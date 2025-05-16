#' Add Row Data to Limma Analysis Results
#'
#' Combines the row data (feature annotations) from a QFeatures object with the statistical
#' results from a limma analysis. This function merges feature metadata with differential
#' expression results, making it easier to interpret and filter the analysis results
#' based on feature annotations.
#'
#' @param limma_stats A data frame containing the results from limma analysis, typically
#'   the output from \code{\link[perform_limma_analysis]{perform_limma_analysis}}. Must
#'   contain feature IDs as row names.
#' @param qf A QFeatures object containing the original data and feature annotations.
#'   The rowData from the specified assay will be merged with the limma results.
#' @param assay_name Character string specifying which assay's rowData to use for
#'   merging with the limma results.
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item All columns from the original limma results (e.g., logFC, P.Value, adj.P.Val)
#'     \item All columns from the feature annotations (rowData)
#'     \item Feature IDs in a column named "id"
#'   }
#'   The rows are matched by feature ID, and the order is preserved from the limma results.
#'
#' @export
#'
#' @examples
#' # Basic usage with limma results:
#' # results <- perform_limma_analysis(qfeatures_obj, "protein", ~group, "groupB - groupA")
#' # annotated_results <- add_rowdata_to_limma_results(results$top_table, qfeatures_obj, "protein")
#' 
#' # Filter results based on annotations:
#' # significant_proteins <- annotated_results |>
#' #   dplyr::filter(adj.P.Val < 0.05, !is.na(gene_name))
#' 
#' # Group results by annotation:
#' # results_by_pathway <- annotated_results |>
#' #   dplyr::group_by(pathway) |>
#' #   dplyr::summarise(n_significant = sum(adj.P.Val < 0.05))
#'
#' @note
#' This function:
#' \itemize{
#'   \item Uses a right join to preserve all rows from the limma results
#'   \item Converts row names to a column named "id" for merging
#'   \item Maintains the original order of the limma results
#' }
#' If some features in the limma results don't have corresponding annotations,
#' those rows will still be included in the output with NA values for the
#' annotation columns.
#'
#' @seealso \code{\link[perform_limma_analysis]{perform_limma_analysis}} for generating
#'   the limma results that this function processes
add_rowdata_to_limma_results = function(limma_stats, qf, assay_name){

  # Extracting rowData from the assay
  rowData = SummarizedExperiment::rowData(qf[[assay_name]]) |>
    as.data.frame() |>
    tibble::rownames_to_column("id")

  # Adding id column to limma_stats
  limma_stats = limma_stats |>
    tibble::rownames_to_column("id")

  # Combining
  out = dplyr::right_join(rowData, limma_stats, by = "id")

  return(out)
}
