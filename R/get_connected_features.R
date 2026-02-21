#' Get Feature IDs Connected to a Given Feature in an Adjacency Matrix
#'
#' From a QFeatures assay that has an adjacency matrix in its metadata, returns
#' the feature IDs (column names) that are connected to the row specified by
#' `id` (i.e. non-zero entries in that row).
#'
#' @param qf A `QFeatures` object.
#' @param i Character. Name of the assay to use (default: `"protein_groups"`).
#' @param am_name Character. Name of the metadata element containing the
#'   adjacency matrix (default: `"adjacencyMatrix"`).
#'
#' @return Character vector of feature IDs connected to the row `id` in the
#'   adjacency matrix. Note: `id` must be available in the calling environment
#'   (e.g. passed from a caller that has it).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After building a QFeatures object with adjacency matrix metadata
#' ids <- get_connected_features(qf, i = "protein_groups")
#' # Use a specific feature id from the assay:
#' feature_id <- rownames(qf[["protein_groups"]])[1]
#' connected <- get_connected_features(qf)
#' # (id would need to be set in scope for the function body to work)
#' }
get_connected_features <- function(qf,i = "protein_groups",
                                   am_name = "adjacencyMatrix"){
  # Extract the adjacency matrix
  adj_mat <- adjacencyMatrix(qf[[i]],am_name)
  # Figure out what corresponds with id
  out <- colnames(adj_mat)[adj_mat[id, ] > 0]
  return(out)
}
