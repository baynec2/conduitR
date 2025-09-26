#' get_connected_features
#'
#' @param qf QFeatures object
#' @param i SE Name
#' @param id id that you would like to get connected features for
#'
#' @returns
#' @export
#'
#' @examples
get_connected_features <- function(qf,i = "protein_groups",
                                   am_name = "adjacencyMatrix"){
  # Extract the adjacency matrix
  adj_mat <- adjacencyMatrix(qf[[i]],am_name)
  # Figure out what corresponds with id
  out <- colnames(adj_mat)[adj_mat[id, ] > 0]
  return(out)
}
