#' Title
#'
#' @param kegg_pathway_ids
#'
#' @returns
#' @export
#'
#' @examples
keggpath_to_map <- function(kegg_pathway_ids){
  # Replacing organism code with map.
  kegg_map_ids <- sub("^[a-z]{3}","map",kegg_pathway_ids)
  return(kegg_map_ids)
}
