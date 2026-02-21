#' Convert KEGG Organism Pathway IDs to Reference Map IDs
#'
#' Replaces the three-letter organism code in KEGG pathway IDs (e.g. `hsa`, `eco`)
#' with `map` to obtain the reference pathway map ID (e.g. `hsa04110` -> `map04110`).
#' Useful for pathway visualization or lookup in the reference KEGG map.
#'
#' @param kegg_pathway_ids Character vector of KEGG pathway IDs in form
#'   `org:number` or `org12345` (e.g. `hsa04110`, `eco00010`).
#'
#' @return Character vector of the same length with organism code replaced by
#'   `map`.
#'
#' @export
#'
#' @examples
#' keggpath_to_map(c("hsa04110", "eco00010"))
#' # "map04110" "map00010"
keggpath_to_map <- function(kegg_pathway_ids){
  # Replacing organism code with map.
  kegg_map_ids <- sub("^[a-z]{3}","map",kegg_pathway_ids)
  return(kegg_map_ids)
}
