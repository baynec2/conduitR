#' Convert KEGG Pathway IDs to KO Pathway IDs
#'
#' Replaces the organism code or \code{map} prefix in KEGG pathway IDs with
#' \code{ko} to obtain the KO (KEGG Orthology) pathway ID.  KO pathways are
#' accessible via the public KEGG REST API and their nodes carry KO IDs
#' (e.g. \code{K01320}), making them suitable for joining against the
#' \code{kegg_orthology} column produced by the conduitR annotation pipeline.
#'
#' @param kegg_pathway_ids Character vector of KEGG pathway IDs with an
#'   organism prefix (e.g. \code{"hsa04610"}) or a \code{map} prefix
#'   (e.g. \code{"map04610"}).
#'
#' @return Character vector of the same length with the prefix replaced by
#'   \code{ko}.
#'
#' @export
#'
#' @examples
#' keggpath_to_ko(c("hsa04610", "map00010"))
#' # "ko04610" "ko00010"
keggpath_to_ko <- function(kegg_pathway_ids) {
  sub("^([a-z]+)", "ko", kegg_pathway_ids)
}
