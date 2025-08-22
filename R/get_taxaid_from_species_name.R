#' Get NCBI taxa id corresponding to species name
#'
#' This function will search a string of species names and return the NCBI
#' taxonomy id for that string. If a match isn't found, an informative message
#' and NA is returned. If more than one match is found, an informative message
#' and the first match is returned.
#'
#' @param organism species name (string)
#'
#' @returns NCBI taxonomy id if found, or NA and informative message.
#' @export
#'
#' @examples
#' get_taxaid_from_species_name("homo sapiens")
get_taxid_from_species_name <- function(organism) {
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=", URLencode(organism), "&retmode=json")
  res <- jsonlite::fromJSON(url)
  ids <- res$esearchresult$idlist
  if (length(ids) == 0) {
    message(paste("No match found for:", organism))
    return(NA)
  } else if (length(ids) > 1) {
    message(paste("Multiple matches for:", organism, "- using first hit:", ids[1]))
  }
  return(ids[1])
}
