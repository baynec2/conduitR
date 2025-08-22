#' Get NCBI taxonomy ids from species names
#'
#' Get the NCBI taxonomy ids corresponding to a vector of species names in a
#' convient tibble.
#'
#' @param species_names a charachter vector of NCBI taxonomy ids.
#'
#' @returns a tibble of the species name and ncbi taxonomy id
#' @export
#'
#' @examples
#' get_taxaids_from_species_names(c("Homo sapiens", "Clostridium symbiosum"))
get_taxaids_from_species_names <- function(species_names) {
  taxids <- purrr::map_chr(species_names, get_taxid_from_species_name)
  result <- tibble::tibble(
    species_name = species_names,
    taxonomy_id = taxids
  )
  return(result)
}



