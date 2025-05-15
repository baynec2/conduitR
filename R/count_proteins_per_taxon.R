#' count_proteins_per_taxon
#'
#' count the number of proteins per taxon
#'
#' @param tibble a tibble with columns organism_type, domain, kingdom,
#' phylum, class, order, family, genus, and species
#'
#' @returns
#' @export
#'
#' @examples
count_proteins_per_taxon = function(tibble){

  sum = tibble |>
    dplyr::group_by(.data$organism_type,.data$domain,.data$kingdom,
                    .data$phylum,.data$class,.data$order,.data$family,
                    .data$genus,.data$species) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::ungroup()

  return(sum)

}
