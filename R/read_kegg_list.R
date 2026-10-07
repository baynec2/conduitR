#' Read a KEGG REST \code{list} File into a Term-Name Dictionary
#'
#' Reads a two-column TSV as produced by the KEGG REST \code{list} endpoint
#' (e.g. \code{https://rest.kegg.jp/list/ko}, \code{.../list/pathway},
#' \code{.../list/module}, \code{.../list/brite}) and returns a term-name
#' dictionary. The database prefix on each identifier (e.g. \code{ko:},
#' \code{path:}, \code{md:}, \code{br:}) is stripped so the resulting \code{term}
#' matches the accessions stored by the conduitR annotation pipeline
#' (e.g. \code{K00001}, \code{map00010}, \code{M00001}, \code{ko00000}).
#'
#' The result is intended to be tagged with the appropriate
#' \code{annotation_type} (e.g. \code{"kegg_orthology"}, \code{"kegg_map_pathway"},
#' \code{"kegg_module"}, \code{"brite"}) and passed to [add_term_descriptions()].
#'
#' @param path Path to a KEGG \code{list} TSV file (id in column 1, name in
#'   column 2).
#'
#' @return A tibble with columns \code{term} (the prefix-stripped identifier) and
#'   \code{description} (the KEGG name).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ko_dict <- read_kegg_list("resources/annotation/dictionaries/kegg_ko.tsv")
#' ko_dict$annotation_type <- "kegg_orthology"
#' }
read_kegg_list <- function(path) {
  readr::read_tsv(
    path,
    col_names = c("id", "description"),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    # Strip the leading database prefix (lowercase letters up to the first colon).
    dplyr::transmute(
      term        = sub("^[A-Za-z]+:", "", .data$id),
      description = .data$description
    ) |>
    dplyr::filter(!is.na(.data$term), .data$term != "",
                  !is.na(.data$description), .data$description != "") |>
    dplyr::distinct(.data$term, .keep_all = TRUE)
}
