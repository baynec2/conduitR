#' Parse a GO OBO File into a Term-Name Dictionary
#'
#' Reads a Gene Ontology \code{.obo} file (e.g. \code{go-basic.obo}) and returns a
#' tibble mapping every GO identifier to its term name. Both primary identifiers
#' (\code{id:}) and secondary identifiers (\code{alt_id:}) are emitted, each
#' mapped to the primary term name, so annotations that reference a secondary GO
#' ID are still resolved. Obsolete terms are retained because legacy annotations
#' may still reference them.
#'
#' The result is intended to be tagged with \code{annotation_type = "go"} and
#' passed to [add_term_descriptions()].
#'
#' @param obo_path Path to a GO OBO file (e.g. a date-pinned \code{go-basic.obo}).
#'
#' @return A tibble with columns \code{term} (the GO ID, e.g. \code{"GO:0008150"})
#'   and \code{description} (the term name, e.g. \code{"biological_process"}).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' go_dict <- parse_go_obo("resources/annotation/dictionaries/go-basic.obo")
#' go_dict$annotation_type <- "go"
#' }
parse_go_obo <- function(obo_path) {
  lines <- readLines(obo_path, warn = FALSE)

  # Stanza headers look like "[Term]", "[Typedef]", etc. Number every line by the
  # stanza it belongs to, then keep only the body lines of [Term] stanzas.
  is_header <- grepl("^\\[.*\\]$", lines)
  stanza_id <- cumsum(is_header)
  term_stanzas <- stanza_id[is_header & lines == "[Term]"]

  keep <- stanza_id %in% term_stanzas & !is_header
  body <- tibble::tibble(stanza = stanza_id[keep], line = lines[keep])

  ids <- body |>
    dplyr::filter(startsWith(.data$line, "id: GO:")) |>
    dplyr::transmute(.data$stanza, term = sub("^id: ", "", .data$line))

  names_ <- body |>
    dplyr::filter(startsWith(.data$line, "name: ")) |>
    dplyr::transmute(.data$stanza, description = sub("^name: ", "", .data$line))

  alts <- body |>
    dplyr::filter(startsWith(.data$line, "alt_id: GO:")) |>
    dplyr::transmute(.data$stanza, term = sub("^alt_id: ", "", .data$line))

  stanza_names <- dplyr::inner_join(ids, names_, by = "stanza")

  primary <- stanza_names |>
    dplyr::transmute(.data$term, .data$description)

  secondary <- alts |>
    dplyr::inner_join(
      dplyr::select(stanza_names, "stanza", "description"),
      by = "stanza"
    ) |>
    dplyr::transmute(.data$term, .data$description)

  dplyr::bind_rows(primary, secondary) |>
    dplyr::filter(!is.na(.data$description), .data$description != "") |>
    dplyr::distinct(.data$term, .keep_all = TRUE)
}
