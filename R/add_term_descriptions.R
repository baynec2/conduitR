#' Fill Missing Annotation Descriptions from Authoritative Dictionaries
#'
#' Populates the \code{description} column of a long-format annotation table by
#' joining against one or more authoritative term-name dictionaries (e.g. the GO
#' OBO, KEGG list files, ENZYME, Pfam clans). Only blank descriptions are filled;
#' descriptions that are already present are never overwritten, so the source of
#' each description stays unambiguous. There is intentionally **no** fallback to
#' descriptions found elsewhere in the same table (e.g. borrowing a UniProt name
#' for an eggNOG term) — a term is described only if its own dictionary names it,
#' otherwise its description remains \code{NA}.
#'
#' Matching is exact on the pair \code{(annotation_type, term)}. Build the
#' dictionary with the same \code{annotation_type} values used in
#' \code{annotations} (e.g. tag the GO OBO with \code{annotation_type = "go"}).
#'
#' @param annotations A data frame / tibble in conduit long annotation format.
#'   Must contain \code{annotation_type} and \code{term}; a \code{description}
#'   column is added if absent. Any other columns (e.g. \code{protein_id},
#'   \code{Protein.Group}) are preserved.
#' @param dictionary A data frame / tibble with columns \code{annotation_type},
#'   \code{term} and \code{description} giving the authoritative name for each
#'   term. Rows with a blank dictionary description are ignored, and duplicate
#'   \code{(annotation_type, term)} pairs keep their first description.
#'
#' @return \code{annotations} with the same columns and row order, with blank
#'   \code{description} values filled where the dictionary provides a name.
#'
#' @export
#'
#' @examples
#' annotations <- tibble::tibble(
#'   protein_id      = c("P1", "P2"),
#'   annotation_type = c("go", "go"),
#'   term            = c("GO:0008150", "GO:0003674"),
#'   description     = c(NA_character_, NA_character_)
#' )
#' dictionary <- tibble::tibble(
#'   annotation_type = "go",
#'   term            = c("GO:0008150", "GO:0003674"),
#'   description     = c("biological_process", "molecular_function")
#' )
#' add_term_descriptions(annotations, dictionary)
add_term_descriptions <- function(annotations, dictionary) {
  if (!all(c("annotation_type", "term") %in% names(annotations))) {
    stop("`annotations` must contain `annotation_type` and `term` columns.")
  }
  if (!all(c("annotation_type", "term", "description") %in% names(dictionary))) {
    stop("`dictionary` must contain `annotation_type`, `term` and `description` columns.")
  }

  # Ensure a description column exists and treat empty strings as missing.
  if (!"description" %in% names(annotations)) {
    annotations[["description"]] <- NA_character_
  }
  annotations <- annotations |>
    dplyr::mutate(description = dplyr::na_if(as.character(.data$description), ""))

  # One name per (annotation_type, term); drop blank dictionary entries.
  dict <- dictionary |>
    dplyr::transmute(
      annotation_type    = as.character(.data$annotation_type),
      term               = as.character(.data$term),
      .dict_description   = as.character(.data$description)
    ) |>
    dplyr::filter(!is.na(.data$.dict_description), .data$.dict_description != "") |>
    dplyr::distinct(.data$annotation_type, .data$term, .keep_all = TRUE)

  annotations |>
    dplyr::mutate(
      annotation_type = as.character(.data$annotation_type),
      term            = as.character(.data$term)
    ) |>
    dplyr::left_join(dict, by = c("annotation_type", "term")) |>
    # Keep any existing description; only fall back to the dictionary name.
    dplyr::mutate(description = dplyr::coalesce(.data$description, .data$.dict_description)) |>
    dplyr::select(-".dict_description")
}
