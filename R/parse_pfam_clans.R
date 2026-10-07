#' Parse a Pfam-A Clans File into a Pfam Term-Name Dictionary
#'
#' Reads the Pfam \code{Pfam-A.clans.tsv} file and returns a tibble mapping each
#' Pfam accession to its description. The file is a headerless five-column TSV:
#' Pfam accession, clan accession, clan ID, Pfam ID and Pfam description.
#'
#' The result is intended to be tagged with \code{annotation_type = "pfam"} and
#' passed to [add_term_descriptions()].
#'
#' @param path Path to a \code{Pfam-A.clans.tsv} file.
#'
#' @return A tibble with columns \code{term} (the Pfam accession, e.g.
#'   \code{"PF00001"}) and \code{description} (the Pfam family description).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pfam_dict <- parse_pfam_clans("resources/annotation/dictionaries/Pfam-A.clans.tsv")
#' pfam_dict$annotation_type <- "pfam"
#' }
parse_pfam_clans <- function(path) {
  readr::read_tsv(
    path,
    col_names = c("pfam_acc", "clan_acc", "clan_id", "pfam_id", "description"),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    dplyr::transmute(term = .data$pfam_acc, description = .data$description) |>
    dplyr::filter(!is.na(.data$term), .data$term != "",
                  !is.na(.data$description), .data$description != "") |>
    dplyr::distinct(.data$term, .keep_all = TRUE)
}
