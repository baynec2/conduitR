#' Parse an ENZYME \code{enzyme.dat} File into an EC Number Dictionary
#'
#' Reads the ExPASy ENZYME database flat file (\code{enzyme.dat}) and returns a
#' tibble mapping each EC number to its accepted enzyme name. Each record begins
#' with an \code{ID} line carrying the EC number and a \code{DE} line carrying the
#' description; the first \code{DE} line of each record is used and a trailing
#' period is removed.
#'
#' The result is intended to be tagged with \code{annotation_type = "ec_number"}
#' and passed to [add_term_descriptions()].
#'
#' @param path Path to an ENZYME \code{enzyme.dat} flat file.
#'
#' @return A tibble with columns \code{term} (the EC number, e.g. \code{"1.1.1.1"})
#'   and \code{description} (the enzyme name, e.g. \code{"Alcohol dehydrogenase"}).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ec_dict <- parse_enzyme_dat("resources/annotation/dictionaries/enzyme.dat")
#' ec_dict$annotation_type <- "ec_number"
#' }
parse_enzyme_dat <- function(path) {
  lines <- readLines(path, warn = FALSE)

  is_id <- startsWith(lines, "ID   ")
  record <- cumsum(is_id)

  keep <- record > 0
  body <- tibble::tibble(record = record[keep], line = lines[keep])

  ids <- body |>
    dplyr::filter(startsWith(.data$line, "ID   ")) |>
    dplyr::transmute(.data$record, term = trimws(sub("^ID   ", "", .data$line)))

  des <- body |>
    dplyr::filter(startsWith(.data$line, "DE   ")) |>
    # Keep only the first DE line of each record.
    dplyr::distinct(.data$record, .keep_all = TRUE) |>
    dplyr::transmute(
      .data$record,
      description = sub("\\.$", "", trimws(sub("^DE   ", "", .data$line)))
    )

  dplyr::inner_join(ids, des, by = "record") |>
    dplyr::filter(!is.na(.data$description), .data$description != "") |>
    dplyr::transmute(.data$term, .data$description) |>
    dplyr::distinct(.data$term, .keep_all = TRUE)
}
