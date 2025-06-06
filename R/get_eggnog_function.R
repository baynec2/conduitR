#' Get Functional Annotations from EggNOG Database
#'
#' Retrieves Gene Ontology (GO) terms and other functional annotations for proteins
#' using their EggNOG ortholog group IDs. This function queries the EggNOG API to
#' obtain detailed functional information including GO terms, evidence codes, and
#' annotation statistics.
#'
#' @param eggnog_ids Character vector containing one or more EggNOG ortholog group IDs
#'   (e.g., "ENOG502S5P5"). These IDs can be obtained from EggNOG-mapper or other
#'   orthology prediction tools.
#'
#' @return A tibble containing functional annotations with the following columns:
#'   \itemize{
#'     \item eggnog_id: The input EggNOG ortholog group ID
#'     \item go_type: Type of GO term (e.g., "biological_process", "molecular_function")
#'     \item go_id: The GO term identifier
#'     \item go_description: Description of the GO term
#'     \item evidence_codes: GO evidence codes supporting the annotation
#'     \item seq_count: Number of sequences in the ortholog group
#'     \item frequency: Frequency of the GO term in the ortholog group
#'     \item annotation_count: Number of times the GO term is annotated
#'   }
#'   If an EggNOG ID is not found or returns an error, that ID will be omitted from
#'   the results with a warning message.
#'
#' @export
#'
#' @examples
#' # Get functional annotations for a single protein:
#' # insulin_annotations <- get_eggnog_function("ENOG502S5P5")
#' 
#' # Get annotations for multiple proteins:
#' # protein_annotations <- get_eggnog_function(
#' #   c("ENOG502S5P5", "ENOG410XNJK", "ENOG410XNJL")
#' # )
#' 
#' # Use the results for downstream analysis:
#' # - Filter by GO type
#' # - Analyze most common functions
#' # - Compare functional profiles
#'
#' @note
#' This function:
#' \itemize{
#'   \item Requires an internet connection to access the EggNOG API
#'   \item May take some time for multiple IDs due to API rate limiting
#'   \item Returns NA for any fields that cannot be retrieved
#'   \item Handles API errors gracefully with warning messages
#' }
#' The EggNOG API is free to use but has rate limits. For bulk queries, consider
#' implementing appropriate delays between requests.
get_eggnog_function <- function(eggnog_ids) {
  # Initialize an empty results tibble
  results <- tibble::tibble(
    eggnog_id = character(), go_type = character(), go_id = character(),
    go_description = character(), evidence_codes = character(),
    seq_count = numeric(), frequency = numeric(), annotation_count = numeric()
  )

  # Loop through each EggNOG ID
  for (id in eggnog_ids) {
    url <- paste0("http://eggnogapi5.embl.de/nog_data/json/go_terms/", id)
    response <- httr::GET(url)

    if (httr::status_code(response) == 200) {
      data <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"))

      if ("go_terms" %in% names(data)) {
        col_names <- c(
          "go_id", "go_description", "evidence_codes",
          "seq_count", "frequency", "annotation_count"
        )

        # Extract GO terms per category
        go_data <- list(
          MF = data$go_terms$`Molecular Function`,
          BP = data$go_terms$`Biological Process`,
          CC = data$go_terms$`Cellular Component`
        )

        go_df <- dplyr::bind_rows(lapply(names(go_data), function(go_type) {
          if (!is.null(go_data[[go_type]])) {
            as.data.frame(go_data[[go_type]]) |>
              setNames(col_names) |>
              dplyr::mutate(
                eggnog_id = id, go_type = go_type, .before = 1,
                seq_count = as.numeric(seq_count),
                frequency = as.numeric(frequency),
                annotation_count = as.numeric(annotation_count)
              )
          } else {
            tibble::tibble(
              eggnog_id = id, go_type = go_type,
              go_id = NA_character_, go_description = NA_character_,
              evidence_codes = NA_character_, seq_count = NA_real_,
              frequency = NA_real_, annotation_count = NA_real_
            )
          }
        }))

        results <- dplyr::bind_rows(results, go_df)
        cat("GO terms downloaded for Eggnog ID:", id, "\n")
      } else {
        # No GO terms found, add empty row
        results <- dplyr::bind_rows(results, tibble(
          eggnog_id = id, go_type = NA_character_, go_id = NA_character_,
          go_description = NA_character_, evidence_codes = NA_character_,
          seq_count = NA_real_, frequency = NA_real_, annotation_count = NA_real_
        ))
        cat("Warning: No GO terms found for Eggnog ID:", id, "\n")
      }
    } else {
      # Failed API request, add empty row
      results <- bind_rows(results, tibble(
        eggnog_id = id, go_type = NA_character_, go_id = NA_character_,
        go_description = NA_character_, evidence_codes = NA_character_,
        seq_count = NA_real_, frequency = NA_real_, annotation_count = NA_real_
      ))
      cat("Error: Request failed for Eggnog ID:", id, "\n")
    }
  }
  out <- tibble::as_tibble(results)
}
