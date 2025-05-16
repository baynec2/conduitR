#' Create Hierarchical Taxonomy Count Structure
#'
#' Generates a hierarchical data structure for taxonomic counts, creating parent-child
#' relationships between different taxonomic levels. This function is particularly
#' useful for preparing data for hierarchical visualizations like sunburst plots or
#' taxonomic trees.
#'
#' @param data A data frame or tibble containing taxonomic information. Must have
#'   columns corresponding to the taxonomic levels specified in the `levels` parameter.
#'   Each row should represent a unique taxon at the most specific level (e.g., species).
#' @param levels A character vector specifying the taxonomic levels to include in the
#'   hierarchy, ordered from most general to most specific (default: c("domain",
#'   "kingdom", "phylum", "class", "order", "family", "genus", "species")). The
#'   columns in `data` must match these level names.
#'
#' @return A tibble containing:
#'   \itemize{
#'     \item parent: The taxonomic name of the parent level
#'     \item labels: The taxonomic name of the current level
#'     \item ids: Unique identifiers for each taxon (same as labels)
#'     \item count: Number of entries for each taxon
#'   }
#'   The output is structured hierarchically, with:
#'   \itemize{
#'     \item A root level (empty parent) containing counts of unique taxa at the
#'       most general level
#'     \item Subsequent levels showing parent-child relationships between taxa
#'     \item Counts representing the number of entries at each taxonomic level
#'   }
#'
#' @export
#'
#' @examples
#' # Create hierarchy from taxonomy data:
#' # hierarchy <- create_hierarchy_taxa_count(taxonomy_data)
#' 
#' # Use custom taxonomic levels:
#' # hierarchy <- create_hierarchy_taxa_count(taxonomy_data,
#' #   levels = c("phylum", "class", "order", "family"))
#' 
#' # The resulting hierarchy can be used with visualization functions:
#' # plot_sunburst(hierarchy)
#' # plot_taxa_tree(hierarchy)
#'
#' @note
#' This function is particularly useful for:
#' \itemize{
#'   \item Preparing data for hierarchical visualizations
#'   \item Analyzing taxonomic distributions
#'   \item Creating interactive taxonomic plots
#' }
#' The function assumes that:
#' \itemize{
#'   \item Taxonomic levels are ordered from general to specific
#'   \item Each taxon at a specific level belongs to exactly one parent taxon
#'   \item The input data contains all necessary taxonomic columns
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[plot_sunburst]{plot_sunburst}} for creating interactive
#'     sunburst plots from the hierarchy
#'   \item \code{\link[plot_taxa_tree]{plot_taxa_tree}} for creating taxonomic
#'     tree visualizations
#'   \item \code{\link[count_proteins_per_taxon]{count_proteins_per_taxon}} for
#'     counting proteins at each taxonomic level
#' }
create_hierarchy_taxa_count <- function(data, levels = c("domain",
                                                         "kingdom",
                                                         "phylum",
                                                         "class",
                                                         "order",
                                                         "family",
                                                         "genus",
                                                         "species")) {
  hierarchy_list <- list()

  # Loop through levels to create parent-child relationships
  for (i in seq_along(levels)[-length(levels)]) {
    parent_col <- levels[i]
    child_col <- levels[i + 1]

    df <- data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(parent_col, child_col)))) |>
      dplyr::summarize(count = dplyr::n(), .groups = "drop") |>
      dplyr::rename(parent = !!parent_col, labels = !!child_col) |>
      dplyr::mutate(ids = labels)

    hierarchy_list[[i]] <- df
  }

  # Create top-level (root) aggregation - count of unique taxa in the first level
  root_level <- data |>
    dplyr::group_by(!!rlang::sym(levels[1])) |>
    dplyr::summarize(count = dplyr::n_distinct(!!rlang::sym(levels[1])), .groups = "drop") |>
    dplyr::rename(labels = !!rlang::sym(levels[1])) |>
    dplyr::mutate(parent = "", ids = labels)

  hierarchy_list <- base::append(list(root_level), hierarchy_list)

  # Combine all levels into one data frame
  out = dplyr::bind_rows(hierarchy_list)

  return(out)
}
