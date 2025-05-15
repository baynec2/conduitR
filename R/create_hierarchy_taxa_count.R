#' create_hierachy_taxa_count
#'
#' @param data
#' @param levels
#'
#' @returns
#' @export
#'
#' @examples
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
