#' plot_taxa_tree
#'
#' plot a taxa tree using taxonomy data.
#'
#'
#' @param taxonomy_data taxonomy data as a tibble. Each observation is a row,
#' columns are.
#' @param filter_taxa_rank what rank to filter taxa by.
#' @param node_size what variable to color node size by.
#'
#' @returns
#' @export
#'
#' @examples
plot_taxa_tree <- function(taxonomy_data,
                           filter_taxa_rank = "species",
                           node_size = "n_obs",
                           ...) {
  # Adding a string with all levels of taxonomy concatenated
  taxonomy <- taxonomy_data |>
    dplyr::mutate(taxonomy = paste0(
      "organism_type__", organism_type, ";",
      "superkingdom__", superkingdom, ";",
      "kindom__", kingdom, ";",
      "phylum__", phylum, ";",
      "class__", class, ";",
      "order__", order, ";",
      "family__", family, ";",
      "genus__", genus, ";",
      "species__", species
    ))

  taxmap <- metacoder::parse_tax_data(taxonomy,
    class_cols = "taxonomy",
    class_sep = ";", # The character used to separate taxa in the classification
    class_key = c(
      tax_rank = "taxon_rank", # A key describing each regex capture group
      tax_name = "taxon_name"
    ),
    class_regex = "^(.+)__(.+)$"
  ) # Regex identifying where the data for each taxon is


  ft <- taxmap |>
    metacoder::filter_taxa(taxon_ranks == filter_taxa_rank, supertaxa = TRUE)

  # Defining plot
  if (is.null(node_size)) {
    p1 <- heat_tree(ft,
      node_label = taxon_names,
      make_node_legend = FALSE,
      node_label_size_range = c(0.005, 0.05),
      ...
    )
  } else {
    if (node_size == "n_obs") {
      p1 <- metacoder::heat_tree(ft,
        node_label = taxon_names,
        node_size = n_obs,
        ...
      )
    } else {
      cat("Node size must be n_obj or NULL")
    }
  }
  return(p1)
}

