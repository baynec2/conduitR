#' Create Taxonomic Tree Visualization
#'
#' Generates a hierarchical tree visualization of taxonomic relationships using
#' the metacoder package. The tree shows the relationships between different
#' taxonomic levels, with node sizes and colors representing various metrics.
#'
#' @param taxonomy_data A tibble containing taxonomic information. Must have columns:
#'   \itemize{
#'     \item organism_type: Type of organism (e.g., "Bacteria", "Archaea")
#'     \item superkingdom: Taxonomic superkingdom
#'     \item kingdom: Taxonomic kingdom
#'     \item phylum: Taxonomic phylum
#'     \item class: Taxonomic class
#'     \item order: Taxonomic order
#'     \item family: Taxonomic family
#'     \item genus: Taxonomic genus
#'     \item species: Taxonomic species
#'   }
#' @param filter_taxa_rank Character string specifying which taxonomic rank to filter by
#'   (default: "species"). The tree will show all taxa at or above this rank.
#' @param node_size Character string specifying what variable to use for node sizing
#'   (default: "n_obs"). Currently only supports "n_obs" (number of observations)
#'   or NULL (uniform node sizes).
#' @param ... Additional arguments passed to metacoder::heat_tree() for customizing
#'   the appearance of the tree.
#'
#' @return A metacoder heat tree plot object with:
#'   \itemize{
#'     \item Hierarchical tree structure showing taxonomic relationships
#'     \item Node sizes representing the specified metric
#'     \item Node labels showing taxonomic names
#'     \item Color gradients indicating node values
#'     \item Interactive hover information
#'   }
#'
#' @export
#'
#' @examples
#' # Create a basic taxonomic tree:
#' # plot_taxa_tree(taxonomy_data)
#' 
#' # Filter to show only phylum level and above:
#' # plot_taxa_tree(taxonomy_data, filter_taxa_rank = "phylum")
#' 
#' # Create a tree with uniform node sizes:
#' # plot_taxa_tree(taxonomy_data, node_size = NULL)
#' 
#' # Customize the tree appearance:
#' # plot_taxa_tree(taxonomy_data,
#' #   node_label_size_range = c(0.01, 0.05),
#' #   make_node_legend = TRUE
#' # )
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

