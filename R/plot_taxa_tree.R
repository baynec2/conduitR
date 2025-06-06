#' Create Taxonomic Tree Visualization
#'
#' Generates a hierarchical tree visualization of taxonomic relationships using
#' the metacoder package. The tree shows the relationships between different
#' taxonomic levels, with node sizes and colors representing various metrics.
#'
#' @param taxonomy_data A tibble containing taxonomic information. Must have columns:
#'   \itemize{
#'     \item organism_type: Type of organism (e.g., "Bacteria", "Archaea")
#'     \item domain: Taxonomic domain
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
                           node_color = "phylum",
                           ...) {
  if (nrow(taxonomy_data) == 1) {
    cat("There is only one taxa in the data.")
    taxonomy_data <- tidyr::pivot_longer(taxonomy_data,
      cols = c(
        "organism_type", "domain", "kingdom", "phylum", "class", "order",
        "family", "genus", "species"
      ),
      names_to = "taxon_rank",
      values_to = "taxon"
    ) |>
      dplyr::mutate(taxon_rank = factor(taxon_rank,
        levels = c(
          "organism_type", "domain",
          "kingdom", "phylum", "class",
          "order", "family", "genus",
          "species"
        )
      )) |>
      dplyr::arrange(taxon_rank)

    p1 <- ggplot2::ggplot(taxonomy_data, ggplot2::aes(x = 1:9, y = 0.5)) +
      ggplot2::geom_point(
        shape = 21, size = 15,
        fill = "lightblue", color = "black"
      ) +
      ggplot2::geom_text(ggplot2::aes(label = taxon), size = 4) +
      ggplot2::xlim(0, 9) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_void()

    return(p1)
  } else {
    # Adding a string with all levels of taxonomy concatenated
    taxonomy <- taxonomy_data |>
      dplyr::mutate(taxonomy = paste0(
        "organism_type__", organism_type, ";",
        "domain__", domain, ";",
        "kingdom__", kingdom, ";",
        "phylum__", phylum, ";",
        "class__", class, ";",
        "order__", order, ";",
        "family__", family, ";",
        "genus__", genus, ";",
        "species__", species
      ))

    taxmap <- metacoder::parse_tax_data(taxonomy,
      class_cols = "taxonomy",
      class_sep = ";", # The character used to separate taxa
      class_key = c(
        tax_rank = "taxon_rank", # A key describing each regex capture group
        tax_name = "taxon_name"
      ),
      class_regex = "^(.+)__(.+)$"
    ) # Regex identifying where the data for each taxon is


    ft <- taxmap |>
      metacoder::filter_taxa(taxon_ranks == filter_taxa_rank, supertaxa = TRUE)

    # Defining labels
    lab <- taxonomy_data |>
      dplyr::mutate(
        node_color_lab = as.factor(!!dplyr::sym(node_color)),
        node_color = as.numeric(node_color_lab)
      ) |>
      dplyr::select(species, node_color_lab, node_color)

    # mapping node color
    node_color_df <- ft$data$class_data |>
      dplyr::left_join(lab, by = c("tax_name" = "species")) |>
      # There will always be NA colored nodes. Changing to 0
      dplyr::mutate(node_color = dplyr::case_when(
        is.na(node_color) ~ 0,
        TRUE ~ node_color
      )) |>
      # Adding 1 to NA (0) nodes so we can select color in vector later.
      dplyr::mutate(node_color = node_color + 1)

    # Assigning node color to ft object
    ft$data$class_data <- node_color_df

    # Function to color heat tree
    heat_tree_color <- function(number_vector) {
      out <- viridis::magma(length(unique(number_vector)))
      return(out)
    }

    # Generating color pallete
    pal <- heat_tree_color(ft$data$class_data$node_color)

    # Defining plot
    if (is.null(node_size)) {
      p1 <- metacoder::heat_tree(ft,
        node_label = taxon_names,
        make_node_legend = FALSE,
        node_label_size_range = c(0.005,0.005),
        node_color = node_color,
        node_color_range = pal,
        node_label_max = 5000,
        repel_labels = TRUE,
        node_label_color = "#1E90FF",
        ...
      )
    } else {
      if (node_size == "n_obs") {
        p1 <- metacoder::heat_tree(ft,
          node_color = node_color,
          node_label = taxon_names,
          node_size = n_obs,
          node_color_range = pal,
          node_label_max = 5000,
          repel_labels = TRUE,
          node_label_color = "#1E90FF",
          make_node_legend = FALSE,
          ...
        )
      } else {
        cat("Node size must be n_obj or NULL \n")
      }
    }

    # Making the legend
    legend_df <- node_color_df |>
      dplyr::filter(!is.na(node_color_lab)) |>
      dplyr::select(node_color_lab, node_color) |>
      dplyr::distinct()

    # dealing with NAs
    zero <- data.frame(
      node_color_lab = "NA",
      node_color = 1
    )

    # Assembling into final
    final_legend_df <- dplyr::bind_rows(legend_df, zero) |>
      dplyr::arrange(node_color) |>
      dplyr::mutate(color = pal[node_color])

    # building a legend plot
    legend_plot <- ggplot2::ggplot(final_legend_df, ggplot2::aes(x = 1, y = 1, color = node_color_lab)) +
      ggplot2::geom_point(size = 5) +
      ggplot2::scale_color_manual(values = setNames(final_legend_df$color, final_legend_df$node_color_lab)) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right")

    # Extracting the legend
    legend_only <- cowplot::get_plot_component(legend_plot, "guide-box-right", return_all = TRUE)

    # Adding the legend to plot
    p2 <- cowplot::plot_grid(p1, legend_only, ncol = 2, rel_widths = c(1, 0.3))

    # Returning final plot
    return(p2)
  }
}
