#' Create Interactive Sunburst Plot of Taxonomic Distribution
#'
#' Generates an interactive sunburst plot showing the hierarchical distribution
#' of proteins across taxonomic levels. The plot allows users to explore the
#' taxonomic composition of the dataset by clicking on segments to zoom in
#' and see more detailed information.
#'
#' @param taxonomy A tibble containing taxonomic information for detected proteins,
#'   with columns for organism_type, domain, kingdom, phylum, class, order,
#'   family, genus, and species (e.g. \code{conduit_obj@@taxonomy}).
#'
#' @return An interactive plotly sunburst plot object with the following features:
#'   \itemize{
#'     \item Hierarchical segments representing taxonomic levels
#'     \item Segment sizes proportional to protein counts
#'     \item Interactive zooming and clicking functionality
#'     \item Hover information showing taxonomic details and counts
#'     \item Color-coded segments for different taxonomic levels
#'   }
#'
#' @export
#'
#' @examples
#' # Create a basic sunburst plot:
#' # plot_sunburst(conduit_obj)
#' 
#' # Create a plot counting only uniquely identified proteins:
#' # plot_sunburst(conduit_obj, type = "one_protein_in_group")
#' 
#' # Customize the plot appearance:
#' # plot_sunburst(conduit_obj,
#' #   colors = c("red", "blue", "green"),
#' #   hoverinfo = "text",
#' #   textinfo = "label+value"
#' # )
#' 
#' # The plot can be used alongside other taxonomic visualizations:
#' # p1 <- plot_sunburst(conduit_obj)
#' # p2 <- plot_percent_detected_taxa_tree(conduit_obj)
#' # plotly::subplot(p1, p2)
plot_sunburst = function(taxonomy){

  #Summarizing taxonomy by some metric.
  sum = taxonomy |>
    dplyr::group_by(organism_type,domain,kingdom,phylum,class,order,family,genus,species) |>
    dplyr::summarise(size = dplyr::n(),.groups = "drop") |>
    dplyr::mutate(domain = as.character(domain)) |>
    tidyr::replace_na(list(
      domain = "Unknown",
      kingdom = "Unknown",
      phylum = "Unknown",
      class = "Unknown",
      order = "Unknown",
      family = "Unknown",
      genus = "Unknown",
      species = "Unknown"
    )) |> dplyr::ungroup() |>
    as.data.frame()
  # Converting into nested json object
  tree <- d3r::d3_nest(sum, value_cols = "size")
  # Plotting with sunburstR
  plot <- sunburstR::sunburst(tree, width="100%", height=400,legend = FALSE)

  return(plot)
}
