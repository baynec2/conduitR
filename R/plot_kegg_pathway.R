#' Plot Data on KEGG Pathway Map
#'
#' Visualizes protein or gene expression data on a KEGG pathway map using the
#' ggkegg package.
#'
#' This function creates a pathway diagram where nodes (proteins/genes) are colored according
#' to their expression values or other quantitative data.
#'
#' @param stats_results some sort of statistics results (intended to be the
#' toptable produced by the perfom_limma_stats function). Must have a column
#' named xref_kegg containing the kegg id corresponding to each protein.
#' @param kegg_pathway_id Character string specifying the KEGG pathway ID to
#' plot
#' @param fill_by the metric to fill the color of the boxes by.
#' @export
#'
#' @examples
#' # Note: KEGG pathway IDs can be found at https://www.genome.jp/kegg/pathway.html
plot_kegg_pathway <- function(stats_results,
                              kegg_pathway_id,
                              fill_by = logFC) {

  # Downloading kegg pathway of interest
  graph <- ggkegg::pathway(kegg_pathway_id, use_cache=TRUE)

  # Extracting the nodes and edges from the pathway.
  node_df <- data.frame(tidygraph::activate(graph, "nodes"))
  edges_df <- data.frame(tidygraph::activate(graph,"edges"))

  # Some Protein.Groups are linked to multiple KEGG IDs.
  # Split them so each KEGG ID can be visualized individually.

  stats_results <- stats_results |>
    tidyr::unnest("xref_kegg")

  # Modifying the nodes to include stats results
  new_node_df <- dplyr::left_join(node_df,
                                  stats_results,
                                  by = c("name" = "xref_kegg"))

  #Creating a new graph with the new node_df and the old edges.
  new_graph <- tidygraph::tbl_graph(nodes = new_node_df,edges_df)

  # Plotting the keggg pathway
  plot <- new_graph |>
    ggraph::ggraph(layout="manual", x=x, y=y)+
        ggkegg::geom_node_rect(ggplot2::aes(
                      filter=type == "gene",
                      fill = {{fill_by}}))+
        ggkegg::overlay_raw_map()+
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
        ggplot2::theme_void()

  return(plot)

}
