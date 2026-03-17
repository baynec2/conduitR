#' Plot Data on KEGG Pathway Map
#'
#' Visualizes protein or gene expression data on a KEGG pathway map using the
#' ggkegg package.  Two pathway types are supported:
#'
#' \describe{
#'   \item{KO pathways}{IDs with prefix \code{ko} (e.g. \code{"ko04610"}) or
#'     any prefix that will be converted via \code{\link{keggpath_to_ko}}.
#'     Nodes carry KO IDs (\code{"ko:K01320"}); the function joins against
#'     \code{kegg_col = "kegg_orthology"} by default.  KO pathways are
#'     accessible via the public KEGG REST API and work with the
#'     \code{kegg_orthology} column produced by the conduitR annotation
#'     pipeline.}
#'   \item{Organism-specific pathways}{IDs with a three-letter organism prefix
#'     (e.g. \code{"hsa04610"}).  Nodes carry gene IDs (\code{"hsa:1234"});
#'     the function joins against \code{kegg_col = "xref_kegg"} by default.
#'     Requires \code{xref_kegg} to be present in the rowData — add it with
#'     \code{\link{add_kegg_ids_to_qf}}.}
#' }
#'
#' @param stats_results A data frame of statistics results (e.g. the
#'   \code{top_table} from \code{\link{perform_limma_analysis}}).
#' @param kegg_pathway_id Character string specifying the KEGG pathway ID.
#'   Organism-specific IDs (e.g. \code{"hsa04610"}) are automatically converted
#'   to their KO equivalent (\code{"ko04610"}) when \code{use_ko = TRUE}.
#' @param kegg_col Name of the list-column in \code{stats_results} to join
#'   against pathway node names.  Defaults to \code{"kegg_orthology"} for KO
#'   pathways and \code{"xref_kegg"} for organism-specific pathways.  Set
#'   explicitly to override.
#' @param use_ko Logical.  When \code{TRUE} (the default), the pathway ID is
#'   converted to a KO pathway ID and nodes are joined on KO IDs.  Set to
#'   \code{FALSE} to use the pathway ID and join column as provided.
#' @param fill_by The metric to use for node fill colour. Defaults to
#'   \code{logFC}.
#' @export
#'
#' @examples
#' # Note: KEGG pathway IDs can be found at https://www.genome.jp/kegg/pathway.html
plot_kegg_pathway <- function(stats_results,
                              kegg_pathway_id,
                              kegg_col = NULL,
                              use_ko   = TRUE,
                              fill_by  = logFC) {

  # ── Resolve pathway ID and join column ──────────────────────────────────────
  if (use_ko) {
    kegg_pathway_id <- keggpath_to_ko(kegg_pathway_id)
    if (is.null(kegg_col)) kegg_col <- "kegg_orthology"
  } else {
    if (is.null(kegg_col)) kegg_col <- "xref_kegg"
  }

  if (!kegg_col %in% colnames(stats_results)) {
    stop(
      "Column '", kegg_col, "' not found in stats_results. ",
      if (kegg_col == "xref_kegg") {
        "Add KEGG gene IDs to the QFeatures rowData with add_kegg_ids_to_qf() before running perform_limma_analysis()."
      } else {
        "Ensure the annotation pipeline has been run and the column is present in the QFeatures rowData."
      }
    )
  }

  # ── Download and parse pathway graph ────────────────────────────────────────
  graph    <- ggkegg::pathway(kegg_pathway_id, use_cache = TRUE)
  node_df  <- data.frame(tidygraph::activate(graph, "nodes"))
  edges_df <- data.frame(tidygraph::activate(graph, "edges"))

  # For KO pathways, node names are "ko:K01320" — strip the "ko:" prefix so
  # they match the bare K-numbers stored in kegg_orthology.
  if (use_ko) {
    node_df$name <- sub("^ko:", "", node_df$name)
  }

  # ── Join expression data onto nodes ─────────────────────────────────────────
  # Some Protein.Groups map to multiple KEGG IDs; unnest so each gets one row.
  stats_results <- tidyr::unnest(stats_results, tidyr::all_of(kegg_col))

  new_node_df <- dplyr::left_join(
    node_df,
    stats_results,
    by = stats::setNames(kegg_col, "name")
  )

  # ── Build and plot ───────────────────────────────────────────────────────────
  new_graph <- tidygraph::tbl_graph(nodes = new_node_df, edges_df)

  node_type <- if (use_ko) "ortholog" else "gene"

  plot <- new_graph |>
    ggraph::ggraph(layout = "manual", x = x, y = y) +
      ggkegg::geom_node_rect(ggplot2::aes(
        filter = type == node_type,
        fill   = {{fill_by}}
      )) +
      ggkegg::overlay_raw_map() +
      ggplot2::scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0
      ) +
      ggplot2::theme_void()

  return(plot)
}
