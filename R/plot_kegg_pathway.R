#' Plot Data on KEGG Pathway Map
#'
#' Visualizes protein or gene expression data on a KEGG pathway map using the pathview package.
#' This function creates a pathway diagram where nodes (proteins/genes) are colored according
#' to their expression values or other quantitative data.
#'
#' @param qf A QFeatures object containing the data to visualize
#' @param assay_name Character string specifying which assay contains the IDs to map
#'   (default: "log2_imputed")
#' @param kegg_pathway_id Character string specifying the KEGG pathway ID to plot
#'   (default: "03018", which is RNA degradation)
#' @param ids_to_color Named numeric vector where names are KEGG IDs and values are
#'   the data to use for coloring (e.g., expression values, fold changes)
#'
#' @return Invisibly returns NULL. The function displays a pathway map where:
#'   \itemize{
#'     \item Nodes represent proteins/genes in the pathway
#'     \item Node colors indicate the values in ids_to_color
#'     \item The pathway structure shows relationships between components
#'   }
#'
#' @export
#'
#' @examples
#' # Plot a specific pathway with expression data:
#' # plot_kegg_pathway(qfeatures_obj, "protein", "00010", 
#' #                  ids_to_color = c(K00001 = 2.5, K00002 = -1.5))
#' 
#' # Plot RNA degradation pathway with custom IDs:
#' # plot_kegg_pathway(qfeatures_obj, "log2_imputed", "03018",
#' #                  ids_to_color = c(K04077 = 1, K04078 = -1))
#' 
#' # Note: KEGG pathway IDs can be found at https://www.genome.jp/kegg/pathway.html
plot_kegg_pathway <- function(qf,
                              assay_name = "log2_imputed",
                              kegg_pathway_id = "03018",
                              ids_to_color = c(K04077 = 1)) {

  # Use a temporary directory
  temp_dir <- tempdir()
  original_wd <- getwd()
  setwd(temp_dir)

  # Run pathview and generate image
  pathview::pathview(
    gene.data = ids_to_color,
    pathway.id = kegg_pathway_id,
    species = "ko"
  )

  # Find the PNG file (returns the first match)
  file <- list.files(temp_dir, pattern = ".*pathview.png$", full.names = TRUE)

  # Load and plot image
  img <- png::readPNG(file[1])
  grid::grid.raster(img)

  # Restore original working directory
  setwd(original_wd)

  invisible(NULL)

}
