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
#'   the data to use for coloring (e.g., expression values, fold changes). If NULL,
#'   will attempt to extract data from the QFeatures object.
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
#' # plot_kegg_pathway(qfeatures_obj, "protein", "00010")
#' 
#' # Plot RNA degradation pathway with custom IDs:
#' # plot_kegg_pathway(qfeatures_obj, "log2_imputed", "03018",
#' #                  ids_to_color = c(K04077 = 1, K04078 = -1))
#' 
#' # Note: KEGG pathway IDs can be found at https://www.genome.jp/kegg/pathway.html
plot_kegg_pathway <- function(qf,
                             assay_name,
                             kegg_pathway_id,
                             ids_to_color = NULL) {
  
  # Input validation
  if (!inherits(qf, "QFeatures")) {
    stop("'qf' must be a QFeatures object")
  }
  
  if (!assay_name %in% names(qf)) {
    stop("'assay_name' not found in QFeatures object")
  }
  
  # Extract data from QFeatures if ids_to_color not provided
  if (is.null(ids_to_color)) {
    assay_data <- assay(qf, assay_name)
    # Assuming row names are KEGG IDs - you might need to adjust this
    ids_to_color <- setNames(as.numeric(assay_data[,1]), rownames(assay_data))
  }
  
  if (!is.numeric(ids_to_color) || is.null(names(ids_to_color))) {
    stop("'ids_to_color' must be a named numeric vector")
  }
  
  # Create a temporary directory for this specific operation
  temp_dir <- tempfile("pathview_")
  dir.create(temp_dir)
  on.exit({
    # Clean up temporary files and directory
    unlink(temp_dir, recursive = TRUE, force = TRUE)
  })
  
  # Run pathview and generate image
  tryCatch({
    pathview::pathview(
      gene.data = ids_to_color,
      pathway.id = kegg_pathway_id,
      species = "ko",
      out.suffix = "pathview",
      kegg.dir = temp_dir
    )
    
    # Find the PNG file
    png_files <- list.files(temp_dir, pattern = ".*pathview.png$", full.names = TRUE)
    if (length(png_files) == 0) {
      stop("No pathway image was generated")
    }
    
    # Load and plot image
    img <- png::readPNG(png_files[1])
    grid::grid.raster(img)
    
  }, error = function(e) {
    stop("Error generating pathway view: ", e$message)
  })
  
  invisible(NULL)
}
