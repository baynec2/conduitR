#' plot_kegg_pathway
#'
#' plot values in data on kegg pathway.
#'
#' @param qf q features object
#' @param assay_name assay name with IDs
#' @param kegg_pathway_id kegg pathway ID
#' @param ids_to_color ids to plot.
#'
#' @returns
#' @export
#'
#' @examples
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
