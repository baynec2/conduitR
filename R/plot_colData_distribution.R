#' Plot Distribution of Column Metadata
#'
#' Creates a histogram showing the distribution of a variable from the column metadata
#' (colData) of a conduit object. This is useful for exploring the distribution of
#' sample characteristics and experimental variables.
#'
#' @param conduit_obj A conduit object containing the QFeatures data
#' @param x_axis Character string specifying which column from colData to plot
#'
#' @return A ggplot object containing:
#'   \itemize{
#'     \item Histogram showing the distribution of the specified variable
#'     \item X-axis showing the variable values
#'     \item Y-axis showing the count/frequency
#'   }
#'
#' @export
#'
#' @examples
#' # Plot distribution of a categorical variable:
#' # plot_colData_distribution(conduit_obj, "group")
#'
#' # Plot distribution of a numeric variable:
#' # plot_colData_distribution(conduit_obj, "concentration")
#'
#' # Note: For numeric variables, you might want to use geom_density() instead
#' # of the default geom_histogram() by modifying the returned plot
plot_colData_distribution = function(conduit_obj,
                                     x_axis){

  p1 =  SummarizedExperiment::colData(slot(conduit_obj,"QFeatures")) |>
    as.data.frame() |>
    ggplot2::ggplot(ggplot2::aes(.data[[x_axis]]))+
    ggplot2::geom_histogram(stat = "count")+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  return(p1)

}
