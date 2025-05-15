#' plot_colData_distribution'
#'
#' plot the distribution of a variable present in the colData.
#'
#' @param conduit_obj
#' @param x_axis
#'
#' @returns
#' @export
#'
#' @examples
plot_colData_distribution = function(conduit_obj,
                                     x_axis){

  p1 =  SummarizedExperiment::colData(slot(conduit_obj,"QFeatures")) |>
    as.data.frame() |>
    ggplot2::ggplot(ggplot2::aes(.data[[x_axis]]))+
    ggplot2::geom_histogram(stat = "count")

  return(p1)

}
