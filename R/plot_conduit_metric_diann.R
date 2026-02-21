#' Plot DIA-NN Sample-Level Metrics from a Conduit Object
#'
#' Plots a chosen metric from the DIA-NN stats stored in the conduit object
#' (e.g. MS1 signal, retention time) by sample (File.Name). Points are
#' connected by a line to show trend across files.
#'
#' @param conduit_obj A conduit object with a `metrics` slot containing
#'   `diann_stats` (tibble with sample-level metrics).
#' @param column_name Character. Name of the column in `diann_stats` to plot
#'   on the y-axis (default: `"MS1.Signal"`).
#'
#' @return A ggplot object: x = File.Name, y = selected metric, with points
#'   and a line (group = 1); x-axis labels rotated 90 degrees.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot_conduit_metric_diann(conduit_obj, column_name = "MS1.Signal")
#' plot_conduit_metric_diann(conduit_obj, column_name = "RT.Mean")
#' }
plot_conduit_metric_diann = function(conduit_obj,column_name = "MS1.Signal"){

  diann_metrics <- slot(conduit_obj, "metrics")$diann_stats

  p1 = diann_metrics |>
    ggplot2::ggplot(ggplot2::aes(x = File.Name, y = !!rlang::sym(column_name)))+
    ggplot2::geom_point()+
    ggplot2::geom_line(ggplot2::aes(group=1))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  return(p1)
}
