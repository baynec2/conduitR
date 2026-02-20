#' Plot Diann Sample Level Metrics for Conduit Object
#'
#' @param conduit_obj conduit object
#' @param column_name column name to plot
#'
#' @returns
#' @export
#'
#' @examples
plot_conduit_metric_diann = function(conduit_obj,column_name = "MS1.Signal"){

  diann_metrics <- slot(conduit_obj, "metrics")$diann_stats

  p1 = diann_metrics |>
    ggplot2::ggplot(ggplot2::aes(x = File.Name, y = !!rlang::sym(column_name)))+
    ggplot2::geom_point()+
    ggplot2::geom_line(ggplot2::aes(group=1))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  return(p1)
}
