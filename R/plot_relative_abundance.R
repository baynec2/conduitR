#' plot_relative_abundance
#' Plot the relative abundance for all samples.
#'
#' @param qf q features object
#' @param assay_name assay to get relative abundance for
#' @param facet_formula formula to specify how the plot is faceted
#'
#' @returns
#' @export
#'
#' @examples
plot_relative_abundance <- function(qf,
                                    assay_name,
                                    facet_formula){

  tidy_qf <- conduitR::tidy_conduit(qf,assay_name)

  p1 = tidy_qf |>
    ggplot2::ggplot(ggplot2::aes(sample_name,
                                 value,
                                 fill = rowid))+
    ggplot2::geom_col()+
    ggplot2::facet_grid(facet_formula,scales = "free_x",space = "free_x")+
    ggplot2::theme(legend.position = "right",
                   axis.text.x = ggplot2::element_blank())

  # Add faceting if requested
  if (!is.null(facet_formula)) {
    p1 <- p1 + ggplot2::facet_wrap(facet_formula)
  }

return(p1)

}
