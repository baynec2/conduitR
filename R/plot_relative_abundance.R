#' Plot Relative Abundance of Features Across Samples
#'
#' Creates a stacked bar plot showing the relative abundance of features (e.g., proteins or species)
#' across all samples in a QFeatures object. Each bar represents a sample, with the
#' height of each segment representing the relative abundance of a feature. The function
#' automatically calculates relative abundances by normalizing each sample to 100%,
#' making it easy to compare the composition of features across samples.
#'
#' @param qf A QFeatures object containing the abundance data. The assay should contain
#'   non-negative abundance values (e.g., protein intensities or species counts).
#' @param assay_name Character string specifying which assay to use for plotting.
#'   The assay should contain the raw abundance values that will be converted to
#'   relative abundances.
#' @param facet_formula Optional formula for faceting the plot (e.g., ~group or ~treatment + timepoint).
#'   This allows for comparing relative abundances across different experimental groups
#'   or conditions.
#'
#' @return A ggplot object containing a stacked bar plot with:
#'   \itemize{
#'     \item X-axis showing samples (labels hidden to prevent overcrowding)
#'     \item Y-axis showing relative abundance (0-100%)
#'     \item Stacked bars colored by feature (e.g., protein or species)
#'     \item Optional faceting based on the provided formula
#'     \item Legend showing feature names
#'     \item Free x-axis scales within facets to optimize space usage
#'   }
#'   The plot is designed to show the proportional composition of features
#'   within each sample, making it easy to identify dominant features and
#'   compare patterns across samples or groups.
#'
#' @export
#'
#' @examples
#' # Basic relative abundance plot:
#' # plot_relative_abundance(qfeatures_obj, "protein")
#'
#' # With faceting by experimental group:
#' # plot_relative_abundance(qfeatures_obj, "protein", ~group)
#'
#' # With faceting by multiple variables:
#' # plot_relative_abundance(qfeatures_obj, "protein", ~treatment + timepoint)
#'
#' # The resulting plot can be customized further:
#' # p <- plot_relative_abundance(qfeatures_obj, "protein", ~group)
#' # p + ggplot2::theme_minimal() +
#' #     ggplot2::scale_fill_viridis_d() +
#' #     ggplot2::labs(title = "Protein Relative Abundance by Group")
#'
#' @note
#' The input data should contain non-negative abundance values. For large datasets
#' (e.g., >1000 features), consider the following:
#' \itemize{
#'   \item Filter to the most abundant features before plotting to improve readability
#'   \item Use facet_formula to split the visualization into smaller, more manageable plots
#'   \item Consider using a subset of samples if the dataset is very large
#' }
#' The function uses tidy_conduit() internally to reshape the data, which may require
#' significant memory for large datasets. For optimal performance, pre-filter the data
#' to include only the features of interest.
#'
#' @seealso \code{\link[tidy_conduit]{tidy_conduit}} for the data reshaping function used internally
plot_relative_abundance <- function(qf,
                                    assay_name,
                                    facet_formula){

  tidy_qf <- conduitR::tidy_conduit(qf,assay_name)

  p1 = tidy_qf |>
    ggplot2::ggplot(ggplot2::aes(file,
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
