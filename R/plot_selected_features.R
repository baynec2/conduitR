#' Plot Selected Features Across Samples
#'
#' Creates a scatter plot showing the expression/abundance values of selected features
#' (e.g., proteins) across samples, with options for customizing the visualization
#' through faceting, coloring, and point shapes.
#'
#' @param qf A QFeatures object containing the data to plot
#' @param assay_name Character string specifying which assay to use for plotting
#' @param features Character vector of feature IDs to plot (e.g., protein IDs)
#' @param x_axis Character string specifying the column name to use for the x-axis
#'   (typically a sample or group identifier)
#' @param facet_formula Optional formula for faceting the plot (e.g., ~group)
#' @param color_by Optional character string specifying a column to use for point colors
#' @param shape Optional character string specifying a column to use for point shapes
#'
#' @return A ggplot object containing a scatter plot with:
#'   \itemize{
#'     \item X-axis showing the specified variable
#'     \item Y-axis showing feature values
#'     \item Points colored and/or shaped by specified variables
#'     \item Optional faceting based on the provided formula
#'   }
#'
#' @export
#'
#' @examples
#' # Basic plot of selected proteins:
#' # plot_selected_features(qfeatures_obj, "protein", 
#' #                       features = c("P12345", "P67890"),
#' #                       x_axis = "sample")
#' 
#' # With grouping and faceting:
#' # plot_selected_features(qfeatures_obj, "protein",
#' #                       features = c("P12345", "P67890"),
#' #                       x_axis = "timepoint",
#' #                       color_by = "group",
#' #                       facet_formula = ~treatment)
#' 
#' # With custom point shapes:
#' # plot_selected_features(qfeatures_obj, "protein",
#' #                       features = c("P12345", "P67890"),
#' #                       x_axis = "sample",
#' #                       shape = "replicate")
plot_selected_features = function(qf,
                                  assay_name,
                                  features,
                                  x_axis,
                                  facet_formula = NULL,
                                  color_by = NULL,
                                  shape = NULL){

  # Getting tidy format
  tidy = conduitR::tidy_conduit(qf, assay_name) |>
    dplyr::filter(rowid %in% features)

  # Dealing with weird shiny "" input problems
  if (is.null(facet_formula) || identical(facet_formula, "")) facet_formula <- NULL
  if (is.null(color_by) || identical(color_by, "")) color_by <- NULL
  if (is.null(shape) || identical(shape, "")) color_by <- NULL

  # Initialize ggplot
  p1 = tidy |>
    ggplot2::ggplot(ggplot2::aes(x = .data[[x_axis]], y = value))

  # Add geom_point with color/shape logic
  aes_args = list()
  if (!is.null(color_by)) aes_args$color = rlang::sym(color_by)
  if (!is.null(shape)) aes_args$shape = rlang::sym(shape)

  p1 = p1 + ggplot2::geom_point(mapping = do.call(ggplot2::aes, aes_args))

  # Add faceting if provided
  if (!is.null(facet_formula)) {
    p1 = p1 + ggplot2::facet_wrap(as.formula(facet_formula))
  }

  return(p1)
}
