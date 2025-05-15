#' plot_selected_features
#'
#' plot the selected features
#'
#' @param qf
#' @param assay_name
#' @param features
#' @param x_axis
#' @param facet_formula
#' @param color_by
#' @param shape
#'
#' @returns
#' @export
#'
#' @examples
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
