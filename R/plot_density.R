#' Plot Density Distributions of Transformed Data
#'
#' Creates density plots showing the distribution of log-transformed and imputed data
#' across samples or groups. This is useful for assessing the effect of data
#' transformation and imputation on the distribution of values.
#'
#' @param qf_norm A QFeatures object containing normalized and transformed data
#' @param assay_name Character string specifying the base assay name (without transformation suffixes)
#' @param log_base Numeric value specifying the base for log transformation (default: 2)
#' @param color Character string specifying which column to use for coloring the density curves
#'   (typically a sample or group identifier)
#'
#' @return A ggplot object containing:
#'   \itemize{
#'     \item Density curves for each sample/group
#'     \item Separate facets for log-transformed and imputed data
#'     \item X-axis showing log-transformed intensity values
#'     \item Y-axis showing density
#'   }
#'
#' @export
#'
#' @examples
#' # Basic density plot colored by sample:
#' # plot_density(qfeatures_obj, "protein", color = "sample")
#' 
#' # Using log10 transformation:
#' # plot_density(qfeatures_obj, "protein", log_base = 10, color = "group")
#' 
#' # The function will automatically look for assays named:
#' # - {assay_name}_log{log_base}
#' # - {assay_name}_log{log_base}_imputed
plot_density <- function(qf_norm,
                         assay_name,
                         log_base = 2,
                         color) {
  # Log transformed values
  log <- tidy_conduit(qf_norm, paste0(assay_name, "_log", log_base)) |>
    dplyr::mutate(conduit_transformation = paste0("log", log_base))

  # Log transformed and imputed values
  imputed <- tidy_conduit(qf_norm, paste0(assay_name, "_log", log_base, "_imputed")) |>
    dplyr::mutate(conduit_transformation = paste0("log", log_base, "_imputed"))

  # Combining
  all <- dplyr::bind_rows(log, imputed)


  p1 <- all |>
    ggplot2::ggplot(ggplot2::aes(value,
      color = .data[[color]]
    )) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(~conduit_transformation) +
    ggplot2::xlab(paste0("log", log_base, "(intensity)"))

  p1
}
