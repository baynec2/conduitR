#' Plot Density Distributions of Transformed Data
#'
#' Creates density plots showing the distribution of log-transformed and imputed data
#' across samples or groups. This is useful for assessing the effect of data
#' transformation and imputation on the distribution of values.
#'
#' @param qf_norm A QFeatures object containing the transformed data.
#' @param log_assay Character string specifying the name of the log-transformed assay
#'   (e.g. "protein_groups_log2").
#' @param imputed_assay Character string specifying the name of the imputed assay
#'   (e.g. "protein_groups_log2_MinDet").
#' @param color Character string specifying which colData column to use for coloring
#'   the density curves.
#'
#' @return A ggplot object with density curves faceted by transformation stage.
#'
#' @export
#'
#' @examples
#' # qf <- add_log_imputed_norm_assay(qf, "protein_groups")
#' # plot_density(qf,
#' #   log_assay     = "protein_groups_log2",
#' #   imputed_assay = "protein_groups_log2_MinDet",
#' #   color         = "group")
plot_density <- function(qf_norm,
                         log_assay,
                         imputed_assay,
                         color) {
  # Log transformed values
  log <- tidy_conduit(qf_norm, log_assay) |>
    dplyr::mutate(conduit_transformation = log_assay)

  # Log transformed and imputed values
  imputed <- tidy_conduit(qf_norm, imputed_assay) |>
    dplyr::mutate(conduit_transformation = imputed_assay)

  # Combining
  all <- dplyr::bind_rows(log, imputed)


  p1 <- all |>
    ggplot2::ggplot(ggplot2::aes(value,
      color = .data[[color]]
    )) +
    ggplot2::geom_density() +
    ggplot2::facet_wrap(~conduit_transformation) +
    ggplot2::xlab("log-transformed intensity")

  p1
}
