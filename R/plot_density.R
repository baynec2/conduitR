#' plot_density
#'
#' @param qf_norm
#' @param log_base
#' @param assay_name
#' @param color
#'
#' @returns
#' @export
#'
#' @examples
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
