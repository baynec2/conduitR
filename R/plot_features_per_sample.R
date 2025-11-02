#' Title
#'
#' @param qf
#' @param assay
#'
#' @returns
#' @export
#'
#' @examples
plot_features_per_sample <- function(qf, assay = "protein_groups") {
  se <- qf[[assay]]

  # Convert to long format and binarize
  df <- assay(se) |>
    data.frame() |>
    tibble::rownames_to_column("rowname") |>
    tidyr::pivot_longer(
      cols = -rowname,
      names_to = "ID",
      values_to = "bin"
    ) |>
    dplyr::mutate(bin = ifelse(is.na(bin) | bin == 0, 0, 1))

  # Count number of samples each feature is detected in
  stat <- df |>
    dplyr::group_by(rowname) |>
    dplyr::summarize(sum = sum(bin), .groups = "drop")

  # Frequency table of detection counts
  table <- stat |>
    dplyr::count(sum) |>
    dplyr::arrange(dplyr::desc(sum))

  # Plot
  p <- ggplot2::ggplot(table, ggplot2::aes(x = sum, y = n, fill = sum)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_gradient(name = "Detected\nin n\nSamples") +
    ggplot2::ylab(paste0(assay, " count")) +
    ggplot2::xlab("Number of samples detected in")

  return(p)
}
