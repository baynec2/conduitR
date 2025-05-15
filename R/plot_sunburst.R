#' plot_sunburst
#'
#' Create a sunburst plot.
#'
#' @param taxonomy tibble with taxonomy information. Must have columns
#' superkingdom,kindom,phylum,class,order,family,genus,species.
#'
#' @returns
#' @export
#'
#' @examples
plot_sunburst = function(taxonomy){

  #Summarizing taxonomy by some metric.
  sum = taxonomy |>
    dplyr::group_by(organism_type,superkingdom,kingdom,phylum,class,order,family,genus,species) |>
    dplyr::summarise(size = dplyr::n(),.groups = "drop") |>
    dplyr::mutate(superkingdom = as.character(superkingdom)) |>
    tidyr::replace_na(list(
      superkingdom = "Unknown",
      kingdom = "Unknown",
      phylum = "Unknown",
      class = "Unknown",
      order = "Unknown",
      family = "Unknown",
      genus = "Unknown",
      species = "Unknown"
    )) |> dplyr::ungroup() |>
    as.data.frame()
  # Converting into nested json object
  tree <- d3r::d3_nest(sum, value_cols = "size")
  # Plotting with sunburstR
  plot <- sunburstR::sunburst(tree, width="100%", height=400,legend = FALSE)

  return(plot)
}
