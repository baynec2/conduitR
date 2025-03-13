plot_sunburst <- function(data,
                          taxonomic_columns = c("superkingdom", "kingdom",
                                                "phylum", "class", "order",
                                                "family", "genus", "species")){
  # Putting the data in a format that will work with sunburst plot.
  hierarchy = create_hierarchy_taxa_count(data) |>
    dplyr::filter(count > 5)
# Creating the plot
   plot =  plotly::plot_ly(data = hierarchy,
            type = "sunburst",
            labels = ~labels,
            parents = ~parent,
            values = ~count,
            branchvalues = "total")

   return(plot)

}

