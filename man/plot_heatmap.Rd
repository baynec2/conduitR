% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmap.R
\name{plot_heatmap}
\alias{plot_heatmap}
\title{Create Static Heatmap with Annotations}
\usage{
plot_heatmap(
  qf,
  assay_name,
  col_color_variables = NULL,
  row_color_variables = NULL,
  scale = TRUE,
  ...
)
}
\arguments{
\item{qf}{A QFeatures object containing the data to visualize}

\item{assay_name}{Character string specifying which assay to use for the heatmap}

\item{col_color_variables}{Optional character vector specifying which columns from
colData to use for column annotations (e.g., c("group", "batch"))}

\item{row_color_variables}{Optional character vector specifying which columns from
rowData to use for row annotations (e.g., c("protein_class", "pathway"))}

\item{scale}{Logical indicating whether to scale the data (default: TRUE)}

\item{...}{Additional arguments passed to sechm::sechm()}
}
\value{
A ComplexHeatmap object containing:
\itemize{
\item Main heatmap showing the data matrix
\item Column annotations using viridis color palette
\item Row annotations using magma color palette
\item Hierarchical clustering of both rows and columns
\item Color scale from viridis palette
}
}
\description{
Generates a static heatmap visualization using the sechm package,
with options for row and column annotations and hierarchical clustering.
This function provides a simpler alternative to plot_heatmaply for
non-interactive visualizations.
}
\examples{
# Basic heatmap:
# plot_heatmap(qfeatures_obj, "protein")

# With sample annotations:
# plot_heatmap(qfeatures_obj, "protein",
#             col_color_variables = c("group", "batch"))

# With both sample and feature annotations:
# plot_heatmap(qfeatures_obj, "protein",
#             col_color_variables = c("treatment", "replicate"),
#             row_color_variables = c("protein_class"))

# Without scaling:
# plot_heatmap(qfeatures_obj, "protein", scale = FALSE)
}
