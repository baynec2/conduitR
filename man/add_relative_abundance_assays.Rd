% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_relative_abundance_assays.R
\name{add_relative_abundance_assays}
\alias{add_relative_abundance_assays}
\title{Add Relative Abundance Assays to QFeatures Object}
\usage{
add_relative_abundance_assays(
  qf,
  assay_names = c("domain", "kingdom", "phylum", "class", "order", "family", "genus",
    "species")
)
}
\arguments{
\item{qf}{A QFeatures object containing the abundance data. The assays should contain
non-negative abundance values (e.g., protein intensities or species counts).}

\item{assay_names}{A character vector specifying which assays to convert to relative
abundances. By default, includes all taxonomic levels: domain, kingdom, phylum,
class, order, family, genus, and species. Each assay should exist in the QFeatures
object.}
}
\value{
A QFeatures object with additional assays named "{assay_name}_rel_abundance"
for each input assay. Each new assay contains:
\itemize{
\item Values normalized to percentages (0-100\%)
\item Same dimensions as the input assay
\item NA values preserved from the input
\item Column sums equal to 100\% for each sample (excluding NA values)
}
}
\description{
Adds new assays containing relative abundance values (as percentages) to a QFeatures object
for each specified taxonomic level. This function calculates the relative abundance of
features within each sample by normalizing the raw abundance values to percentages,
making it easier to compare the composition of features across samples.
}
\note{
The input assays should contain non-negative abundance values. The function:
\itemize{
\item Preserves NA values from the input data
\item Handles each assay independently
\item Normalizes each sample separately
}
For large datasets, consider processing one assay at a time to manage memory usage.
}
\examples{
# Add relative abundance assays for all default taxonomic levels:
# qf_with_rel_abundance <- add_relative_abundance_assays(qfeatures_obj)

# Add relative abundance for specific assays:
# qf_with_rel_abundance <- add_relative_abundance_assays(qfeatures_obj,
#   assay_names = c("genus", "species"))

# The resulting assays can be used with plotting functions:
# plot_relative_abundance(qf_with_rel_abundance, "genus_rel_abundance")

# Or for downstream analysis:
# perform_limma_analysis(qf_with_rel_abundance, "species_rel_abundance", ~group, "groupB - groupA")

}
\seealso{
\code{\link[plot_relative_abundance]{plot_relative_abundance}} for visualizing
the relative abundance data
}
