% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conduit.R
\docType{class}
\name{conduit-class}
\alias{conduit-class}
\title{Conduit Class for Proteomics Data Analysis}
\description{
A comprehensive class for storing and analyzing proteomics data, combining QFeatures
objects with taxonomic and detection metrics. This class serves as the central data
structure for the conduit package, enabling integrated analysis of protein abundance,
taxonomy, and detection statistics.
}
\section{Slots}{

\describe{
\item{\code{QFeatures}}{A QFeatures object containing proteomics data, including:
\itemize{
\item Protein abundance measurements across samples
\item Sample metadata (colData)
\item Feature metadata (rowData)
\item Multiple assays (e.g., raw, normalized, log-transformed data)
}}

\item{\code{combined_metrics}}{A tibble containing comprehensive metrics about protein detection,
including:
\itemize{
\item Number of proteins detected at each taxonomic level
\item Detection rates and coverage statistics
\item Quality metrics for protein identification
}}

\item{\code{database_protein_taxonomy}}{A tibble containing taxonomy information for all proteins
in the database, including:
\itemize{
\item Complete taxonomic classification (domain to species)
\item Organism type and source information
\item Protein identifiers and annotations
}}

\item{\code{detected_protein_taxonomy}}{A tibble containing taxonomy information specifically for
detected proteins, including:
\itemize{
\item Taxonomic classification of detected proteins
\item Detection statistics for each taxon
\item Quality metrics for identified proteins
}}
}}

\note{
The conduit class is designed to work seamlessly with other functions in the package:
\itemize{
\item Use \code{\link[create_conduit_obj]{create_conduit_obj}} to create new instances
\item Use \code{\link[calc_percent_proteins_detected]{calc_percent_proteins_detected}}
for detection rate analysis
\item Use \code{\link[plot_taxa_tree]{plot_taxa_tree}} and
\code{\link[plot_sunburst]{plot_sunburst}} for taxonomic visualization
\item Use \code{\link[perform_limma_analysis]{perform_limma_analysis}} for
differential expression analysis
}
For large datasets, consider:
\itemize{
\item Using appropriate filtering before creating the object
\item Managing memory usage when working with multiple assays
\item Utilizing the package's visualization functions for efficient data exploration
}
}
\examples{
# Create a conduit object from workflow outputs:
# conduit_obj <- create_conduit_obj(
#   QFeatures.rds = "results/qfeatures.rds",
#   combined_metrics.tsv = "results/metrics.tsv",
#   database_protein_taxonomy.tsv = "results/database_taxonomy.tsv",
#   detected_protein_taxonomy.tsv = "results/detected_taxonomy.tsv"
# )

# Access and analyze the data:
# - View summary of the object
# show(conduit_obj)

# - Access QFeatures data
# assays <- assays(conduit_obj@QFeatures)
# metadata <- colData(conduit_obj@QFeatures)

# - Analyze protein detection rates
# rates <- calc_percent_proteins_detected(conduit_obj)

# - Create taxonomic visualizations
# plot_taxa_tree(conduit_obj@database_protein_taxonomy)
# plot_sunburst(conduit_obj@detected_protein_taxonomy)

}
\seealso{
\itemize{
\item \code{\link[create_conduit_obj]{create_conduit_obj}} for creating new instances
\item \code{\link[QFeatures]{QFeatures}} for details about the QFeatures class
\item \code{\link[plot_taxa_tree]{plot_taxa_tree}} and
\code{\link[plot_sunburst]{plot_sunburst}} for visualization
\item \code{\link[calc_percent_proteins_detected]{calc_percent_proteins_detected}}
for detection rate analysis
}
}
