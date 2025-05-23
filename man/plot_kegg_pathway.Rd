% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_kegg_pathway.R
\name{plot_kegg_pathway}
\alias{plot_kegg_pathway}
\title{plot_kegg_pathway}
\usage{
plot_kegg_pathway(
  qf,
  assay_name = "log2_imputed",
  kegg_pathway_id = "03018",
  ids_to_color = c(K04077 = 1)
)
}
\arguments{
\item{qf}{q features object}

\item{assay_name}{assay name with IDs}

\item{kegg_pathway_id}{kegg pathway ID}

\item{ids_to_color}{ids to plot.}
}
\description{
plot values in data on kegg pathway.
}
