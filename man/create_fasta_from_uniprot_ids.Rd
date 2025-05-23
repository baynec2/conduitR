% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_fasta_from_uniprot_ids.R
\name{create_fasta_from_uniprot_ids}
\alias{create_fasta_from_uniprot_ids}
\title{create_fasta_from_uniprot_ids
Creates a .fasta database from a list of uniprot ids.
file from the first proteome returned via uniprots API.}
\usage{
create_fasta_from_uniprot_ids(
  uniprot_ids,
  fasta_out_fp = paste0(getwd(), "/", Sys.Date(), ".fasta")
)
}
\arguments{
\item{uniprot_ids}{= vector of uniprot ids}

\item{fasta_out_fp}{= file path to write the .fasta database too}
}
\value{
a downloaded concatinated
}
\description{
create_fasta_from_uniprot_ids
Creates a .fasta database from a list of uniprot ids.
file from the first proteome returned via uniprots API.
}
\examples{
# Human and mouse insulin, will pull down the human and mouse proteome.
create_fasta_from_uniprot_ids(c("P01325", "P01308"), fasta_out_fp = "test.fasta")
}
