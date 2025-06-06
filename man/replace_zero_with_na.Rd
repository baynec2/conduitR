% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/replace_zero_with_na.R
\name{replace_zero_with_na}
\alias{replace_zero_with_na}
\title{Replace Zero Values with NA in QFeatures Object}
\usage{
replace_zero_with_na(qf)
}
\arguments{
\item{qf}{A QFeatures object containing the data to process}
}
\value{
A QFeatures object with all zero values replaced by NA in all assays.
The structure and metadata of the object remain unchanged.
}
\description{
Replaces all zero values with NA in all assays of a QFeatures object.
This is typically done as a preprocessing step before log transformation
and imputation, as zeros in proteomics data often represent missing values
rather than true zero abundances.
}
\examples{
# Replace zeros with NA in all assays:
# qf_no_zeros <- replace_zero_with_na(qfeatures_obj)

# This is often used as a preprocessing step:
# qf_processed <- qfeatures_obj |>
#   replace_zero_with_na() |>
#   add_log_imputed_norm_assays()

# The resulting object can be used with plotting functions:
# plot_missing_val_heatmap(qf_no_zeros, "protein")
}
