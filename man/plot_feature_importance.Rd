% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_feature_importance.R
\name{plot_feature_importance}
\alias{plot_feature_importance}
\title{Plot Feature Importance from Classification or Regression Model}
\usage{
plot_feature_importance(predict_classification_list, start = 1, end = 10)
}
\arguments{
\item{predict_classification_list}{A list containing model results, including:
\itemize{
\item importance: A data frame with columns 'feature' and 'importance' containing
feature importance scores
}}

\item{start}{Integer specifying the starting rank of features to plot (default: 1)}

\item{end}{Integer specifying the ending rank of features to plot (default: 10)}
}
\value{
A ggplot object containing a horizontal bar plot with:
\itemize{
\item Y-axis showing feature names, ordered by importance
\item X-axis showing importance scores
\item Horizontal bars representing feature importance
}
}
\description{
Creates a horizontal bar plot showing the importance of features (e.g., proteins)
from a trained machine learning model. Features are ordered by their importance score.
}
\examples{
# Plot top 10 most important features:
# plot_feature_importance(model_results)

# Plot features ranked 5-15:
# plot_feature_importance(model_results, start = 5, end = 15)

# Plot only top 5 features:
# plot_feature_importance(model_results, end = 5)
}
