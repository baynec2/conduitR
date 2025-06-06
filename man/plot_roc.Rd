% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_roc.R
\name{plot_roc}
\alias{plot_roc}
\title{Plot ROC Curve for Classification Results}
\usage{
plot_roc(predict_classification_list, data_set = "test")
}
\arguments{
\item{predict_classification_list}{A list containing prediction results, including:
\itemize{
\item test_predictions: A data frame with prediction probabilities and true outcomes for test data
\item training_predictions: A data frame with prediction probabilities and true outcomes for training data
\item outcome: The name of the outcome variable
}}

\item{data_set}{Character string specifying which dataset to plot. Must be either "test" or "training".
Defaults to "test".}
}
\value{
A ggplot object containing the ROC curve plot with AUC score in the subtitle.
}
\description{
Creates a Receiver Operating Characteristic (ROC) curve plot for classification model predictions,
showing the trade-off between true positive rate and false positive rate.
}
\examples{
# Assuming you have a classification model output:
# plot_roc(model_predictions, data_set = "test")
# plot_roc(model_predictions, data_set = "training")
}
