#' Plot Confusion Matrix for Classification Results
#'
#' Creates a heatmap visualization of the confusion matrix from a classification model,
#' showing the distribution of predicted versus true class labels.
#'
#' @param predict_classification_list A list containing model results, including:
#'   \itemize{
#'     \item confusion_matrix: A confusion matrix object (typically from yardstick::conf_mat)
#'       containing the cross-tabulation of predicted and true class labels
#'   }
#'
#' @return A ggplot object containing a heatmap with:
#'   \itemize{
#'     \item X-axis showing predicted class labels
#'     \item Y-axis showing true class labels
#'     \item Color intensity representing the count/frequency of each prediction
#'     \item Title and axis labels
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- predict_classification(qf, "protein_groups", "group")
#' plot_confusion_matrix(res)
#' }
plot_confusion_matrix <- function(predict_classification_list) {

  confusion_matrix <- predict_classification_list$confusion_matrix

  p1 <- ggplot2::autoplot(confusion_matrix,type = "heatmap") +
    ggplot2::labs(
      title = "Confusion Matrix",
      x = "Predicted Class",
      y = "True Class",
      fill = "Count"
    )

  return(p1)
}
