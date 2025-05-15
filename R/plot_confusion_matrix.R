#' plot a confusion matrix
#'
#' @param predict_classification_list
#'
#' @returns
#' @export
#'
#' @examples
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
