#' Plot Precision-Recall Curve for Classification Results
#'
#' Creates a precision-recall curve plot for classification model predictions,
#' showing the trade-off between precision and recall at different classification thresholds.
#' This is particularly useful for imbalanced classification problems.
#'
#' @param predict_classification_list A list containing prediction results, including:
#'   \itemize{
#'     \item test_predictions: A data frame with prediction probabilities and true outcomes for test data
#'     \item training_predictions: A data frame with prediction probabilities and true outcomes for training data
#'     \item outcome: The name of the outcome variable
#'   }
#' @param data_set Character string specifying which dataset to plot. Must be either "test" or "training".
#'   Defaults to "test".
#'
#' @return A ggplot object containing the precision-recall curve plot with:
#'   \itemize{
#'     \item X-axis showing recall
#'     \item Y-axis showing precision
#'     \item Title indicating the dataset used
#'     \item Subtitle showing the area under the precision-recall curve (AUC-PR)
#'   }
#'
#' @export
#'
#' @examples
#' # Plot precision-recall curve for test data:
#' # plot_precision_recall(model_predictions, data_set = "test")
#' 
#' # Plot precision-recall curve for training data:
#' # plot_precision_recall(model_predictions, data_set = "training")
plot_precision_recall <- function(predict_classification_list,data_set = "test") {

  if(data_set == "test") {
    # Generate Precision-Recall data for test set
    pr_data <- yardstick::pr_curve(
      data = predict_classification_list$test_predictions,
      truth = predict_classification_list$outcome,
      names(predict_classification_list$test_predictions)[2]  # assumes binary classification, and `.pred_1` is the prob for class 1
    )
    pr_auc <- yardstick::pr_auc( data = predict_classification_list$test_predictions,
                                 truth = predict_classification_list$outcome,
                                 names(predict_classification_list$test_predictions)[2] )$.estimate
  } else if(data_set == "training") {
    # Generate Precision-Recall data for training set
    pr_data <- yardstick::pr_curve(
      data = predict_classification_list$training_predictions,
      truth =  predict_classification_list$outcome,
      names(predict_classification_list$training_predictions)[2]  # assumes binary classification, and `.pred_1` is the prob for class 1
    )
    pr_auc <- yardstick::pr_auc(
      data = predict_classification_list$training_predictions,
      truth =  predict_classification_list$outcome,
      names(predict_classification_list$training_predictions)[2]  # assumes binary classification, and `.pred_1` is the prob for class 1
    )$.estimate
  } else {
    stop("type must be 'test' or 'training'")
  }

  # Plot via autoplot
  p1 = ggplot2::autoplot(pr_data)+
    ggplot2::labs(
      title = paste("Precision Recall Curve -", data_set, "Data"),
      subtitle = paste("AUC = ", round(pr_auc, 2)
      )
    )

  return(p1)
}
