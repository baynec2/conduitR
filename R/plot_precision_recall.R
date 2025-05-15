#' plot precision recall curves from predicted classifications
#'
#' @param predict_classification_list
#' @param data_set
#'
#' @returns
#' @export
#'
#' @examples
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
