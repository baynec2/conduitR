#' plot_roc
#'
#' @param predict_classification_list
#' @param type
#'
#' @returns
#' @export
#'
#' @examples
plot_roc = function(predict_classification_list,
                    data_set = "test"){

  if(data_set == "test"){
  # Generate ROC data
  roc_data <- yardstick::roc_curve(
    data = predict_classification_list$test_predictions,
    truth = predict_classification_list$outcome,
    names(predict_classification_list$test_predictions)[2]  # assumes binary classification with levels 0 and 1, and `.pred_1` is the prob for class 1
  )

  roc_auc <- yardstick::roc_auc(
    data = predict_classification_list$test_predictions,
    truth = predict_classification_list$outcome,
    names(predict_classification_list$test_predictions)[2]  # assumes binary classification with levels 0 and 1, and `.pred_1` is the prob for class 1
  )$.estimate

  } else if(data_set == "training"){
    # Generate ROC data
    roc_data <- yardstick::roc_curve(
      data = predict_classification_list$training_predictions,
      truth = predict_classification_list$outcome,
      names(predict_classification_list$training_predictions)[2]  # assumes binary classification with levels 0 and 1, and `.pred_1` is the prob for class 1
    )

    roc_auc <- yardstick::roc_auc(
      data = predict_classification_list$training_predictions,
      truth = predict_classification_list$outcome,
      names(predict_classification_list$training_predictions)[2]  # assumes binary classification with levels 0 and 1, and `.pred_1` is the prob for class 1
    )$.estimate


  } else{
    stop("data_set must be test or training")

  }

  # Plot via autoplot
  p1 = ggplot2::autoplot(roc_data)+
    ggplot2::labs(
    title = paste("ROC Curve -", data_set, "Data"),
    subtitle = paste("AUC = ", round(roc_auc, 2)
                     )
    )

  return(p1)
}
