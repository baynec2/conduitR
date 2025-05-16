#' Predict Numeric Outcomes Using Machine Learning Models
#'
#' Performs regression analysis on QFeatures data using various machine learning models.
#' Supports lasso regression, random forest, and XGBoost models with optional hyperparameter tuning.
#' This function handles data preprocessing, model training, and prediction in a single workflow.
#'
#' @param qf A QFeatures object containing the data to analyze
#' @param assay_name Character string specifying which assay to use for prediction
#' @param outcome Character string specifying the column name of the numeric outcome variable
#' @param train_percent Numeric value between 0 and 100 specifying the percentage of data to use for training
#'   (default: 70)
#' @param model_type Character string specifying the type of model to use. Must be one of:
#'   \itemize{
#'     \item "lasso_regression": L1-regularized linear regression
#'     \item "random_forest": Random forest regression
#'     \item "xgboost": Gradient boosting regression
#'   }
#' @param v Integer specifying the number of folds for cross-validation when tuning (default: 5)
#' @param grid_size Integer specifying the number of parameter combinations to try during tuning (default: 20)
#'
#' @return A list containing:
#'   \itemize{
#'     \item fit: The trained model object
#'     \item training_predictions: Predictions on the training data
#'     \item test_predictions: Predictions on the test data
#'     \item outcome: Name of the outcome variable
#'     \item test: The test dataset
#'     \item importance: Feature importance scores (varies by model type)
#'   }
#'
#' @export
#'
#' @examples
#' # Basic usage with lasso regression:
#' # results <- predict_regression(qfeatures_obj, "protein", "concentration")
#' 
#' # Using random forest with tuning:
#' # results <- predict_regression(qfeatures_obj, "protein", "concentration",
#' #                             model_type = "random_forest", v = 5)
#' 
#' # Using XGBoost with custom training split:
#' # results <- predict_regression(qfeatures_obj, "protein", "concentration",
#' #                             model_type = "xgboost", train_percent = 80)
#'
#' @note
#' This function requires the following packages:
#' \itemize{
#'   \item glmnet for lasso regression
#'   \item ranger for random forest
#'   \item xgboost for gradient boosting
#'   \item caret for model training and tuning
#' }
#' Computation time varies by model type:
#' \itemize{
#'   \item Lasso regression: Fastest, suitable for large datasets
#'   \item Random forest: Moderate, scales with number of trees and features
#'   \item XGBoost: Can be slow for large datasets, but often provides best performance
#' }
#' For large datasets, consider using a subset of features or increasing the removeVar
#' parameter in preprocessing steps to improve performance.
predict_regression <- function(qf,
                               assay_name,
                               outcome,
                               train_percent = 70,
                               model_type = "lasso_regression",
                               v = 5,
                               grid_size = 20) {

  supported_models <- c("lasso_regression", "random_forest", "xgboost")

  if (!(model_type %in% supported_models)) {
    stop("model_type must be one of ", paste(supported_models, collapse = ", "))
  } else {
    message(model_type, " model_type is supported")
  }

  if (v > 1) {
    tune <- TRUE
  } else if (v == 1) {
    tune <- FALSE
  } else {
    stop("v must be >= 1")
  }

  if (!inherits(qf, "QFeatures")) stop("qf must be a QFeatures object.")
  if (!assay_name %in% names(qf)) stop("Invalid assay name.")

  tidy <- conduitR::tidy_conduit(qf, assay_name)

  selected_data <- tidy |>
    dplyr::select(sample, !!rlang::sym(outcome), rowid, value) |>
    tidyr::pivot_wider(names_from = rowid, values_from = value)

  split <- rsample::initial_split(selected_data, prop = train_percent / 100,
                                  strata = outcome)
  train <- rsample::training(split)
  test <- rsample::testing(split)

  if (!is.numeric(train[[outcome]])) {
    stop("This function only supports numeric outcomes for regression.")
  }

  rec <- recipes::recipe(as.formula(paste(outcome, "~ .")), data = train) |>
    recipes::update_role(sample, new_role = "id") |>
    recipes::update_role(outcome, new_role = "outcome")

  if (model_type == "lasso_regression") {
    rec <- rec |>
      recipes::step_impute_knn(recipes::all_numeric_predictors(), skip = TRUE)
  }

  rec <- rec |>
    recipes::step_zv(recipes::all_numeric_predictors())

  if (model_type %in% c("lasso_regression", "xgboost")) {
    rec <- rec |>
      recipes::step_normalize(recipes::all_numeric_predictors())
  }

  # Model spec with or without tuning
  if (tune) {
    model_spec <- switch(model_type,
                         "lasso_regression" = parsnip::linear_reg(penalty = tune(), mixture = 1) |>
                           parsnip::set_engine("glmnet") |>
                           parsnip::set_mode("regression"),

                         "random_forest" = parsnip::rand_forest(mtry = tune(), min_n = tune(), trees = tune()) |>
                           parsnip::set_engine("ranger", importance = "impurity") |>
                           parsnip::set_mode("regression"),

                         "xgboost" = parsnip::boost_tree(
                             trees = 1000,
                             tree_depth = tune(),
                             learn_rate = tune(),
                             mtry = tune(),
                             loss_reduction = tune(),
                             sample_size = tune()
                           ) |>
                             parsnip::set_engine("xgboost") |>
                             parsnip::set_mode("regression")
    )
  } else {
    model_spec <- switch(model_type,
                         "lasso_regression" = parsnip::linear_reg(penalty = 0.1, mixture = 1) |>
                           parsnip::set_engine("glmnet") |>
                           parsnip::set_mode("regression"),

                         "random_forest" = parsnip::rand_forest(trees = 500) |>
                           parsnip::set_engine("ranger", importance = "impurity") |>
                           parsnip::set_mode("regression"),

                         "xgboost" = parsnip::boost_tree(trees = 100, mtry = 3) |>
                           parsnip::set_engine("xgboost") |>
                           parsnip::set_mode("regression")
    )
  }

  wf <- workflows::workflow() |>
    workflows::add_recipe(rec) |>
    workflows::add_model(model_spec)

  if (tune) {
    cv_folds <- rsample::vfold_cv(train, v = v)

    # Define model-specific parameters
    if (model_type == "lasso_regression") {
      # Lasso Regression
      tuning_params <- dials::parameters(
        penalty(range = c(0.001, 1))
        )
      tuning_grid <- dials::grid_regular(tuning_params, levels = 10)  # Generate grid with 10 levels for each parameter
    } else if (model_type == "random_forest") {
      # Random Forest
      tuning_params <- dials::parameters(
        mtry(range = c(1, 10)),      # Tune mtry (number of variables to try at each split)
        min_n(range = c(1, 20)),     # Tune min_n (minimum number of data points required for leaf)
        trees(range = c(50, 2000))   # Tune number of trees
      )
      tuning_grid <- dials::grid_regular(tuning_params, levels = 5)  # Generate grid with 5 levels for each parameter
    } else if (model_type == "xgboost") {
      # XGBoost
      tuning_params <- dials::parameters(
        tree_depth(range = c(3, 10)),
        learn_rate(range = c(0.01, 0.3)),
        mtry(range = c(2, 10)),
        sample_prop(range = c(0.5, 1.0)),
        loss_reduction(range = c(0, 5))
      )
      tuning_grid <- dials::grid_regular(tuning_params, levels = 5)  # Generate grid with 5 levels for each parameter
    } else {
      stop("Unknown model type. Supported types are: lasso_regression, random_forest, xgboost")
    }

    # Finalize the grid using training data (train)
    finalized_params <- dials::finalize(tuning_params,train)

    # Create the grid using space filling
    grid <- dials::grid_space_filling(finalized_params, size = grid_size)

    # Perform the grid search using tune_grid
    tuned <- tune::tune_grid(
      wf,
      resamples = cv_folds,
      grid = grid,
      metrics = yardstick::metric_set(yardstick::rmse)
    )

    # Get the best parameters based on RMSE
    best_params <- tune::select_best(tuned, metric = "rmse")
    wf <- tune::finalize_workflow(wf, best_params)
  }


  fit <- parsnip::fit(wf, data = train)

  # Predictions
  preds <- dplyr::bind_cols(
    predict(fit, test),
    test[, outcome, drop = FALSE]
  )

  training_predictions <- dplyr::bind_cols(
    predict(fit, train),
    train[, outcome, drop = FALSE]
  )

  # Variable importance
  if (model_type == "random_forest") {
    rf_fit <- workflows::extract_fit_parsnip(fit)
    importance <- tibble::enframe(rf_fit$fit$variable.importance,
                                  name = "feature", value = "importance") |>
      dplyr::arrange(desc(importance))
  } else if (model_type == "xgboost") {
    xgb_fit <- workflows::extract_fit_parsnip(fit)
    importance <- tibble::as_tibble(xgboost::xgb.importance(model = xgb_fit$fit)) |>
      dplyr::rename(feature = Feature, importance = Gain) |>
      dplyr::arrange(desc(importance))
  } else if (model_type == "lasso_regression") {
    lasso_fit <- workflows::extract_fit_parsnip(fit)
    coefs <- as.matrix(coef(lasso_fit$fit, s = lasso_fit$spec$args$penalty))[, 1]
    importance <- tibble::tibble(feature = names(coefs), importance = coefs) |>
      dplyr::filter(feature != "(Intercept)") |>
      dplyr::arrange(desc(abs(importance)))
  }

  output <- list(
    fit = fit,
    training_predictions = training_predictions,
    test_predictions = preds,
    outcome = outcome,
    test = test,
    importance = importance
  )

  return(output)
}
