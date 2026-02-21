#' Classification Model Prediction with Optional Tuning
#'
#' Fits a classification model using lasso regression, random forest, or XGBoost
#' on a `QFeatures` object to a specified outcome variable in colData.
#' Optionally performs hyperparameter tuning using cross-validation.
#'
#' @param qf A `QFeatures` object.
#' @param assay_name Name of the assay within `qf` to use.
#' @param outcome Column name of the classification outcome variable.
#' @param train_percent Percentage of data to use for training. Default is 70.
#' @param model_type One of `"lasso_regression"`, `"random_forest"`, or `"xgboost"`.
#' @param v Number of folds for cross-validation. Default is 5. If set to 1, tuning is skipped.
#' @param grid_size Number of tuning grid points. Default is 20.
#'
#' @return A list with fitted model, test and training predictions, confusion matrix, and variable importance.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- predict_classification(qf, "protein_groups", "group", model_type = "random_forest")
#' plot_confusion_matrix(res)
#' plot_feature_importance(res)
#' }
predict_classification <- function(qf,
                                   assay_name,
                                   outcome,
                                   train_percent = 70,
                                   model_type = "random_forest",
                                   v = 5,
                                   grid_size = 20) {
  # Supported models
  supported_models <- c("lasso_regression", "random_forest", "xgboost")
  if (!model_type %in% supported_models) {
    stop("model_type must be one of: ", paste(supported_models, collapse = ", "))
  } else {
    message(model_type, " model_type is supported.")
  }

  # Check cross-validation folds
  if (v < 1) stop("v must be >= 1")
  tune <- v > 1

  # Check QFeatures input and assay
  if (!inherits(qf, "QFeatures")) stop("qf must be a QFeatures object.")
  if (!assay_name %in% names(qf)) stop("Invalid assay name: ", assay_name)

  # Tidy data
  tidy <- conduitR::tidy_conduit(qf, assay_name)

  selected_data <- tidy |>
    dplyr::select(file, !!rlang::sym(outcome), rowid, value) |>
    tidyr::pivot_wider(names_from = rowid, values_from = value) |>
    dplyr::mutate(!!rlang::sym(outcome) := as.factor(!!rlang::sym(outcome)))

  # Stopping if there are NAs an a lasso regression model.
  if (model_type == "lasso_regression" && anyNA(tidy$value)) {
    stop(
      "Lasso regression does not support missing values (NA). ",
      "Please impute missing data or choose a model like XGBoost that can handle
      them."
    )
  }

  # Split train/test
  split <- rsample::initial_split(
    selected_data,
    prop = train_percent / 100,
    strata = !!rlang::sym(outcome)
  )
  train <- rsample::training(split)
  test <- rsample::testing(split)

  if (is.numeric(train[[outcome]])) stop("This function only supports classification.")

  # Define recipe
  rec <- recipes::recipe(as.formula(paste(outcome, "~ .")), data = train) |>
    recipes::update_role(file, new_role = "id") |>
    recipes::step_rm(file) |>
    recipes::step_zv()

  if (model_type %in% c("lasso_regression", "xgboost")) {
    rec <- rec |> recipes::step_normalize(recipes::all_numeric_predictors())
  }

  # Model specification
  if (tune) {
    model_spec <- switch(model_type,
                         "lasso_regression" = parsnip::logistic_reg(penalty = tune(), mixture = 1) |>
                           parsnip::set_engine("glmnet") |>
                           parsnip::set_mode("classification"),
                         "random_forest" = parsnip::rand_forest(mtry = tune(), min_n = tune(), trees = 500) |>
                           parsnip::set_engine("ranger", importance = "impurity") |>
                           parsnip::set_mode("classification"),
                         "xgboost" = parsnip::boost_tree(
                           trees = 1000,
                           tree_depth = tune(),
                           learn_rate = tune(),
                           mtry = tune(),
                           loss_reduction = tune(),
                           sample_size = tune()
                         ) |>
                           parsnip::set_engine("xgboost") |>
                           parsnip::set_mode("classification")
    )
  } else {
    model_spec <- switch(model_type,
                         "lasso_regression" = parsnip::logistic_reg(penalty = 0.1, mixture = 1) |>
                           parsnip::set_engine("glmnet") |>
                           parsnip::set_mode("classification"),
                         "random_forest" = parsnip::rand_forest(trees = 500) |>
                           parsnip::set_engine("ranger", importance = "impurity") |>
                           parsnip::set_mode("classification"),
                         "xgboost" = parsnip::boost_tree(trees = 100, mtry = 3) |>
                           parsnip::set_engine("xgboost") |>
                           parsnip::set_mode("classification")
    )
  }

  # Workflow
  wf <- workflows::workflow() |>
    workflows::add_recipe(rec) |>
    workflows::add_model(model_spec)

  # Choose metric set
  if (nlevels(train[[outcome]]) == 2) {
    metrics <- yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy)
  } else {
    metrics <- yardstick::metric_set(yardstick::accuracy)
  }

  # Tuning if enabled
  if (tune) {
    cv_folds <- rsample::vfold_cv(train, v = v, strata = !!rlang::sym(outcome))

    finalized_params <- hardhat::extract_parameter_set_dials(model_spec) |>
      dials::finalize(train)

    grid <- dials::grid_space_filling(finalized_params, size = grid_size)

    tuned <- tune::tune_grid(
      wf,
      resamples = cv_folds,
      grid = grid,
      metrics = metrics
    )

    best_params <- tune::select_best(tuned, metric = ifelse("roc_auc" %in% names(tuned$.metrics[[1]]), "roc_auc", "accuracy"))
    wf <- tune::finalize_workflow(wf, best_params)
  }

  # Fit workflow
  fit <- workflows::fit(wf, data = train)

  # Predictions
  preds <- dplyr::bind_cols(
    predict(fit, test, type = "prob"),
    predict(fit, test),
    test[, outcome, drop = FALSE]
  )

  training_predictions <- dplyr::bind_cols(
    predict(fit, train, type = "prob"),
    predict(fit, train),
    train[, outcome, drop = FALSE]
  )

  # Feature importance
  importance <- NULL
  if (model_type == "random_forest") {
    rf_fit <- workflows::extract_fit_parsnip(fit)
    if (!is.null(rf_fit$fit$variable.importance)) {
      importance <- tibble::enframe(rf_fit$fit$variable.importance,
                                    name = "feature", value = "importance"
      ) |> dplyr::arrange(desc(importance))
    }
  } else if (model_type == "xgboost") {
    xgb_fit <- workflows::extract_fit_parsnip(fit)
    importance <- tryCatch(
      tibble::as_tibble(xgboost::xgb.importance(model = xgb_fit$fit)) |>
        dplyr::rename(feature = Feature, importance = Gain) |>
        dplyr::arrange(desc(importance)),
      error = function(e) NULL
    )
  } else if (model_type == "lasso_regression") {
    lasso_fit <- workflows::extract_fit_parsnip(fit)
    lambda <- tryCatch(lasso_fit$fit$lambdaOpt, error = function(e) lasso_fit$spec$args$penalty)
    coefs <- as.matrix(coef(lasso_fit$fit, s = lambda))[, 1]
    importance <- tibble::tibble(feature = names(coefs), importance = coefs) |>
      dplyr::filter(feature != "(Intercept)") |>
      dplyr::arrange(desc(abs(importance)))
  }

  # Confusion matrix
  confusion_matrix <- yardstick::conf_mat(
    data = preds,
    truth = !!rlang::sym(outcome),
    estimate = !!rlang::sym(".pred_class")
  )

  # Metrics
  performance <- yardstick::metrics(
    data = preds,
    truth = !!rlang::sym(outcome),
    estimate = !!rlang::sym(".pred_class")
  )

  # Return output
  output <- list(
    fit = fit,
    confusion_matrix = confusion_matrix,
    metrics = performance,
    training_predictions = training_predictions,
    test_predictions = preds,
    outcome = outcome,
    test = test,
    importance = importance
  )

  return(output)
}
