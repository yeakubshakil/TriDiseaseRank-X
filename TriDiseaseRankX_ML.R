############################################################
# TriDiseaseRank-X Stage 1: generic internal 113-model ML engine
#
# Purpose:
#   Run the user's 113 machine-learning pipelines on one disease cohort
#   using only internal nested cross-validation. No external validation,
#   ANN, pathway, or regulatory analysis is performed here.
#
# Candidate genes are supplied dynamically from TD-CCHI/GNN output.
# Model ranking uses internal out-of-fold AUC only.
# Gene-level ML evidence is derived from feature-selection stability among
# the top internally performing selector-based pipelines.
############################################################

required_ml_packages <- c(
  "caret", "pROC", "glmnet", "ranger", "mboost", "MASS",
  "plsRglm", "gbm", "xgboost", "naivebayes", "kernlab",
  "data.table", "matrixStats", "dplyr", "tidyr", "future", "future.apply"
)
missing_ml_packages <- required_ml_packages[!vapply(required_ml_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_ml_packages)) {
  stop("Missing ML packages: ", paste(missing_ml_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(caret); library(pROC); library(glmnet); library(ranger); library(mboost)
  library(MASS); library(plsRglm); library(gbm); library(xgboost); library(naivebayes)
  library(kernlab); library(data.table); library(matrixStats); library(dplyr); library(tidyr)
  library(future); library(future.apply)
})

options(stringsAsFactors = FALSE, scipen = 999)

# Global defaults used by the inherited model helper functions.
case_label <- "Case"
control_label <- "Control"
inner_folds_requested <- 5L
tuning_length <- 5L

model_definitions <- c(
  "Lasso + Stepglm [both]",
  "SVM",
  "glmBoost + SVM",
  "Ridge",
  "Lasso + SVM",
  "glmBoost + Ridge",
  "Enet [alpha=0.1]",
  "glmBoost + Enet [alpha=0.1]",
  "Enet [alpha=0.2]",
  "Enet [alpha=0.3]",
  "glmBoost + Enet [alpha=0.3]",
  "glmBoost + Enet [alpha=0.2]",
  "Enet [alpha=0.4]",
  "glmBoost + Enet [alpha=0.4]",
  "Lasso + glmBoost",
  "Enet [alpha=0.5]",
  "glmBoost",
  "glmBoost + Enet [alpha=0.5]",
  "Enet [alpha=0.6]",
  "glmBoost + Enet [alpha=0.6]",
  "glmBoost + Enet [alpha=0.7]",
  "glmBoost + Enet [alpha=0.8]",
  "Enet [alpha=0.8]",
  "Enet [alpha=0.9]",
  "Lasso",
  "Enet [alpha=0.7]",
  "glmBoost + Enet [alpha=0.9]",
  "glmBoost + Lasso",
  "Lasso + plsRglm",
  "glmBoost + plsRglm",
  "glmBoost + Stepglm [forward]",
  "Lasso + Stepglm [forward]",
  "RF + SVM",
  "Stepglm [forward]",
  "plsRglm",
  "RF + Ridge",
  "RF + Enet [alpha=0.1]",
  "RF + plsRglm",
  "RF + Stepglm [forward]",
  "RF + Enet [alpha=0.2]",
  "RF + Enet [alpha=0.3]",
  "RF + Enet [alpha=0.6]",
  "RF + Lasso",
  "RF + Enet [alpha=0.7]",
  "RF + Enet [alpha=0.5]",
  "RF + glmBoost",
  "RF + Enet [alpha=0.9]",
  "RF + Enet [alpha=0.4]",
  "RF + Enet [alpha=0.8]",
  "RF + Stepglm [both]",
  "RF + Stepglm [backward]",
  "Stepglm [both] + Ridge",
  "Stepglm [backward] + Ridge",
  "Stepglm [both] + plsRglm",
  "Stepglm [backward] + plsRglm",
  "Stepglm [both] + Enet [alpha=0.9]",
  "Stepglm [backward] + Enet [alpha=0.9]",
  "Stepglm [both] + Enet [alpha=0.1]",
  "Stepglm [backward] + Enet [alpha=0.1]",
  "Stepglm [both] + Enet [alpha=0.8]",
  "Stepglm [backward] + Enet [alpha=0.8]",
  "Stepglm [both] + Enet [alpha=0.2]",
  "Stepglm [backward] + Enet [alpha=0.2]",
  "Stepglm [both] + Lasso",
  "Stepglm [backward] + Lasso",
  "Stepglm [both] + Enet [alpha=0.6]",
  "Stepglm [backward] + Enet [alpha=0.6]",
  "glmBoost + GBM",
  "Stepglm [both] + Enet [alpha=0.7]",
  "Stepglm [backward] + Enet [alpha=0.7]",
  "Lasso + Stepglm [backward]",
  "Stepglm [both]",
  "Stepglm [backward]",
  "glmBoost + Stepglm [both]",
  "glmBoost + Stepglm [backward]",
  "Stepglm [both] + Enet [alpha=0.4]",
  "Stepglm [backward] + Enet [alpha=0.4]",
  "Stepglm [both] + Enet [alpha=0.3]",
  "Stepglm [backward] + Enet [alpha=0.3]",
  "Stepglm [both] + glmBoost",
  "Stepglm [backward] + glmBoost",
  "Stepglm [both] + Enet [alpha=0.5]",
  "Stepglm [backward] + Enet [alpha=0.5]",
  "glmBoost + RF",
  "RF",
  "Lasso + GBM",
  "RF + GBM",
  "GBM",
  "Stepglm [both] + SVM",
  "Stepglm [backward] + SVM",
  "Lasso + RF",
  "Stepglm [both] + GBM",
  "Stepglm [backward] + GBM",
  "Stepglm [both] + RF",
  "LDA",
  "glmBoost + LDA",
  "RF+LDA",
  "Stepglm [both]+LDA",
  "Stepglm [backward]+LDA",
  "Lasso + LDA",
  "Stepglm [backward] + RF",
  "XGBoost",
  "Lasso+XGBoost",
  "glmBoost+XGBoost",
  "RF+XGBoost",
  "Stepglm[both]+XGBoost",
  "Stepglm[backward]+XGBoost",
  "NaiveBayes",
  "Lasso+NaiveBayes",
  "glmBoost+NaiveBayes",
  "RF+NaiveBayes",
  "Stepglm[both]+NaiveBayes",
  "Stepglm[backward]+NaiveBayes"
)

stopifnot(length(model_definitions) == 113)


# 4. Helper functions: model parsing and metrics
# ============================================================

canonical_token <- function(token) {
  token <- gsub("[[:space:]]+", "", token)
  lower <- tolower(token)

  if (grepl("^enet\\[alpha=", lower)) {
    return(
      paste0(
        "Enet[alpha=",
        sub(
          "^.*alpha=([0-9.]+)\\].*$",
          "\\1",
          lower
        ),
        "]"
      )
    )
  }

  token_map <- c(
    "none" = "None",
    "lasso" = "Lasso",
    "ridge" = "Ridge",
    "svm" = "SVM",
    "glmboost" = "glmBoost",
    "plsrglm" = "plsRglm",
    "rf" = "RF",
    "gbm" = "GBM",
    "lda" = "LDA",
    "xgboost" = "XGBoost",
    "naivebayes" = "NaiveBayes",
    "stepglm[both]" = "Stepglm[both]",
    "stepglm[backward]" = "Stepglm[backward]",
    "stepglm[forward]" = "Stepglm[forward]"
  )

  if (lower %in% names(token_map)) {
    return(unname(token_map[lower]))
  }

  stop("Unsupported algorithm token: ", token)
}

parse_pipeline <- function(model_name) {
  compact <- gsub("[[:space:]]+", "", model_name)
  parts <- strsplit(compact, "\\+", fixed = FALSE)[[1]]

  if (length(parts) == 1) {
    list(
      selector = "None",
      classifier = canonical_token(parts[1])
    )
  } else if (length(parts) == 2) {
    list(
      selector = canonical_token(parts[1]),
      classifier = canonical_token(parts[2])
    )
  } else {
    stop("Unable to parse model definition: ", model_name)
  }
}

safe_divide <- function(numerator, denominator) {
  ifelse(denominator == 0, NA_real_, numerator / denominator)
}

compute_auc_object <- function(y, probability) {
  pROC::roc(
    response = y,
    predictor = probability,
    levels = c(control_label, case_label),
    direction = "<",
    quiet = TRUE
  )
}

compute_gene_auc_object <- function(y, expression_value) {
  pROC::roc(
    response = y,
    predictor = expression_value,
    levels = c(control_label, case_label),
    direction = "auto",
    quiet = TRUE
  )
}

derive_training_threshold <- function(y, probability) {
  roc_object <- compute_auc_object(y, probability)

  threshold <- pROC::coords(
    roc_object,
    x = "best",
    best.method = "youden",
    ret = "threshold",
    transpose = FALSE
  )

  as.numeric(threshold)[1]
}

compute_metrics <- function(
  y,
  probability,
  threshold,
  dataset_label
) {
  roc_object <- compute_auc_object(y, probability)
  auc_value <- as.numeric(pROC::auc(roc_object))
  auc_ci <- as.numeric(pROC::ci.auc(roc_object, method = "delong"))

  truth_case <- y == case_label
  predicted_case <- probability >= threshold

  tp <- sum(predicted_case & truth_case)
  tn <- sum(!predicted_case & !truth_case)
  fp <- sum(predicted_case & !truth_case)
  fn <- sum(!predicted_case & truth_case)

  sensitivity <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  precision <- safe_divide(tp, tp + fp)
  recall <- sensitivity
  accuracy <- safe_divide(tp + tn, length(y))
  balanced_accuracy <- mean(
    c(sensitivity, specificity),
    na.rm = TRUE
  )
  f1 <- if (
    is.na(precision) ||
      is.na(recall) ||
      precision + recall == 0
  ) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }

  brier <- mean(
    (as.numeric(truth_case) - probability)^2
  )

  data.frame(
    Dataset = dataset_label,
    AUC = auc_value,
    AUC_CI_Lower = auc_ci[1],
    AUC_CI_Upper = auc_ci[3],
    Threshold = threshold,
    Accuracy = accuracy,
    Balanced_Accuracy = balanced_accuracy,
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    Brier = brier,
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn,
    stringsAsFactors = FALSE
  )
}

univariate_fallback <- function(X, y, number = 2) {
  y_binary <- y == case_label

  effect_score <- vapply(
    X,
    function(values) {
      case_values <- values[y_binary]
      control_values <- values[!y_binary]

      pooled_sd <- sqrt(
        (
          stats::var(case_values) +
            stats::var(control_values)
        ) / 2
      )

      if (!is.finite(pooled_sd) || pooled_sd == 0) {
        return(0)
      }

      abs(
        mean(case_values) -
          mean(control_values)
      ) / pooled_sd
    },
    numeric(1)
  )

  names(
    sort(
      effect_score,
      decreasing = TRUE
    )
  )[seq_len(min(number, ncol(X)))]
}

# ============================================================
# 4B. Cross-cohort harmonization and balanced resampling
# ============================================================

standardize_expression_cohort <- function(X) {
  X <- as.data.frame(X, check.names = FALSE)

  means <- vapply(
    X,
    mean,
    numeric(1),
    na.rm = TRUE
  )

  standard_deviations <- vapply(
    X,
    stats::sd,
    numeric(1),
    na.rm = TRUE
  )

  standard_deviations[
    !is.finite(standard_deviations) |
      standard_deviations <= 1e-10
  ] <- 1

  scaled <- sweep(
    sweep(
      as.matrix(X),
      2,
      means,
      "-"
    ),
    2,
    standard_deviations,
    "/"
  )

  scaled <- as.data.frame(
    scaled,
    check.names = FALSE
  )

  rownames(scaled) <- rownames(X)
  scaled
}

dataset_from_rownames <- function(X) {
  sample_names <- rownames(X)

  if (
    is.null(sample_names) ||
      !all(grepl("__", sample_names, fixed = TRUE))
  ) {
    return(rep("SingleCohort", nrow(X)))
  }

  sub("__.*$", "", sample_names)
}

group_from_rownames <- function(X) {
  sample_names <- rownames(X)

  if (is.null(sample_names)) {
    return(paste0("Sample_", seq_len(nrow(X))))
  }

  parts <- strsplit(
    sample_names,
    "__",
    fixed = TRUE
  )

  vapply(
    seq_along(parts),
    function(index) {
      current <- parts[[index]]

      if (length(current) >= 3) {
        paste(
          current[1],
          current[2],
          sep = "__"
        )
      } else {
        sample_names[index]
      }
    },
    character(1)
  )
}

make_balanced_folds <- function(
  X,
  y,
  requested_folds,
  return_train = FALSE
) {
  dataset_id <- dataset_from_rownames(X)
  group_id <- group_from_rownames(X)

  unique_groups <- unique(group_id)

  group_dataset <- vapply(
    unique_groups,
    function(group_value) {
      unique(dataset_id[group_id == group_value])[1]
    },
    character(1)
  )

  group_outcome <- vapply(
    unique_groups,
    function(group_value) {
      observed <- unique(
        as.character(y[group_id == group_value])
      )

      if (length(observed) > 1) {
        "Paired"
      } else {
        observed[1]
      }
    },
    character(1)
  )

  group_strata <- interaction(
    group_dataset,
    group_outcome,
    drop = TRUE,
    lex.order = TRUE
  )

  minimum_stratum <- min(table(group_strata))
  number_folds <- min(
    requested_folds,
    minimum_stratum
  )

  if (number_folds < 3) {
    group_strata <- factor(group_outcome)
    number_folds <- min(
      requested_folds,
      min(table(group_strata))
    )
  }

  number_folds <- max(3, number_folds)

  group_assessment_folds <- caret::createFolds(
    group_strata,
    k = number_folds,
    returnTrain = FALSE
  )

  sample_assessment_folds <- lapply(
    group_assessment_folds,
    function(group_indices) {
      selected_groups <- unique_groups[group_indices]
      which(group_id %in% selected_groups)
    }
  )

  if (return_train) {
    all_indices <- seq_len(nrow(X))

    return(
      lapply(
        sample_assessment_folds,
        function(assessment_indices) {
          setdiff(all_indices, assessment_indices)
        }
      )
    )
  }

  sample_assessment_folds
}

# ============================================================
# 5. Feature-selection functions
# ============================================================

select_lasso_features <- function(X, y, seed) {
  set.seed(seed)

  inner_fold_list <- make_balanced_folds(
    X,
    y,
    requested_folds = inner_folds_requested,
    return_train = FALSE
  )

  fold_id <- rep(
    NA_integer_,
    nrow(X)
  )

  for (fold_index in seq_along(inner_fold_list)) {
    fold_id[inner_fold_list[[fold_index]]] <- fold_index
  }

  if (anyNA(fold_id)) {
    stop("Failed to assign all samples to inner LASSO folds.")
  }

  fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = as.numeric(y == case_label),
    family = "binomial",
    alpha = 1,
    foldid = fold_id,
    type.measure = "auc",
    standardize = TRUE,
    nlambda = 50
  )

  coefficients <- as.matrix(
    coef(fit, s = "lambda.1se")
  )

  selected <- rownames(coefficients)[
    abs(coefficients[, 1]) > 1e-10
  ]

  selected <- setdiff(selected, "(Intercept)")
  selected <- intersect(selected, colnames(X))

  if (!length(selected)) {
    coefficients <- as.matrix(
      coef(fit, s = "lambda.min")
    )

    selected <- rownames(coefficients)[
      abs(coefficients[, 1]) > 1e-10
    ]

    selected <- setdiff(selected, "(Intercept)")
    selected <- intersect(selected, colnames(X))
  }

  if (!length(selected)) {
    selected <- univariate_fallback(X, y, number = 2)
  }

  selected
}

select_rf_features <- function(X, y, seed) {
  set.seed(seed)

  model_data <- data.frame(
    Class = y,
    X,
    check.names = FALSE
  )

  fit <- ranger::ranger(
    Class ~ .,
    data = model_data,
    probability = TRUE,
    importance = "permutation",
    num.trees = 1000,
    mtry = max(1, floor(sqrt(ncol(X)))),
    seed = seed
  )

  importance <- sort(
    fit$variable.importance,
    decreasing = TRUE
  )

  positive <- importance[is.finite(importance) & importance > 0]

  if (length(positive)) {
    cutoff <- stats::median(positive)
    selected <- names(positive)[positive >= cutoff]
  } else {
    selected <- character(0)
  }

  maximum_selected <- min(10, ncol(X))

  if (length(selected) > maximum_selected) {
    selected <- names(importance)[seq_len(maximum_selected)]
  }

  if (length(selected) < 2) {
    selected <- names(importance)[
      seq_len(min(2, length(importance)))
    ]
  }

  if (!length(selected)) {
    selected <- univariate_fallback(X, y, number = 2)
  }

  selected
}

select_glmboost_features <- function(X, y, seed) {
  set.seed(seed)

  inner_index <- make_balanced_folds(
    X,
    y,
    requested_folds = inner_folds_requested,
    return_train = TRUE
  )

  control <- caret::trainControl(
    method = "cv",
    number = length(inner_index),
    index = inner_index,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    allowParallel = FALSE
  )

  fit <- caret::train(
    x = X,
    y = y,
    method = "glmboost",
    metric = "ROC",
    trControl = control,
    preProcess = c("center", "scale"),
    tuneLength = max(3, tuning_length)
  )

  importance <- caret::varImp(
    fit,
    scale = FALSE
  )$importance

  if ("Overall" %in% colnames(importance)) {
    score <- importance$Overall
  } else {
    score <- rowMeans(importance, na.rm = TRUE)
  }

  names(score) <- rownames(importance)
  score <- sort(score, decreasing = TRUE)
  score <- score[is.finite(score)]

  positive <- score[score > 0]

  if (length(positive)) {
    cutoff <- stats::median(positive)
    selected <- names(positive)[positive >= cutoff]
  } else {
    selected <- character(0)
  }

  maximum_selected <- min(10, ncol(X))

  if (length(selected) > maximum_selected) {
    selected <- names(score)[seq_len(maximum_selected)]
  }

  if (length(selected) < 2) {
    selected <- names(score)[seq_len(min(2, length(score)))]
  }

  if (!length(selected)) {
    selected <- univariate_fallback(X, y, number = 2)
  }

  intersect(selected, colnames(X))
}

fit_stepwise_glm <- function(
  X,
  y,
  direction = "both"
) {
  means <- vapply(X, mean, numeric(1), na.rm = TRUE)
  standard_deviations <- vapply(
    X,
    stats::sd,
    numeric(1),
    na.rm = TRUE
  )

  standard_deviations[
    !is.finite(standard_deviations) |
      standard_deviations == 0
  ] <- 1

  X_scaled <- sweep(
    sweep(
      as.matrix(X),
      2,
      means,
      "-"
    ),
    2,
    standard_deviations,
    "/"
  )

  model_data <- data.frame(
    Outcome = as.numeric(y == case_label),
    X_scaled,
    check.names = FALSE
  )

  predictor_names <- colnames(X_scaled)

  full_formula <- reformulate(
    predictor_names,
    response = "Outcome"
  )

  full_model <- suppressWarnings(
    glm(
      full_formula,
      data = model_data,
      family = binomial(),
      control = glm.control(maxit = 100)
    )
  )

  if (direction == "forward") {
    null_model <- glm(
      Outcome ~ 1,
      data = model_data,
      family = binomial()
    )

    selected_model <- suppressWarnings(
      MASS::stepAIC(
        null_model,
        scope = list(
          lower = ~1,
          upper = full_formula
        ),
        direction = "forward",
        trace = FALSE
      )
    )
  } else {
    selected_model <- suppressWarnings(
      MASS::stepAIC(
        full_model,
        direction = direction,
        trace = FALSE
      )
    )
  }

  list(
    model = selected_model,
    center = means,
    scale = standard_deviations,
    selected = attr(
      terms(selected_model),
      "term.labels"
    )
  )
}

select_stepglm_features <- function(X, y, direction) {
  fit <- fit_stepwise_glm(
    X,
    y,
    direction = direction
  )

  selected <- intersect(
    fit$selected,
    colnames(X)
  )

  if (!length(selected)) {
    selected <- univariate_fallback(X, y, number = 2)
  }

  selected
}

select_features <- function(selector, X, y, seed) {
  selector <- canonical_token(selector)

  if (selector == "None") {
    return(colnames(X))
  }

  if (selector == "Lasso") {
    return(select_lasso_features(X, y, seed))
  }

  if (selector == "RF") {
    return(select_rf_features(X, y, seed))
  }

  if (selector == "glmBoost") {
    return(select_glmboost_features(X, y, seed))
  }

  if (selector == "Stepglm[both]") {
    return(select_stepglm_features(X, y, "both"))
  }

  if (selector == "Stepglm[backward]") {
    return(select_stepglm_features(X, y, "backward"))
  }

  stop("Unsupported feature selector: ", selector)
}

# ============================================================
# 6. Classifier fitting and prediction
# ============================================================

fit_caret_classifier <- function(
  classifier,
  X,
  y,
  seed
) {
  set.seed(seed)

  if (ncol(X) < 1) {
    stop("Classifier received zero predictors.")
  }

  zero_variance <- caret::nearZeroVar(
    X,
    saveMetrics = FALSE
  )

  if (length(zero_variance)) {
    X <- X[, -zero_variance, drop = FALSE]
  }

  if (ncol(X) < 1) {
    stop("All predictors had near-zero variance.")
  }

  inner_index <- make_balanced_folds(
    X,
    y,
    requested_folds = inner_folds_requested,
    return_train = TRUE
  )

  train_control <- caret::trainControl(
    method = "cv",
    number = length(inner_index),
    index = inner_index,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    allowParallel = FALSE
  )

  classifier <- canonical_token(classifier)
  lower <- tolower(classifier)

  tree_models <- c("rf", "gbm", "xgboost")
  preprocessing <- if (lower %in% tree_models) {
    NULL
  } else {
    c("center", "scale")
  }

  common_arguments <- list(
    x = X,
    y = y,
    metric = "ROC",
    trControl = train_control,
    preProcess = preprocessing
  )

  if (classifier == "SVM") {
    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "svmRadial",
          tuneLength = tuning_length
        )
      )
    )
  } else if (classifier == "Ridge") {
    tuning_grid <- expand.grid(
      alpha = 0,
      lambda = 10^seq(-4, 1, length.out = 30)
    )

    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "glmnet",
          tuneGrid = tuning_grid
        )
      )
    )
  } else if (classifier == "Lasso") {
    tuning_grid <- expand.grid(
      alpha = 1,
      lambda = 10^seq(-4, 1, length.out = 30)
    )

    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "glmnet",
          tuneGrid = tuning_grid
        )
      )
    )
  } else if (grepl("^Enet\\[alpha=", classifier)) {
    alpha_value <- as.numeric(
      sub(
        "^Enet\\[alpha=([0-9.]+)\\]$",
        "\\1",
        classifier
      )
    )

    tuning_grid <- expand.grid(
      alpha = alpha_value,
      lambda = 10^seq(-4, 1, length.out = 30)
    )

    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "glmnet",
          tuneGrid = tuning_grid
        )
      )
    )
  } else if (classifier == "glmBoost") {
    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "glmboost",
          tuneLength = max(3, tuning_length)
        )
      )
    )
  } else if (classifier == "plsRglm") {
    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "plsRglm",
          tuneLength = min(
            max(2, tuning_length),
            ncol(X)
          )
        )
      )
    )
  } else if (classifier == "RF") {
    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "ranger",
          tuneLength = tuning_length,
          importance = "permutation",
          num.trees = 1000
        )
      )
    )
  } else if (classifier == "GBM") {
    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "gbm",
          tuneLength = tuning_length,
          verbose = FALSE
        )
      )
    )
  } else if (classifier == "LDA") {
    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(method = "lda")
      )
    )
  } else if (classifier == "XGBoost") {
    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "xgbTree",
          tuneLength = tuning_length,
          verbose = FALSE
        )
      )
    )
  } else if (classifier == "NaiveBayes") {
    fit <- do.call(
      caret::train,
      c(
        common_arguments,
        list(
          method = "naive_bayes",
          tuneLength = tuning_length
        )
      )
    )
  } else {
    stop("Unsupported caret classifier: ", classifier)
  }

  list(
    engine = "caret",
    classifier = classifier,
    features = colnames(X),
    model = fit
  )
}

fit_stepglm_classifier <- function(
  classifier,
  X,
  y
) {
  classifier <- canonical_token(classifier)

  direction <- switch(
    classifier,
    "Stepglm[both]" = "both",
    "Stepglm[backward]" = "backward",
    "Stepglm[forward]" = "forward",
    stop("Unsupported stepwise classifier: ", classifier)
  )

  fit <- fit_stepwise_glm(
    X,
    y,
    direction = direction
  )

  list(
    engine = "stepglm",
    classifier = classifier,
    features = colnames(X),
    model = fit$model,
    center = fit$center,
    scale = fit$scale
  )
}

fit_classifier <- function(
  classifier,
  X,
  y,
  seed
) {
  classifier <- canonical_token(classifier)

  if (grepl("^Stepglm\\[", classifier)) {
    fit_stepglm_classifier(
      classifier,
      X,
      y
    )
  } else {
    fit_caret_classifier(
      classifier,
      X,
      y,
      seed
    )
  }
}

predict_classifier <- function(
  fitted_model,
  new_X
) {
  required_features <- fitted_model$features

  missing_features <- setdiff(
    required_features,
    colnames(new_X)
  )

  if (length(missing_features)) {
    stop(
      "Prediction data are missing features: ",
      paste(missing_features, collapse = ", ")
    )
  }

  new_X <- new_X[
    ,
    required_features,
    drop = FALSE
  ]

  if (fitted_model$engine == "caret") {
    probability <- predict(
      fitted_model$model,
      newdata = new_X,
      type = "prob"
    )[, case_label]
  } else if (fitted_model$engine == "stepglm") {
    scaled <- sweep(
      sweep(
        as.matrix(new_X),
        2,
        fitted_model$center[required_features],
        "-"
      ),
      2,
      fitted_model$scale[required_features],
      "/"
    )

    probability <- suppressWarnings(
      predict(
        fitted_model$model,
        newdata = data.frame(
          scaled,
          check.names = FALSE
        ),
        type = "response"
      )
    )
  } else {
    stop("Unknown model engine.")
  }

  probability <- as.numeric(probability)
  probability <- pmin(pmax(probability, 1e-6), 1 - 1e-6)
  probability
}


# -----------------------------------------------------------------------------
# Safer grouped fold construction for generic cohorts.
# This definition intentionally overrides the inherited make_balanced_folds().
# -----------------------------------------------------------------------------
make_balanced_folds <- function(X, y, requested_folds, return_train = FALSE) {
  dataset_id <- dataset_from_rownames(X)
  group_id <- group_from_rownames(X)
  unique_groups <- unique(group_id)

  group_dataset <- vapply(unique_groups, function(g) unique(dataset_id[group_id == g])[1], character(1))
  group_outcome <- vapply(unique_groups, function(g) {
    observed <- unique(as.character(y[group_id == g]))
    if (length(observed) > 1L) "Paired" else observed[1]
  }, character(1))

  if (all(group_outcome == "Paired")) {
    k <- min(as.integer(requested_folds), length(unique_groups))
    if (k < 2L) stop("Too few biological groups for cross-validation.")
    set.seed(1L)
    shuffled <- sample(unique_groups, length(unique_groups), replace = FALSE)
    group_folds <- split(shuffled, rep(seq_len(k), length.out = length(shuffled)))
  } else {
    strata <- interaction(group_dataset, group_outcome, drop = TRUE, lex.order = TRUE)
    strata_counts <- table(strata)
    k <- min(as.integer(requested_folds), min(strata_counts))
    if (k < 2L) {
      strata <- factor(group_outcome)
      strata_counts <- table(strata)
      k <- min(as.integer(requested_folds), min(strata_counts))
    }
    if (k < 2L) stop("Too few independent groups per class for cross-validation.")
    group_folds_idx <- caret::createFolds(strata, k = k, returnTrain = FALSE)
    group_folds <- lapply(group_folds_idx, function(idx) unique_groups[idx])
  }

  sample_assessment_folds <- lapply(group_folds, function(groups) which(group_id %in% groups))
  if (return_train) {
    all_idx <- seq_len(nrow(X))
    return(lapply(sample_assessment_folds, function(test_idx) setdiff(all_idx, test_idx)))
  }
  sample_assessment_folds
}

clean_gene_symbol_tdrx <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NaN", "---", "--", "NULL", ".")] <- NA_character_
  x <- sub("[[:space:]]*(///|//|;|\\|).*$", "", x)
  toupper(trimws(x))
}

load_stage1_expression <- function(expression_file, metadata_file, candidate_gene_file,
                                   dataset_name, sample_col = "Samples",
                                   condition_col = "Condition", pair_col = "Pair_ID") {
  expr_raw <- data.table::fread(expression_file, data.table = FALSE, check.names = FALSE)
  md <- data.table::fread(metadata_file, data.table = FALSE, check.names = FALSE)
  cand <- data.table::fread(candidate_gene_file, data.table = FALSE, check.names = FALSE)

  if (ncol(expr_raw) < 3L) stop(dataset_name, ": expression file must have probe_id, gene_symbol, and samples.")
  if (!all(c(sample_col, condition_col) %in% colnames(md))) {
    stop(dataset_name, ": metadata must contain ", sample_col, " and ", condition_col, ".")
  }

  gene_col <- if ("Gene" %in% colnames(cand)) "Gene" else colnames(cand)[1]
  candidates <- unique(clean_gene_symbol_tdrx(cand[[gene_col]]))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  gene_symbols <- clean_gene_symbol_tdrx(expr_raw[[2]])
  sample_frame <- expr_raw[, -(1:2), drop = FALSE]
  sample_ids <- intersect(colnames(sample_frame), as.character(md[[sample_col]]))
  if (length(sample_ids) < 4L) stop(dataset_name, ": fewer than four matched samples.")

  sample_frame <- sample_frame[, sample_ids, drop = FALSE]
  X_probe <- as.matrix(data.frame(lapply(sample_frame, as.numeric), check.names = FALSE))
  colnames(X_probe) <- sample_ids
  keep <- !is.na(gene_symbols) & gene_symbols %in% candidates & apply(X_probe, 1, function(z) all(is.finite(z)))
  if (sum(keep) < 2L) stop(dataset_name, ": fewer than two candidate genes available.")

  gene_df <- data.frame(gene_symbol = gene_symbols[keep], X_probe[keep, , drop = FALSE], check.names = FALSE)
  gene_split <- split(seq_len(nrow(gene_df)), gene_df$gene_symbol)
  gene_mat <- do.call(rbind, lapply(gene_split, function(idx) colMeans(as.matrix(gene_df[idx, -1, drop = FALSE]))))
  rownames(gene_mat) <- names(gene_split)
  gene_mat <- gene_mat[intersect(candidates, rownames(gene_mat)), , drop = FALSE]

  variances <- matrixStats::rowVars(gene_mat)
  gene_mat <- gene_mat[is.finite(variances) & variances > 1e-10, , drop = FALSE]
  if (nrow(gene_mat) < 2L) stop(dataset_name, ": fewer than two variable candidate genes remain.")

  md <- md[match(sample_ids, md[[sample_col]]), , drop = FALSE]
  cond <- trimws(as.character(md[[condition_col]]))
  if (!all(cond %in% c(control_label, case_label))) {
    stop(dataset_name, ": Condition must be coded exactly as Control/Case in Stage 1 metadata.")
  }
  y <- factor(cond, levels = c(case_label, control_label))

  feature_map <- data.frame(
    Gene_Symbol = rownames(gene_mat),
    Feature_Name = make.names(rownames(gene_mat), unique = TRUE),
    stringsAsFactors = FALSE
  )
  X <- as.data.frame(t(gene_mat), check.names = FALSE)
  colnames(X) <- feature_map$Feature_Name

  if (pair_col %in% colnames(md) && any(!is.na(md[[pair_col]]) & nzchar(trimws(as.character(md[[pair_col]]))))) {
    group <- as.character(md[[pair_col]])
    bad <- is.na(group) | !nzchar(trimws(group))
    group[bad] <- as.character(md[[sample_col]][bad])
  } else {
    group <- as.character(md[[sample_col]])
  }
  rownames(X) <- paste(dataset_name, group, md[[sample_col]], sep = "__")

  list(X = X, y = y, metadata = md, feature_map = feature_map,
       candidates_requested = candidates, genes_used = feature_map$Gene_Symbol)
}

safe_auc_internal <- function(y, probability) {
  tryCatch(as.numeric(pROC::auc(compute_auc_object(y, probability))), error = function(e) NA_real_)
}

run_pipeline_internal <- function(model_index, model_name, training_X, training_y,
                                  checkpoint_dir, global_seed) {
  checkpoint_file <- file.path(checkpoint_dir, sprintf("pipeline_%03d.rds", model_index))
  if (file.exists(checkpoint_file)) {
    cached <- tryCatch(readRDS(checkpoint_file), error = function(e) NULL)
    if (is.list(cached) && identical(cached$model, model_name)) return(cached)
  }

  started <- Sys.time()
  result <- tryCatch({
    parsed <- parse_pipeline(model_name)
    set.seed(global_seed + model_index)
    outer_folds <- make_balanced_folds(training_X, training_y, outer_folds_requested, return_train = FALSE)
    oof <- rep(NA_real_, nrow(training_X))
    fold_selection <- vector("list", length(outer_folds))

    for (fold_index in seq_along(outer_folds)) {
      assess <- outer_folds[[fold_index]]
      analyse <- setdiff(seq_len(nrow(training_X)), assess)
      Xa <- training_X[analyse, , drop = FALSE]
      ya <- training_y[analyse]
      Xt <- training_X[assess, , drop = FALSE]

      selected <- select_features(parsed$selector, Xa, ya,
                                  seed = global_seed + model_index * 1000L + fold_index)
      selected <- intersect(selected, colnames(Xa))
      if (!length(selected)) stop("No features selected in outer fold ", fold_index)
      fold_selection[[fold_index]] <- selected

      fitted <- fit_classifier(parsed$classifier, Xa[, selected, drop = FALSE], ya,
                               seed = global_seed + model_index * 10000L + fold_index)
      oof[assess] <- predict_classifier(fitted, Xt[, selected, drop = FALSE])
    }

    if (anyNA(oof)) stop("OOF predictions contain missing values.")
    auc <- safe_auc_internal(training_y, oof)
    threshold <- tryCatch(derive_training_threshold(training_y, oof), error = function(e) 0.5)

    sel_tab <- sort(table(unlist(fold_selection)), decreasing = TRUE)
    selection_table <- data.frame(
      Model_ID = model_index,
      Model = model_name,
      Feature = names(sel_tab),
      Outer_Folds_Selected = as.integer(sel_tab),
      Selection_Frequency = as.integer(sel_tab) / length(outer_folds),
      stringsAsFactors = FALSE
    )

    full_selected <- select_features(parsed$selector, training_X, training_y,
                                     seed = global_seed + model_index * 100000L)
    full_selected <- intersect(full_selected, colnames(training_X))

    list(status = "OK", error = NA_character_, model_id = model_index, model = model_name,
         selector = parsed$selector, classifier = parsed$classifier,
         selected_features = full_selected, selection_table = selection_table,
         internal_probability = oof, internal_auc = auc, threshold = threshold,
         outer_folds = length(outer_folds), elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")))
  }, error = function(e) {
    list(status = "FAILED", error = conditionMessage(e), model_id = model_index, model = model_name,
         selector = NA_character_, classifier = NA_character_, selected_features = character(),
         selection_table = data.frame(), internal_probability = rep(NA_real_, nrow(training_X)),
         internal_auc = NA_real_, threshold = NA_real_, outer_folds = NA_integer_,
         elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")))
  })
  saveRDS(result, checkpoint_file)
  result
}

run_113ml_internal <- function(expression_file, metadata_file, candidate_gene_file,
                               dataset_name, output_dir,
                               outer_folds = 5L, inner_folds = 4L,
                               tuning_len = 3L, workers = 1L,
                               consensus_top_models = 20L,
                               seed = 20260908L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint_dir <- file.path(output_dir, "checkpoints")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  outer_folds_requested <<- as.integer(outer_folds)
  inner_folds_requested <<- as.integer(inner_folds)
  tuning_length <<- as.integer(tuning_len)

  dat <- load_stage1_expression(expression_file, metadata_file, candidate_gene_file, dataset_name)
  training_X <- dat$X
  training_y <- dat$y
  feature_map <- dat$feature_map

  if (workers > 1L) future::plan(future::multisession, workers = workers) else future::plan(future::sequential)
  pipeline_results <- future.apply::future_lapply(seq_along(model_definitions), function(i) {
    run_pipeline_internal(i, model_definitions[i], training_X, training_y, checkpoint_dir, seed)
  }, future.seed = TRUE)
  future::plan(future::sequential)
  saveRDS(pipeline_results, file.path(output_dir, "all_pipeline_results.rds"))

  performance <- dplyr::bind_rows(lapply(pipeline_results, function(r) {
    selected_original <- feature_map$Gene_Symbol[match(r$selected_features, feature_map$Feature_Name)]
    data.frame(Model_ID = r$model_id, Model = r$model, Selector = r$selector, Classifier = r$classifier,
               Internal_AUC = r$internal_auc, Selected_Genes_n = length(r$selected_features),
               Selected_Genes = paste(na.omit(selected_original), collapse = ", "),
               Status = r$status, Error = r$error, Runtime_Minutes = r$elapsed_minutes,
               stringsAsFactors = FALSE)
  })) %>% arrange(desc(Status == "OK"), desc(Internal_AUC), Selected_Genes_n, Model_ID)
  write.csv(performance, file.path(output_dir, "all_model_performance.csv"), row.names = FALSE)

  eligible <- performance %>% filter(Status == "OK", is.finite(Internal_AUC))
  if (!nrow(eligible)) stop(dataset_name, ": no successful 113-ML pipeline.")
  best <- eligible[1, , drop = FALSE]

  selection <- dplyr::bind_rows(lapply(pipeline_results, function(r) r$selection_table))
  if (nrow(selection)) {
    selection$Gene_Symbol <- feature_map$Gene_Symbol[match(selection$Feature, feature_map$Feature_Name)]
  }
  write.csv(selection, file.path(output_dir, "feature_selection_stability.csv"), row.names = FALSE)

  selector_top <- eligible %>% filter(!is.na(Selector), Selector != "None") %>% slice_head(n = consensus_top_models)
  genes <- feature_map$Gene_Symbol

  if (nrow(selector_top)) {
    model_weights <- pmax(selector_top$Internal_AUC - 0.5, 0)
    if (sum(model_weights) <= 0) model_weights <- rep(1, nrow(selector_top))
    model_weights <- model_weights / sum(model_weights)

    score_mat <- matrix(0, nrow = length(genes), ncol = nrow(selector_top),
                        dimnames = list(genes, selector_top$Model_ID))
    for (j in seq_len(nrow(selector_top))) {
      mid <- selector_top$Model_ID[j]
      sub <- selection %>% filter(Model_ID == mid)
      if (nrow(sub)) {
        idx <- match(sub$Gene_Symbol, genes)
        ok <- !is.na(idx)
        score_mat[idx[ok], j] <- sub$Selection_Frequency[ok]
      }
    }
    s_ml <- as.numeric(score_mat %*% model_weights)
    support_count <- rowSums(score_mat > 0)
    mean_freq <- rowMeans(score_mat)
  } else {
    best_r <- pipeline_results[[best$Model_ID]]
    best_genes <- feature_map$Gene_Symbol[match(best_r$selected_features, feature_map$Feature_Name)]
    s_ml <- as.numeric(genes %in% best_genes)
    support_count <- as.integer(genes %in% best_genes)
    mean_freq <- s_ml
  }

  best_result <- pipeline_results[[best$Model_ID]]
  best_genes <- feature_map$Gene_Symbol[match(best_result$selected_features, feature_map$Feature_Name)]
  best_gene_table <- data.frame(Gene = na.omit(best_genes), stringsAsFactors = FALSE)
  write.csv(best_gene_table, file.path(output_dir, "best_model_selected_genes.csv"), row.names = FALSE)

  gene_scores <- data.frame(
    Gene = genes,
    S_ML_Disease = pmin(pmax(s_ml, 0), 1),
    ML_TopModelSupportCount = support_count,
    ML_MeanSelectionFrequency = mean_freq,
    BestModelSelected = genes %in% best_genes,
    Best_Model_ID = best$Model_ID,
    Best_Model = best$Model,
    Best_Internal_AUC = best$Internal_AUC,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(S_ML_Disease), desc(ML_TopModelSupportCount), Gene)
  write.csv(gene_scores, file.path(output_dir, "gene_ML_scores.csv"), row.names = FALSE)

  report <- c(
    paste0("Dataset: ", dataset_name),
    paste0("Candidate genes requested: ", length(dat$candidates_requested)),
    paste0("Candidate genes used: ", length(dat$genes_used)),
    paste0("Successful pipelines: ", sum(performance$Status == "OK"), "/", length(model_definitions)),
    paste0("Best internal model: ", best$Model, " (OOF AUC=", sprintf("%.4f", best$Internal_AUC), ")"),
    paste0("Top selector-based models used for gene consensus: ", nrow(selector_top)),
    "External validation was intentionally not used in this Stage 1 ML engine."
  )
  writeLines(report, file.path(output_dir, "ML_internal_report.txt"))

  list(performance = performance, best_model = best, selection = selection,
       gene_scores = gene_scores, feature_map = feature_map)
}
