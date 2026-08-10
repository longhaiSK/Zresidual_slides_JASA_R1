#' Cross-validated Z-residuals
#'
#' Compute cross-validated Z-residuals by dispatching to the appropriate
#' model-specific backend. For each cross-validation fold, the model is
#' refitted using the training data, and `Zresidual()` is then called on the
#' held-out test data. The fold-wise residuals are combined into one
#' `cvzresid` object.
#'
#' This function is currently intended for models fitted with the
#' `survival` package, including `coxph` and `survreg` models. The function
#' does not perform Bayesian posterior summarization; it only implements
#' cross-validation for frequentist survival models.
#'
#' @param object A fitted model object. Currently supported classes are
#'   `"coxph"` and `"survreg"`.
#' @param data A data frame containing the variables used to fit `object`.
#'   This must be provided explicitly.
#' @param nfolds Number of cross-validation folds. If `foldlist` is not
#'   supplied and `nfolds` is `NULL`, the default is `min(10, nrow(data))`.
#' @param foldlist Optional list of fold indices. Each element should contain
#'   the row indices for one test fold. The folds must form a non-overlapping
#'   partition of `seq_len(nrow(data_used))`, where `data_used` is the complete
#'   case subset used internally.
#' @param nrep Number of randomized Z-residual replicates to generate for each
#'   observation. Default is `1`.
#' @param randomized Logical. If `TRUE`, randomized Z-residuals are generated
#'   for censored observations. If `FALSE`, only one non-randomized residual is
#'   returned.
#' @param type Optional character string specifying a model subtype or residual
#'   type. For example, `type = "frailty"` can be used for `coxph` models with
#'   frailty terms when the frailty-specific prediction backend is required.
#' @param log_pointpred Optional prediction backend function. If `NULL`, the
#'   backend is selected automatically by the model-specific CV function.
#' @param seed Optional random seed used before fold generation and randomized
#'   residual calculation.
#' @param ... Additional arguments passed to the model-specific CV backend and
#'   ultimately to `Zresidual()`.
#'
#' @return A matrix-like object of class `c("cvzresid", "zresid", "matrix",
#'   "array")`, with one row per observation and one column per residual
#'   replicate. The object contains attributes used by plotting and diagnostic
#'   functions, including `foldlist`, `rsp`, `covariates`, `linear_pred`, and
#'   `zcov`.
#'
#' @details
#' `Zresidual_CV()` is a top-level wrapper. It does not calculate residuals
#' directly. Instead, it identifies the model class and calls a backend such as
#' `Zresidual_CV_survival_coxph()` or `Zresidual_CV_survival_survreg()`.
#'
#' The intended workflow is:
#'
#' 1. split the data into folds;
#' 2. refit the model on the training data for each fold;
#' 3. call `Zresidual()` on the held-out test data;
#' 4. combine the fold-wise residuals into one object.
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' data(lung, package = "survival")
#' lung2 <- na.omit(lung[, c("time", "status", "age", "sex", "ph.ecog")])
#'
#' fit <- coxph(Surv(time, status == 2) ~ age + sex + ph.ecog, data = lung2)
#'
#' zcv <- Zresidual_CV(
#'   object = fit,
#'   data = lung2,
#'   nfolds = 5,
#'   nrep = 10,
#'   seed = 123
#' )
#'
#' plot(zcv)
#' qqnorm(zcv)
#' }
#'
#' @export
Zresidual_CV <- function(object,
                         data,
                         nfolds = NULL,
                         foldlist = NULL,
                         nrep = 1,
                         randomized = TRUE,
                         type = NULL,
                         log_pointpred = NULL,
                         seed = NULL,
                         ...) {
  if (missing(object) || is.null(object)) {
    stop("Zresidual_CV: `object` must be provided.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("Zresidual_CV: `data` must be provided.", call. = FALSE)
  }
  
  if (inherits(object, "coxph")) {
    return(
      Zresidual_CV_survival_coxph(
        object = object,
        data = data,
        nfolds = nfolds,
        foldlist = foldlist,
        nrep = nrep,
        randomized = randomized,
        type = type,
        log_pointpred = log_pointpred,
        seed = seed,
        ...
      )
    )
  }
  
  if (inherits(object, "survreg")) {
    return(
      Zresidual_CV_survival_survreg(
        object = object,
        data = data,
        nfolds = nfolds,
        foldlist = foldlist,
        nrep = nrep,
        randomized = randomized,
        type = type,
        log_pointpred = log_pointpred,
        seed = seed,
        ...
      )
    )
  }
  
  stop(
    sprintf(
      "Zresidual_CV: unsupported model class: %s.",
      paste(class(object), collapse = ", ")
    ),
    call. = FALSE
  )
}