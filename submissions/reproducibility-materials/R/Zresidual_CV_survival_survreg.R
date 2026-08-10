#' Cross-validated Z-residuals for parametric survival regression models
#'
#' Compute cross-validated Z-residuals for a `survival::survreg` model. For
#' each fold, the parametric survival model is refitted using the training
#' data, and `Zresidual()` is called on the held-out test data. The fold-wise
#' residuals are combined into a single `cvzresid` object.
#'
#' This function is the parametric survival regression backend used by
#' `Zresidual_CV()`. Users usually call `Zresidual_CV()` directly rather than
#' calling this backend.
#'
#' @param object A fitted `survival::survreg` object.
#' @param data A data frame containing the variables used in `object`.
#' @param nfolds Number of cross-validation folds.
#' @param foldlist Optional list of fold indices. Each element should contain
#'   row indices for one held-out test fold, relative to the complete-case data
#'   used internally.
#' @param nrep Number of randomized Z-residual replicates.
#' @param randomized Logical. If `TRUE`, randomized residuals are generated.
#' @param type Optional residual or model subtype passed to `Zresidual()`.
#' @param log_pointpred Optional prediction backend function.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to `Zresidual()`.
#'
#' @keywords internal
Zresidual_CV_survival_survreg <- function(object,
                                          data,
                                          nfolds = NULL,
                                          foldlist = NULL,
                                          nrep = 1,
                                          randomized = TRUE,
                                          type = NULL,
                                          log_pointpred = NULL,
                                          seed = NULL,
                                          ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (missing(object) || is.null(object)) {
    stop("Zresidual_CV_survival_survreg: `object` must be provided.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("Zresidual_CV_survival_survreg: `data` must be provided.", call. = FALSE)
  }
  
  nrep <- as.integer(nrep)
  
  if (length(nrep) != 1L || is.na(nrep) || nrep < 1L) {
    stop(
      "Zresidual_CV_survival_survreg: `nrep` must be a positive integer.",
      call. = FALSE
    )
  }
  
  nrep_eff <- if (isTRUE(randomized)) nrep else 1L
  
  mf_all <- stats::model.frame(
    object$terms,
    data = data,
    na.action = stats::na.pass
  )
  
  keep <- stats::complete.cases(mf_all)
  
  if (!all(keep)) {
    warning(
      sprintf(
        "Zresidual_CV_survival_survreg: removed %d rows with missing values before cross-validation.",
        sum(!keep)
      ),
      call. = FALSE
    )
  }
  
  data_used <- data[keep, , drop = FALSE]
  
  mf <- stats::model.frame(
    object$terms,
    data = data_used,
    na.action = stats::na.pass
  )
  
  y <- stats::model.response(mf)
  ymat <- as.matrix(y)
  
  if (ncol(ymat) < 2L) {
    stop(
      "Zresidual_CV_survival_survreg: the survival response must contain a status column.",
      call. = FALSE
    )
  }
  
  status <- as.integer(ymat[, ncol(ymat)])
  time <- if (ncol(ymat) == 2L) ymat[, 1L] else ymat[, 2L]
  covariates <- mf[, -1L, drop = FALSE]
  
  n <- nrow(data_used)
  
  if (n < 2L) {
    stop(
      "Zresidual_CV_survival_survreg: at least two complete observations are required.",
      call. = FALSE
    )
  }
  
  if (is.null(foldlist)) {
    if (is.null(nfolds)) {
      nfolds <- min(10L, n)
    }
    
    nfolds <- as.integer(nfolds)
    
    if (length(nfolds) != 1L || is.na(nfolds) || nfolds < 2L) {
      stop(
        "Zresidual_CV_survival_survreg: `nfolds` must be an integer >= 2.",
        call. = FALSE
      )
    }
    
    if (nfolds > n) {
      stop(
        "Zresidual_CV_survival_survreg: `nfolds` cannot exceed the number of observations.",
        call. = FALSE
      )
    }
    
    foldlist <- make_fold(
      fix_var = covariates,
      y = y,
      k = nfolds,
      censor = status
    )
  } else {
    foldlist <- lapply(foldlist, function(x) sort(unique(as.integer(x))))
    flat <- unlist(foldlist, use.names = FALSE)
    
    if (!identical(sort(flat), seq_len(n))) {
      stop(
        paste0(
          "Zresidual_CV_survival_survreg: `foldlist` must be a non-overlapping ",
          "partition of 1:nrow(data_used), where data_used is the complete-case data."
        ),
        call. = FALSE
      )
    }
  }
  
  z_full <- matrix(NA_real_, nrow = n, ncol = nrep_eff)
  rsp_full <- matrix(NA_real_, nrow = n, ncol = nrep_eff)
  lp_full <- rep(NA_real_, n)
  surv_prob_full <- rep(NA_real_, n)
  
  extra_args <- list(...)
  failed_folds <- integer(0)
  
  for (fid in seq_along(foldlist)) {
    test_idx <- foldlist[[fid]]
    
    train_data <- data_used[-test_idx, , drop = FALSE]
    test_data <- data_used[test_idx, , drop = FALSE]
    
    fit_train <- tryCatch(
      suppressWarnings(stats::update(object, data = train_data)),
      error = function(e) {
        warning(
          sprintf(
            "Zresidual_CV_survival_survreg: fold %s refitting failed: %s",
            fid,
            conditionMessage(e)
          ),
          call. = FALSE
        )
        NULL
      }
    )
    
    if (is.null(fit_train)) {
      failed_folds <- c(failed_folds, fid)
      next
    }
    
    z_call <- c(
      list(
        fit = fit_train,
        data = test_data,
        log_pointpred = log_pointpred,
        type = type,
        randomized = randomized,
        nrep = nrep_eff
      ),
      extra_args
    )
    
    z_fold <- tryCatch(
      do.call(Zresidual, z_call),
      error = function(e) {
        warning(
          sprintf(
            "Zresidual_CV_survival_survreg: fold %s Zresidual calculation failed: %s",
            fid,
            conditionMessage(e)
          ),
          call. = FALSE
        )
        NULL
      }
    )
    
    if (is.null(z_fold)) {
      failed_folds <- c(failed_folds, fid)
      next
    }
    
    z_mat <- as.matrix(z_fold)
    
    if (ncol(z_mat) != nrep_eff) {
      stop(
        sprintf(
          "Zresidual_CV_survival_survreg: fold %s returned %d columns, expected %d.",
          fid,
          ncol(z_mat),
          nrep_eff
        ),
        call. = FALSE
      )
    }
    
    if (nrow(z_mat) != length(test_idx)) {
      stop(
        sprintf(
          "Zresidual_CV_survival_survreg: fold %s returned %d rows, expected %d.",
          fid,
          nrow(z_mat),
          length(test_idx)
        ),
        call. = FALSE
      )
    }
    
    z_full[test_idx, ] <- z_mat
    
    rsp_fold <- attr(z_fold, "rsp")
    
    if (!is.null(rsp_fold)) {
      rsp_mat <- as.matrix(rsp_fold)
      
      if (nrow(rsp_mat) == length(test_idx)) {
        rsp_full[test_idx, ] <- rsp_mat
      }
    }
    
    lp_fold <- attr(z_fold, "linear_pred")
    
    if (is.null(lp_fold)) {
      lp_fold <- attr(z_fold, "linear.pred")
    }
    
    if (!is.null(lp_fold)) {
      lp_full[test_idx] <- as.numeric(lp_fold)
    }
    
    surv_fold <- attr(z_fold, "Survival.Prob")
    
    if (!is.null(surv_fold)) {
      surv_prob_full[test_idx] <- as.numeric(surv_fold)
    }
  }
  
  if (length(failed_folds) > 0L) {
    stop(
      sprintf(
        "Zresidual_CV_survival_survreg: failed fold(s): %s.",
        paste(unique(failed_folds), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (anyNA(z_full)) {
    stop(
      "Zresidual_CV_survival_survreg: combined CV residual matrix contains NA values.",
      call. = FALSE
    )
  }
  
  colnames(z_full) <- paste0("rep", seq_len(nrep_eff))
  class(z_full) <- c("cvzresid", "zresid", "matrix", "array")
  
  attr(z_full, "rsp") <- rsp_full
  attr(z_full, "type") <- "survival"
  attr(z_full, "Survival.Prob") <- surv_prob_full
  attr(z_full, "linear.pred") <- lp_full
  attr(z_full, "linear_pred") <- lp_full
  attr(z_full, "censored.status") <- as.integer(status == 0L)
  attr(z_full, "covariates") <- covariates
  attr(z_full, "object.model.frame") <- mf
  attr(z_full, "foldlist") <- foldlist
  attr(z_full, "original_row_index") <- which(keep)
  
  attr(z_full, "zcov") <- list(
    type = type,
    family = "survreg",
    response_name = names(mf)[1L],
    response = y,
    covariates = covariates,
    linear_pred = lp_full,
    obs_id = seq_len(n),
    y_type = ifelse(status == 1L, 1L, 0L),
    y_type_kind = "censor",
    y_type_levels = c(censored = 0L, event = 1L),
    extra = list(
      time = time,
      status = status,
      Survival.Prob = surv_prob_full,
      cv = TRUE,
      foldlist = foldlist,
      original_row_index = which(keep)
    )
  )
  
  z_full
}