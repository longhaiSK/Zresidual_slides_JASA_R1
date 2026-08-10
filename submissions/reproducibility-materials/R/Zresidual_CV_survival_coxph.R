#' Cross-validated Z-residuals for Cox proportional hazards models
#'
#' Compute cross-validated Z-residuals for a `survival::coxph` model. For each
#' fold, the Cox model is refitted on the training data, and `Zresidual()` is
#' called on the held-out test data. The resulting fold-wise residuals are then
#' combined into a single `cvzresid` object.
#'
#' This function is the survival Cox backend used by `Zresidual_CV()`. Users
#' usually call `Zresidual_CV()` directly rather than calling this backend.
#'
#' @param object A fitted `survival::coxph` object.
#' @param data A data frame containing the variables used in `object`.
#' @param nfolds Number of cross-validation folds.
#' @param foldlist Optional list of fold indices. Each element should contain
#'   row indices for one held-out test fold, relative to the complete-case data
#'   used internally.
#' @param nrep Number of randomized Z-residual replicates.
#' @param randomized Logical. If `TRUE`, randomized residuals are generated.
#' @param type Optional residual or model subtype. Use `type = "frailty"` for
#'   frailty Cox models if needed.
#' @param log_pointpred Optional prediction backend function.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to `Zresidual()`.
#'
#' @keywords internal
Zresidual_CV_survival_coxph <- function(object,
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
    stop("Zresidual_CV_survival_coxph: `object` must be provided.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("Zresidual_CV_survival_coxph: `data` must be provided.", call. = FALSE)
  }
  
  nrep <- as.integer(nrep)
  
  if (length(nrep) != 1L || is.na(nrep) || nrep < 1L) {
    stop(
      "Zresidual_CV_survival_coxph: `nrep` must be a positive integer.",
      call. = FALSE
    )
  }
  
  nrep_eff <- if (isTRUE(randomized)) nrep else 1L
  
  get_frailty_group_var <- function(object) {
    tl <- attr(stats::terms(object), "term.labels")
    lab <- tl[grepl("^frailty\\s*\\(", tl)]
    
    if (length(lab) == 0L) {
      return(NULL)
    }
    
    if (length(lab) > 1L) {
      stop(
        "Zresidual_CV_survival_coxph: only one frailty() term is currently supported.",
        call. = FALSE
      )
    }
    
    m <- regexec("^frailty\\s*\\(([^,\\)]+)", lab)
    g <- regmatches(lab, m)[[1]]
    
    if (length(g) < 2L) {
      stop(
        "Zresidual_CV_survival_coxph: cannot parse frailty() group variable.",
        call. = FALSE
      )
    }
    
    trimws(g[2L])
  }
  
  frailty_group_in_train_ok <- function(foldlist, group) {
    group <- factor(group)
    n <- length(group)
    
    for (fid in seq_along(foldlist)) {
      test_idx <- foldlist[[fid]]
      
      if (any(test_idx < 1L | test_idx > n)) {
        return(FALSE)
      }
      
      train_idx <- setdiff(seq_len(n), test_idx)
      
      missing_group <- setdiff(
        unique(as.character(group[test_idx])),
        unique(as.character(group[train_idx]))
      )
      
      missing_group <- missing_group[!is.na(missing_group)]
      
      if (length(missing_group) > 0L) {
        return(FALSE)
      }
    }
    
    TRUE
  }
  
  check_frailty_group_in_train <- function(foldlist, group) {
    group <- factor(group)
    n <- length(group)
    
    for (fid in seq_along(foldlist)) {
      test_idx <- foldlist[[fid]]
      
      if (any(test_idx < 1L | test_idx > n)) {
        stop(
          sprintf(
            "Zresidual_CV_survival_coxph: fold %s contains row indices outside 1:nrow(data_used).",
            fid
          ),
          call. = FALSE
        )
      }
      
      train_idx <- setdiff(seq_len(n), test_idx)
      
      missing_group <- setdiff(
        unique(as.character(group[test_idx])),
        unique(as.character(group[train_idx]))
      )
      
      missing_group <- missing_group[!is.na(missing_group)]
      
      if (length(missing_group) > 0L) {
        stop(
          sprintf(
            paste0(
              "Zresidual_CV_survival_coxph: fold %s has frailty group(s) ",
              "absent from its training data: %s."
            ),
            fid,
            paste(missing_group, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
    
    invisible(TRUE)
  }
  
  mf_all <- stats::model.frame(
    object$terms,
    data = data,
    na.action = stats::na.pass
  )
  
  keep <- stats::complete.cases(mf_all)
  
  if (!all(keep)) {
    warning(
      sprintf(
        "Zresidual_CV_survival_coxph: removed %d rows with missing values before cross-validation.",
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
      "Zresidual_CV_survival_coxph: the survival response must contain a status column.",
      call. = FALSE
    )
  }
  
  status <- as.integer(ymat[, ncol(ymat)])
  time <- if (ncol(ymat) == 2L) ymat[, 1L] else ymat[, 2L]
  
  covariates <- mf[, -1L, drop = FALSE]
  n <- nrow(data_used)
  
  if (n < 2L) {
    stop(
      "Zresidual_CV_survival_coxph: at least two complete observations are required.",
      call. = FALSE
    )
  }
  
  has_frailty <- {
    specials <- tryCatch(
      attr(object$terms, "specials"),
      error = function(e) NULL
    )
    
    !is.null(specials$frailty) && length(specials$frailty) > 0L
  }
  
  type_is_frailty <- !is.null(type) &&
    identical(tolower(type), "frailty")
  
  is_frailty_model <- has_frailty || type_is_frailty
  
  frailty_group_var <- NULL
  
  if (is_frailty_model) {
    frailty_group_var <- get_frailty_group_var(object)
    
    if (is.null(frailty_group_var)) {
      stop(
        paste0(
          "Zresidual_CV_survival_coxph: frailty model was detected or requested, ",
          "but no frailty() term was found in the model formula."
        ),
        call. = FALSE
      )
    }
    
    if (!frailty_group_var %in% names(data_used)) {
      stop(
        sprintf(
          "Zresidual_CV_survival_coxph: frailty group variable '%s' is not found in `data`.",
          frailty_group_var
        ),
        call. = FALSE
      )
    }
  }
  
  ## Important:
  ## Do NOT add frailty_group_var as a factor to make_fold().
  ## make_fold() also checks factor x censor/status tables.
  ## If frailty id is included there, fold generation can become impossible.
  ## Instead, generate folds using ordinary covariates, then separately check
  ## that each test frailty group appears in the corresponding training data.
  covariates_for_fold <- covariates
  
  if (is_frailty_model && !is.null(frailty_group_var)) {
    if (frailty_group_var %in% names(covariates_for_fold)) {
      covariates_for_fold[[frailty_group_var]] <- NULL
    }
  }
  
  if (is.null(foldlist)) {
    if (is.null(nfolds)) {
      nfolds <- min(10L, n)
    }
    
    nfolds <- as.integer(nfolds)
    
    if (length(nfolds) != 1L || is.na(nfolds) || nfolds < 2L) {
      stop(
        "Zresidual_CV_survival_coxph: `nfolds` must be an integer >= 2.",
        call. = FALSE
      )
    }
    
    if (nfolds > n) {
      stop(
        "Zresidual_CV_survival_coxph: `nfolds` cannot exceed the number of observations.",
        call. = FALSE
      )
    }
    
    if (is_frailty_model && !is.null(frailty_group_var)) {
      fold_ok <- FALSE
      fold_attempt <- 0L
      
      while (!fold_ok && fold_attempt < 500L) {
        fold_attempt <- fold_attempt + 1L
        
        foldlist_try <- make_fold(
          fix_var = covariates_for_fold,
          y = y,
          k = nfolds,
          censor = status
        )
        
        fold_ok <- frailty_group_in_train_ok(
          foldlist = foldlist_try,
          group = data_used[[frailty_group_var]]
        )
      }
      
      if (!fold_ok) {
        stop(
          paste0(
            "Zresidual_CV_survival_coxph: failed to generate folds where every ",
            "test frailty group is also present in the corresponding training data. ",
            "Try reducing `nfolds`, checking group sizes, or providing a valid `foldlist`."
          ),
          call. = FALSE
        )
      }
      
      foldlist <- foldlist_try
    } else {
      foldlist <- make_fold(
        fix_var = covariates_for_fold,
        y = y,
        k = nfolds,
        censor = status
      )
    }
  } else {
    foldlist <- lapply(foldlist, function(x) sort(unique(as.integer(x))))
    flat <- unlist(foldlist, use.names = FALSE)
    
    if (!identical(sort(flat), seq_len(n))) {
      stop(
        paste0(
          "Zresidual_CV_survival_coxph: `foldlist` must be a non-overlapping ",
          "partition of 1:nrow(data_used), where data_used is the complete-case data."
        ),
        call. = FALSE
      )
    }
  }
  
  if (is_frailty_model && !is.null(frailty_group_var)) {
    check_frailty_group_in_train(
      foldlist = foldlist,
      group = data_used[[frailty_group_var]]
    )
  }
  
  if (is.null(log_pointpred) && is_frailty_model) {
    log_pointpred <- get0(
      "log_pointpred_survival_coxph_frailty",
      mode = "function",
      inherits = TRUE
    )
    
    if (!is.function(log_pointpred)) {
      stop(
        "Zresidual_CV_survival_coxph: cannot find `log_pointpred_survival_coxph_frailty`.",
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
            "Zresidual_CV_survival_coxph: fold %s refitting failed: %s",
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
    
    fold_extra_args <- extra_args
    
    if (is_frailty_model && is.null(fold_extra_args$traindata)) {
      fold_extra_args$traindata <- train_data
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
      fold_extra_args
    )
    
    z_fold <- tryCatch(
      do.call(Zresidual, z_call),
      error = function(e) {
        warning(
          sprintf(
            "Zresidual_CV_survival_coxph: fold %s Zresidual calculation failed: %s",
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
          "Zresidual_CV_survival_coxph: fold %s returned %d columns, expected %d.",
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
          "Zresidual_CV_survival_coxph: fold %s returned %d rows, expected %d.",
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
        "Zresidual_CV_survival_coxph: failed fold(s): %s.",
        paste(unique(failed_folds), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (anyNA(z_full)) {
    stop(
      "Zresidual_CV_survival_coxph: combined CV residual matrix contains NA values.",
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
  attr(z_full, "frailty_group_var") <- frailty_group_var
  attr(z_full, "is_frailty_model") <- is_frailty_model
  
  attr(z_full, "zcov") <- list(
    type = type,
    family = "coxph",
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
      original_row_index = which(keep),
      frailty_group_var = frailty_group_var,
      is_frailty_model = is_frailty_model
    )
  )
  
  z_full
}

