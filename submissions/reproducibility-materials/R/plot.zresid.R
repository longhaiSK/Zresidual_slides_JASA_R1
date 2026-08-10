#' Scatterplot diagnostics for Z-residuals
#'
#' Produces scatterplots of Z-residuals against observation index, linear
#' predictors, model covariates, or a user-supplied x-axis variable.
#'
#' Metadata supplied through `zcov` is typically returned by `Zcov()` and is
#' used to provide covariates, linear predictors, response type labels, and
#' plotting groups for hurdle, zero/logistic, count, and survival models.
#'
#' @param x A numeric vector or matrix of Z-residuals, typically returned by
#'   `Zresidual()`. Columns represent residual replications.
#' @param zcov Optional metadata returned by `Zcov()`.
#' @param irep Integer vector specifying which column(s) of `x` to plot.
#' @param ylab Label for the y-axis.
#' @param normality.test Character vector specifying which diagnostic p-values
#'   to display. Supported values are `"SW"`, `"AOV"`, and `"BL"`.
#' @param k.test Integer controlling grouping used by the diagnostic tests.
#' @param x_axis_var Variable used on the x-axis. It may be `"index"`, `"lp"`,
#'   `"covariate"`, a covariate name in `zcov$covariates`, a length-n vector,
#'   or a function `function(z, zcov)` returning either a vector or a list with
#'   components `values` and optionally `label`.
#' @param main.title Main title of the plot. If omitted, a default title is
#'   constructed from the residual type metadata when available.
#' @param outlier.return Logical; if `TRUE`, mark observations with
#'   `|Z| > outlier.value` and non-finite residuals.
#' @param outlier.value Numeric threshold used to define outliers.
#' @param category Optional grouping variable of length n used to modify point
#'   appearance.
#' @param outlier.set A named list of graphical arguments passed to `symbols()`
#'   and `text()` when annotating outliers.
#' @param xlab Label for the x-axis. If `NULL`, an automatic label is used.
#' @param my.mar Numeric vector passed to `par(mar = ...)`.
#' @param add_lowess Logical; if `TRUE`, add a LOWESS smooth when the x-axis is
#'   numeric.
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @return Invisibly returns a list with component `outliers`, containing the
#'   indices of observations flagged as outliers for the last plotted replicate.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30
#'   x <- rnorm(n)
#'   t_event <- rexp(n, rate = exp(0.3 * x))
#'   t_cens  <- rexp(n, rate = 0.4)
#'   status  <- as.integer(t_event <= t_cens)
#'   time    <- pmin(t_event, t_cens)
#'   dat <- data.frame(time = time, status = status, x = x)
#'
#'   fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
#'   z <- Zresidual(fit, data = dat, nrep = 1, seed = 1)
#'   zcov <- Zcov(fit, data = dat)
#'
#'   plot(z, zcov = zcov, x_axis_var = "index")
#'   plot(z, zcov = zcov, x_axis_var = "lp")
#'   plot(z, zcov = zcov, x_axis_var = "x")
#' }
#'
#' @seealso `Zresidual`, `Zcov`
#'
#' @method plot zresid
#' @export
plot.zresid <- function(x,
                        zcov = NULL,
                        irep = 1,
                        ylab = "Z-Residual",
                        normality.test = c("SW", "AOV", "BL"),
                        k.test = 10,
                        x_axis_var = "index",
                        main.title = NULL,
                        outlier.return = TRUE,
                        outlier.value = 3.5,
                        category = NULL,
                        outlier.set = list(),
                        xlab = NULL,
                        my.mar = c(5, 4, 4, 6) + 0.1,
                        add_lowess = FALSE,
                        ...) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  .as_matrix <- function(z) {
    if (is.null(dim(z))) matrix(z, ncol = 1L) else as.matrix(z)
  }
  
  .first_nonempty <- function(...) {
    vals <- list(...)
    for (v in vals) {
      if (!is.null(v) && length(v) > 0L) return(v)
    }
    NULL
  }
  
  .as_data_frame_or_null <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x)) return(x)
    as.data.frame(x)
  }
  
  .get_zcov <- function(z, zcov) {
    if (!is.null(zcov)) return(zcov)
    attr(z, "zcov") %||% attr(z, "Zcov")
  }
  
  .map_type_from_zcov <- function(zcov0) {
    if (is.null(zcov0) || !is.list(zcov0)) return(NULL)
    
    type0 <- zcov0$type
    if (!is.null(type0) && length(type0) >= 1L && nzchar(as.character(type0)[1L])) {
      return(as.character(type0)[1L])
    }
    
    kind <- zcov0$y_type_kind
    if (is.null(kind) || length(kind) == 0L) return(NULL)
    kind <- as.character(kind[1L])
    
    if (identical(kind, "censor")) return("survival")
    if (identical(kind, "trunc")) return("count")
    if (identical(kind, "hurdle")) return("hurdle")
    NULL
  }
  
  .fill_metadata <- function(z, zcov0) {
    if (is.null(zcov0) || !is.list(zcov0)) return(z)
    
    cov0 <- zcov0$covariates %||% zcov0$covariate
    lp0  <- zcov0$linear_pred %||% zcov0$linear.pred
    
    if (is.null(attr(z, "covariates")) && !is.null(cov0)) {
      attr(z, "covariates") <- cov0
    }
    
    if (is.null(attr(z, "linear.pred")) && !is.null(lp0)) {
      attr(z, "linear.pred") <- lp0
    }
    
    if (is.null(attr(z, "linear_pred")) && !is.null(lp0)) {
      attr(z, "linear_pred") <- lp0
    }
    
    if (is.null(attr(z, "type"))) {
      attr(z, "type") <- .map_type_from_zcov(zcov0)
    }
    
    if (is.null(attr(z, "zero_id"))) {
      zid <- NULL
      if (!is.null(zcov0$extra) && !is.null(zcov0$extra$zero_id)) {
        zid <- zcov0$extra$zero_id
      }
      if (is.null(zid)) zid <- zcov0$zero_id
      if (!is.null(zid)) attr(z, "zero_id") <- zid
    }
    
    if (is.null(attr(z, "y_type")) && !is.null(zcov0$y_type)) {
      attr(z, "y_type") <- zcov0$y_type
    }
    
    if (is.null(attr(z, "censored.status"))) {
      if (!is.null(zcov0$censored.status)) {
        attr(z, "censored.status") <- zcov0$censored.status
      } else if (identical(zcov0$y_type_kind, "censor")) {
        yt <- zcov0$y_type
        if (!is.null(yt)) attr(z, "censored.status") <- as.integer(yt == 0L)
      }
    }
    
    z
  }
  
  .resolve_xlab <- function(lbl) {
    if (is.null(lbl)) return(NULL)
    
    lbl_chr <- as.character(lbl)[1L]
    if (startsWith(lbl_chr, "tex(")) {
      if (requireNamespace("latex2exp", quietly = TRUE)) {
        latex_string <- sub("tex\\((.*)\\)", "\\1", lbl_chr)
        return(latex2exp::TeX(latex_string))
      }
      
      warning(
        "The 'latex2exp' package is not installed. LaTeX formatting for xlab will not be used.",
        call. = FALSE
      )
      return(sub("tex\\((.*)\\)", "\\1", lbl_chr))
    }
    
    lbl
  }
  
  .coerce_plot_x <- function(xv) {
    if (is.factor(xv)) return(xv)
    if (is.logical(xv)) return(as.integer(xv))
    xv
  }
  
  .resolve_x <- function(z, X, zcov0, x_axis_expr, xlab_from_function = NULL) {
    n <- NROW(z)
    keywords <- c("index", "lp", "covariate")
    
    if (is.character(X) && length(X) > 1L) X <- X[1L]
    
    if (is.function(X)) {
      xr <- X(z, zcov0)
      if (is.list(xr) && !is.null(xr$values)) {
        return(.resolve_x(
          z = z,
          X = xr$values,
          zcov0 = zcov0,
          x_axis_expr = x_axis_expr,
          xlab_from_function = xr$label
        ))
      }
      return(.resolve_x(
        z = z,
        X = xr,
        zcov0 = zcov0,
        x_axis_expr = x_axis_expr,
        xlab_from_function = xlab_from_function
      ))
    }
    
    is_user_vector <- !is.null(X) &&
      length(X) == n &&
      !(is.character(X) && length(X) == 1L && X %in% keywords)
    
    if (is_user_vector) {
      label0 <- xlab_from_function
      if (is.null(label0)) {
        lab <- paste(deparse(x_axis_expr), collapse = "")
        label0 <- if (identical(lab, "x_axis_var")) "X" else lab
      }
      
      return(list(
        values = .coerce_plot_x(X),
        label = label0,
        test_X = X,
        mode = "user_vector"
      ))
    }
    
    if (!is.character(X) || length(X) != 1L) {
      stop(
        "x_axis_var must be one of 'index', 'lp', 'covariate', a covariate name in zcov$covariates, a length-n vector, or a function(z, zcov).",
        call. = FALSE
      )
    }
    
    if (X == "index") {
      return(list(
        values = seq_len(n),
        label = "Index",
        test_X = "index",
        mode = "index"
      ))
    }
    
    if (X == "lp") {
      lp <- .first_nonempty(
        zcov0$linear_pred,
        zcov0$linear.pred,
        attr(z, "linear_pred"),
        attr(z, "linear.pred")
      )
      
      if (is.null(lp)) {
        stop("plot.zresid: x_axis_var = 'lp' requires `zcov = Zcov(...)` with `linear_pred`.", call. = FALSE)
      }
      
      lp <- as.numeric(lp)
      if (length(lp) != n) {
        stop(
          "plot.zresid: length(zcov$linear_pred) does not match nrow(Zresidual). Make sure Zresidual and Zcov use the same data rows.",
          call. = FALSE
        )
      }
      
      return(list(
        values = lp,
        label = "Linear Predictor",
        test_X = "lp",
        mode = "lp"
      ))
    }
    
    covs <- .first_nonempty(
      zcov0$covariates,
      zcov0$covariate,
      attr(z, "covariates")
    )
    
    covs <- .as_data_frame_or_null(covs)
    if (is.null(covs)) {
      stop(
        "plot.zresid: covariates are missing. Please pass `zcov = Zcov(fit, data = data)` or pass a length-n vector directly to `x_axis_var`.",
        call. = FALSE
      )
    }
    
    if (NROW(covs) != n) {
      stop(
        "plot.zresid: nrow(zcov$covariates) does not match nrow(Zresidual). Make sure Zresidual and Zcov use the same data rows.",
        call. = FALSE
      )
    }
    
    cn <- colnames(covs) %||% names(covs)
    if (is.null(cn) || length(cn) == 0L) {
      stop("plot.zresid: zcov$covariates has no column names.", call. = FALSE)
    }
    
    if (X == "covariate") {
      selected <- cn[1L]
      message(
        "x_axis_var = 'covariate' uses the first covariate: ", selected,
        ". Available covariates: ", paste(cn, collapse = ", ")
      )
    } else {
      selected <- X
      if (!(selected %in% cn)) {
        stop(
          "plot.zresid: covariate `", selected,
          "` was not found in zcov$covariates. Available covariates are: ",
          paste(cn, collapse = ", "),
          call. = FALSE
        )
      }
    }
    
    xv <- covs[[selected]]
    if (length(xv) != n) {
      stop("plot.zresid: selected covariate length does not match nrow(Zresidual).", call. = FALSE)
    }
    
    list(
      values = .coerce_plot_x(xv),
      label = selected,
      test_X = selected,
      mode = "covariate"
    )
  }
  
  .finite_abs_max <- function(z) {
    out <- suppressWarnings(max(abs(z[is.finite(z)]), na.rm = TRUE))
    if (!is.finite(out)) qnorm(0.9999) else out
  }
  
  .sanitize_z_for_plot <- function(zj, outlier.value) {
    original <- zj
    
    id_nan <- which(is.nan(zj))
    id_inf <- which(is.infinite(zj))
    
    value_notfinite <- NULL
    if (length(id_inf) > 0L) {
      value_notfinite <- rep(NA_character_, length(id_inf))
      value_notfinite[zj[id_inf] > 0] <- "Inf"
      value_notfinite[zj[id_inf] < 0] <- "-Inf"
      
      finite_abs_max <- .finite_abs_max(zj)
      zj[id_inf] <- sign(zj[id_inf]) * (finite_abs_max + 0.1)
      
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    
    if (length(id_nan) > 0L) {
      message("NaNs exist! The model or the fitting process has a problem!")
    }
    
    finite_abs_max <- .finite_abs_max(zj)
    ylim0 <- max(qnorm(0.9999), finite_abs_max)
    
    if (!is.finite(ylim0)) {
      stop(
        "plot.zresid: failed to compute finite ylim. Check whether Z-residuals are all NA/NaN/Inf.",
        call. = FALSE
      )
    }
    
    id_outlier <- which(abs(zj) > outlier.value | is.infinite(original))
    
    list(
      z = zj,
      id_nan = id_nan,
      id_inf = id_inf,
      value_notfinite = value_notfinite,
      id_outlier = id_outlier,
      ylim0 = ylim0
    )
  }
  
  .as_numeric_for_smooth <- function(v) {
    if (inherits(v, "POSIXt") || inherits(v, "Date")) return(as.numeric(v))
    if (is.factor(v)) return(as.numeric(v))
    v
  }
  
  .add_lowess_line <- function(xv, yv, args) {
    if (!isTRUE(add_lowess)) return(invisible(NULL))
    
    xv_num <- .as_numeric_for_smooth(xv)
    if (!is.numeric(xv_num)) return(invisible(NULL))
    
    log_opt <- args[["log"]]
    log_x <- !is.null(log_opt) && grepl("x", log_opt)
    
    ok <- is.finite(xv_num) & is.finite(yv)
    if (log_x) ok <- ok & (xv_num > 0)
    if (sum(ok) < 3L) return(invisible(NULL))
    
    if (log_x) {
      lw <- stats::lowess(x = log10(xv_num[ok]), y = yv[ok])
      graphics::lines(x = 10^lw$x, y = lw$y, col = "red", lwd = 3)
    } else {
      lw <- stats::lowess(x = xv_num[ok], y = yv[ok])
      graphics::lines(x = lw$x, y = lw$y, col = "red", lwd = 3)
    }
    
    invisible(NULL)
  }
  
  .draw_legends_outside_right <- function(legend.args, test.legend,
                                          mar_right_min = 6,
                                          gap_factor = 0.8) {
    if (is.null(legend.args) && is.null(test.legend)) return(invisible(NULL))
    
    mar <- graphics::par("mar")
    if (mar[4] < mar_right_min) {
      graphics::par(mar = c(mar[1L], mar[2L], mar[3L], mar_right_min))
    }
    
    old_xpd <- graphics::par("xpd")
    on.exit(graphics::par(xpd = old_xpd), add = TRUE)
    graphics::par(xpd = TRUE)
    
    usr <- graphics::par("usr")
    xlog <- graphics::par("xlog")
    
    if (xlog) {
      x_anchor <- 10^(usr[2L] + (usr[2L] - usr[1L]) * 0.01)
    } else {
      x_anchor <- usr[2L] + (usr[2L] - usr[1L]) * 0.01
    }
    y_top <- usr[4L]
    
    lg1 <- NULL
    if (!is.null(legend.args)) {
      la <- legend.args
      la$plot <- FALSE
      la$x <- x_anchor
      la$y <- y_top
      la$xjust <- 0
      la$yjust <- 1
      lg1 <- do.call(graphics::legend, la)
      
      la$plot <- TRUE
      do.call(graphics::legend, la)
    }
    
    if (!is.null(test.legend)) {
      tl <- test.legend
      tl$xpd <- TRUE
      tl$bg <- "white"
      tl$box.col <- NA
      
      if (!is.null(lg1)) {
        gap <- graphics::strheight("M", units = "user") * gap_factor
        x2 <- x_anchor
        y2 <- lg1$rect$top - lg1$rect$h - gap
      } else {
        x2 <- x_anchor
        y2 <- y_top
      }
      
      do.call(graphics::legend, c(list(x = x2, y = y2, xjust = 0, yjust = 1), tl))
    }
    
    invisible(NULL)
  }
  
  .make_groups <- function(z, zcov0, category, args) {
    n <- NROW(z)
    
    if (!is.null(category)) {
      if (length(category) != n) {
        stop("plot.zresid: length(category) must match nrow(Zresidual).", call. = FALSE)
      }
      
      unique_cats <- unique(category)
      pal <- if (!is.null(args[["col"]])) {
        args[["col"]]
      } else if (length(unique_cats) == 1L) {
        "red"
      } else {
        grDevices::rainbow(max(length(unique_cats), 2L))[seq_along(unique_cats)]
      }
      
      pchs <- if (!is.null(args[["pch"]])) {
        args[["pch"]]
      } else if (length(unique_cats) == 1L) {
        1
      } else {
        (1:25)[seq_along(unique_cats)]
      }
      
      return(list(
        col = pal[match(category, unique_cats)],
        pch = pchs[match(category, unique_cats)],
        legend = unique_cats,
        legend_col = pal,
        legend_pch = pchs,
        legend_title = deparse(substitute(category))
      ))
    }
    
    type <- attr(z, "type")
    if (is.null(type) || length(type) == 0L || !nzchar(as.character(type)[1L])) {
      type <- NULL
    } else {
      type <- as.character(type)[1L]
    }
    
    if (is.null(type)) {
      return(list(col = "black", pch = 1, legend = NULL, legend_col = NULL, legend_pch = NULL, legend_title = NULL))
    }
    
    if (identical(type, "survival")) {
      censored <- attr(z, "censored.status")
      if (is.null(censored) || length(censored) != n) {
        return(list(col = "black", pch = 1, legend = NULL, legend_col = NULL, legend_pch = NULL, legend_title = NULL))
      }
      
      censored <- as.integer(censored)
      idx <- censored + 1L
      idx[!(idx %in% c(1L, 2L))] <- NA_integer_
      
      col <- c("blue", "red")[idx]
      pch <- c(3, 2)[idx]
      col[is.na(col)] <- "black"
      pch[is.na(pch)] <- 1
      
      return(list(
        col = col,
        pch = pch,
        legend = c("Uncensored", "Censored"),
        legend_col = c("blue", "red"),
        legend_pch = c(3, 2),
        legend_title = NULL
      ))
    }
    
    if (identical(type, "hurdle")) {
      zero_id <- attr(z, "zero_id")
      if (is.null(zero_id)) zero_id <- integer(0L)
      
      is_zero <- seq_len(n) %in% zero_id
      return(list(
        col = c("blue", "red")[is_zero + 1L],
        pch = c(3, 2)[is_zero + 1L],
        legend = c("count", "zero"),
        legend_col = c("blue", "red"),
        legend_pch = c(3, 2),
        legend_title = NULL
      ))
    }
    
    if (type %in% c("zero", "logistic")) {
      yt <- .first_nonempty(zcov0$y_type, attr(z, "y_type"))
      if (is.null(yt) || length(yt) != n) {
        stop(
          "plot.zresid (zero/logistic): cannot find y_type. Pass `zcov = Zcov(...)` aligned with this Zresidual.",
          call. = FALSE
        )
      }
      
      fam_info <- NULL
      if (!is.null(zcov0$family)) fam_info <- as.character(zcov0$family)[1L]
      is_hurdle <- !is.null(fam_info) && grepl("^hurdle(_|$)", fam_info)
      
      if (is_hurdle) {
        labels <- c("count", "zero")
        codes <- c(1L, 0L)
      } else {
        labels <- c("0", "1")
        codes <- c(0L, 1L)
      }
      
      idx <- match(as.integer(yt), codes)
      if (anyNA(idx)) {
        stop("plot.zresid (zero/logistic): y_type contains values outside expected codes.", call. = FALSE)
      }
      
      return(list(
        col = c("blue", "red")[idx],
        pch = c(3, 2)[idx],
        legend = labels,
        legend_col = c("blue", "red"),
        legend_pch = c(3, 2),
        legend_title = NULL
      ))
    }
    
    if (identical(type, "count")) {
      return(list(
        col = "blue",
        pch = 3,
        legend = "count",
        legend_col = "blue",
        legend_pch = 3,
        legend_title = NULL
      ))
    }
    
    list(
      col = "red",
      pch = 2,
      legend = type,
      legend_col = "red",
      legend_pch = 2,
      legend_title = NULL
    )
  }
  
  .make_legend_args <- function(groups, args) {
    if (is.null(groups$legend) || length(groups$legend) == 0L) return(NULL)
    if (is.null(groups$legend_col) || is.null(groups$legend_pch)) return(NULL)
    
    default_legend <- list(
      legend = groups$legend,
      col = groups$legend_col,
      pch = groups$legend_pch,
      cex = 0.9,
      xpd = TRUE,
      bty = "n",
      title = groups$legend_title,
      horiz = FALSE,
      y.intersp = 1
    )
    
    user_legend_overrides <- args[!names(args) %in% c("col", "pch")]
    out <- modifyList(default_legend, user_legend_overrides)
    out[names(out) %in% names(formals(graphics::legend))]
  }
  
  .make_test_legend <- function(z, x_spec, zcov0, j, normality.test, k.test) {
    if (length(normality.test) == 0L) return(NULL)
    
    test_map <- c(SW = "sw.test.zresid", AOV = "aov.test.zresid", BL = "bartlett.test.zresid")
    normality.test <- toupper(normality.test)
    normality.test <- normality.test[normality.test %in% names(test_map)]
    if (length(normality.test) == 0L) return(NULL)
    
    pvals <- character(0L)
    
    for (tt in normality.test) {
      fun_name <- unname(test_map[tt])
      if (!exists(fun_name, mode = "function")) {
        pv_str <- "NA"
      } else {
        fun <- get(fun_name, mode = "function")
        pv <- tryCatch(
          {
            if (identical(tt, "SW")) {
              fun(z)
            } else {
              fun(z, X = x_spec, zcov = zcov0, k.anova = k.test, k.bl = k.test)
            }
          },
          error = function(e) rep(NA_real_, NCOL(z))
        )
        
        pv_j <- suppressWarnings(as.numeric(pv[j]))
        pv_str <- if (is.finite(pv_j)) sprintf("%3.2f", pv_j) else "NA"
      }
      
      pvals <- c(pvals, paste(tt, "-", pv_str))
    }
    
    list(
      legend = c(expression(bold("P-value:")), pvals),
      cex = 1,
      bty = "n",
      xpd = TRUE,
      adj = c(0, 0.5)
    )
  }
  
  .draw_outliers <- function(xv, yv, id_outlier, outlier.set) {
    if (!length(id_outlier)) return(invisible(NULL))
    
    default_outlier <- list(
      pos = 4,
      labels = id_outlier,
      cex = 0.8,
      col = "darkolivegreen4",
      add = TRUE,
      inches = FALSE,
      circles = rep((graphics::par("usr")[2L] - graphics::par("usr")[1L]) * 0.03, length(id_outlier)),
      fg = "darkolivegreen4",
      font = 2
    )
    
    outlier_args <- modifyList(default_outlier, outlier.set)
    
    text_fun <- getS3method("text", "default")
    text_args <- outlier_args[names(outlier_args) %in% names(formals(text_fun))]
    symbols_args <- outlier_args[names(outlier_args) %in% names(formals(graphics::symbols))]
    
    symbols_args$add <- NULL
    do.call(
      graphics::symbols,
      c(list(x = xv[id_outlier], y = yv[id_outlier], add = TRUE), symbols_args)
    )
    
    text_x <- xv[id_outlier]
    text_y <- yv[id_outlier]
    text_pos <- outlier_args$pos
    
    if (length(text_pos) == 1L) text_pos <- rep(text_pos, length(id_outlier))
    y_median <- stats::median(yv, na.rm = TRUE)
    
    for (ii in seq_along(id_outlier)) {
      if (!is.na(text_y[ii])) {
        text_pos[ii] <- if (text_y[ii] < y_median) 3 else 1
      }
    }
    
    text_args_no_pos <- text_args[names(text_args) != "pos"]
    do.call(
      graphics::text,
      c(list(x = text_x, y = text_y, pos = text_pos), text_args_no_pos)
    )
    
    invisible(NULL)
  }
  
  Zresidual <- .as_matrix(x)
  n_obs <- NROW(Zresidual)
  
  if (n_obs < 1L || NCOL(Zresidual) < 1L) {
    stop("plot.zresid: `x` must contain at least one observation and one residual column.", call. = FALSE)
  }
  
  irep <- as.integer(irep)
  if (anyNA(irep) || any(irep < 1L | irep > NCOL(Zresidual))) {
    stop("plot.zresid: `irep` contains invalid residual column index.", call. = FALSE)
  }
  
  zcov0 <- .get_zcov(Zresidual, zcov)
  if (is.null(zcov0)) zcov0 <- list()
  
  Zresidual <- .fill_metadata(Zresidual, zcov0)
  
  if (is.null(main.title)) {
    type <- attr(Zresidual, "type")
    main.title <- if (is.null(type) || length(type) == 0L || !nzchar(as.character(type)[1L])) {
      "Z-residual Scatterplot"
    } else {
      paste("Z-residual Scatterplot -", as.character(type)[1L])
    }
  }
  
  x_axis_expr <- substitute(x_axis_var)
  xinfo <- .resolve_x(
    z = Zresidual,
    X = x_axis_var,
    zcov0 = zcov0,
    x_axis_expr = x_axis_expr
  )
  
  xvec <- xinfo$values
  if (length(xvec) != n_obs) {
    stop("plot.zresid: resolved x-axis variable length does not match nrow(Zresidual).", call. = FALSE)
  }
  
  args <- list(...)
  
  groups <- .make_groups(Zresidual, zcov0, category, args)
  col <- if (!is.null(args[["col"]])) args[["col"]] else groups$col
  pch <- if (!is.null(args[["pch"]])) args[["pch"]] else groups$pch
  
  legend.args <- .make_legend_args(groups, args)
  
  op_mar <- graphics::par("mar")
  op_xpd <- graphics::par("xpd")
  on.exit(graphics::par(mar = op_mar, xpd = op_xpd), add = TRUE)
  
  last_outlier <- integer(0L)
  
  for (j in irep) {
    graphics::par(mar = my.mar)
    
    sanitized <- .sanitize_z_for_plot(Zresidual[, j], outlier.value = outlier.value)
    yvec <- sanitized$z
    Zresidual[, j] <- yvec
    last_outlier <- sanitized$id_outlier
    
    current_xlab <- .resolve_xlab(xlab %||% xinfo$label)
    
    default.plot <- modifyList(
      list(
        x = xvec,
        y = yvec,
        col = col,
        pch = pch,
        ylab = ylab,
        ylim = c(-sanitized$ylim0, sanitized$ylim0 + 1),
        main = main.title,
        xlab = current_xlab,
        font.lab = 2
      ),
      args
    )
    
    do.call(graphics::plot, default.plot)
    .add_lowess_line(xvec, yvec, args)
    
    test.legend <- NULL
    if (!identical(xinfo$mode, "index")) {
      test.legend <- .make_test_legend(
        z = Zresidual,
        x_spec = xinfo$test_X,
        zcov0 = zcov0,
        j = j,
        normality.test = normality.test,
        k.test = k.test
      )
    }
    
    .draw_legends_outside_right(legend.args, test.legend)
    
    if (isTRUE(outlier.return) && length(sanitized$id_outlier) > 0L) {
      .draw_outliers(xvec, yvec, sanitized$id_outlier, outlier.set)
    }
    
    if (length(sanitized$id_inf) > 0L) {
      graphics::text(
        x = xvec[sanitized$id_inf],
        y = yvec[sanitized$id_inf],
        labels = sanitized$value_notfinite,
        col = 2,
        pos = 4
      )
    }
    
    hlines <- c(1.96, 3)
    graphics::abline(h = c(hlines, -hlines), lty = 3, col = "grey")
    
    if (isTRUE(outlier.return) && length(sanitized$id_outlier) > 0L) {
      message("Outlier Indices : ", paste(sanitized$id_outlier, collapse = ", "))
    }
  }
  
  invisible(list(outliers = last_outlier))
}