#' Code book for scenario parameters
#' 
#' @param set setting: 1 or 0, 1 = development, 0 = validation
#' @param n sample size per study
#' @param beta: coefficient: _0 = intercept, _b = binary variable(s), _c = continuous variable(s)
#' Note: _c is beta 1, _b is beta 2, _w = continuous variable that is not part of prediction model.
#' @param dir_w Direction of causal effect of w and x.
#' @param prev_b prevalence(s) of binary variable(s)
#' @param mean_c mean(s) of continuous variable(s)
#' @param sd_c sd(s) of continuous variable(s)
#' @param trunc_c At what values should c be truncated? c(-Inf, Inf) for no truncation

make_pars <- function(n       = 2000,
                      set     = 1,
                      beta_0  = -2,
                      beta_c  = c(1),
                      beta_b  = c(1),
                      beta_w  = NULL,
                      dir_w   = "x, w",
                      prev_b  = c(.2),
                      mean_c  = c(0),
                      sd_c    = c(1),
                      trunc_c = c(-Inf, Inf),
                      sce     = 0
                      ) {
  list(n       = n,
       set     = set,
       beta_0  = beta_0,
       beta_b  = beta_b,
       beta_c  = beta_c,
       beta_w  = beta_w,
       dir_w   = dir_w,
       prev_b  = prev_b,
       mean_c  = mean_c,
       sd_c    = sd_c,
       trunc_c = trunc_c,
       sce     = sce
       )
}


make_pars_val <- function(pars, n = pars$n, set = 0,...) {
  dots <- list(...)

  for (item in names(dots))
    if (is.null((pars[[item]] <- dots[[item]])))
      stop("All arguments in ... must be named.")

  pars[["set"]] <- set
  pars[["n"]] <- n
  pars
}

make_outcome_formula <- function(data) {
  rhs <- grep("x", colnames(data), ignore.case = TRUE)
  formula(paste("y ~ ", paste(colnames(data)[rhs], collapse = " + "), sep = ""))
}


#' @example 
#' sce_pars <- make_pars()
#' val_sce_pars <- make_pars_val(sce_pars)
