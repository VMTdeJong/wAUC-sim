#' Generate data
#' 
#' @param sce_pars see sce_pars.R
#' @param set Not used; returned for reference only.
#' @param ... Not used.

#' Code book for data generation
#'
#' @details y outcome
#'          x predictors
#'          set development or validation set. 

generate_data <- function(sce_pars, set, ...) {
  s <- sce_pars
  n_b <- length(s$beta_b)
  n_c <- length(s$beta_c)
  n_p <- n_b + n_c
  n  <- s$n
  
  x_b <- matrix(rbinom(n_b * n, 1, s$prev_b), nrow = n, ncol = n_b, byrow = T)
  x_c <- matrix(rnorm(n_c * n, mean = s$mean_c, sd = s$sd_c), nrow = n, ncol = n_c, byrow = T)
  x <- cbind(x_c, x_b)

  colnames(x) <- paste("x", seq_len(n_p), sep = "")
  
  lp <- s$beta_0 + x_b %*% s$beta_b + x_c %*% s$beta_c
  p  <- inv.logit(lp)  
  y  <- rbinom(n, 1, p)
  
  as.data.frame(cbind(y, x, set))
}

#' @example 
# d <- generate_data(make_pars(), 1)
# head(d)
# str(d)


