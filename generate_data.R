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

generate_data_simple <- function(sce_pars, set, ...) {
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

# Function used in first simulation
# generate_data <- function(sce_pars, set, ...) {
#   if (sce_pars$trunc_c[1] >= sce_pars$trunc_c[2])
#     stop("Lower limit of trunc_c must be smaller than upper limit, or sampling would run forever.")
#   
#   d1 <- generate_data_simple(sce_pars, set, ...)
#   d1 <- d1[d1$x1 >= sce_pars$trunc_c[1] & d1$x1 <= sce_pars$trunc_c[2], ]
# 
#   while (nrow(d1) < sce_pars$n) {
#     d2 <- generate_data_simple(sce_pars, set, ...)
#     d2 <- d2[d2$x1 >= sce_pars$trunc_c[1] & d2$x1 <= sce_pars$trunc_c[2], ]
#     d1 <- rbind(d1, d2)
#   }
#   
#   d1[seq_len(sce_pars$n), ]
# }

generate_data_w <- function(sce_pars, set, ...) {
  s <- sce_pars
  n_b <- length(s$beta_b)
  n_c <- length(s$beta_c)
  n_w <- length(s$beta_w)
  n  <- s$n
  
  if (sce_pars$dir_w == "x, w"){
    x_b <- matrix(rbinom(n_b * n, 1, s$prev_b),                 nrow = n, ncol = n_b, byrow = T)
    x_c <- matrix(rnorm(n_c * n, mean = s$mean_c, sd = s$sd_c), nrow = n, ncol = n_c, byrow = T)
    w   <- matrix(rnorm(n_w * n, mean = 0, sd = 1),             nrow = n, ncol = n_w, byrow = T)
  } else if (sce_pars$dir_w == "x ~ w") {
    x_b <- matrix(rbinom(n_b * n, 1, s$prev_b),                     nrow = n, ncol = n_b, byrow = T)
    w   <- matrix(    rnorm(n_w * n, mean = 0, sd = 1),             nrow = n, ncol = n_w, byrow = T)
    x_c <- matrix(w + rnorm(n_c * n, mean = s$mean_c, sd = s$sd_c), nrow = n, ncol = n_c, byrow = T)
  } else if (sce_pars$dir_w == "w ~ x") {
    x_b <- matrix(rbinom(n_b * n, 1, s$prev_b),                       nrow = n, ncol = n_b, byrow = T)
    x_c <- matrix(      rnorm(n_c * n, mean = s$mean_c, sd = s$sd_c), nrow = n, ncol = n_c, byrow = T)
    w   <- matrix(x_c + rnorm(n_w * n, mean = 0, sd = 1),             nrow = n, ncol = n_w, byrow = T)
  }
  
  lp <- s$beta_0 + x_b %*% s$beta_b + x_c %*% s$beta_c + w %*% s$beta_w
  p  <- inv.logit(lp)  
  y  <- rbinom(n, 1, p)
  
  x <- cbind(x_c, x_b)
  colnames(x) <- paste("x", seq_len(n_b + n_c), sep = "")
  colnames(w) <- paste("w", seq_len(n_w)      , sep = "")
  
  as.data.frame(cbind(y, x, w, set))
}

generate_data <- function(sce_pars, set, ...) {
  if (sce_pars$trunc_c[1] >= sce_pars$trunc_c[2])
    stop("Lower limit of trunc_c must be smaller than upper limit, or sampling would run forever.")
  
  if (is.null(sce_pars$beta_w)) {
    data_generation_fun <- generate_data_simple 
  } else {
    data_generation_fun <- generate_data_w
    }
  
  d1 <- data_generation_fun(sce_pars, set, ...)
  d1 <- d1[d1$x1 >= sce_pars$trunc_c[1] & d1$x1 <= sce_pars$trunc_c[2], ]
  
  while (nrow(d1) < sce_pars$n) {
    d2 <- data_generation_fun(sce_pars, set, ...)
    d2 <- d2[d2$x1 >= sce_pars$trunc_c[1] & d2$x1 <= sce_pars$trunc_c[2], ]
    
    d1 <- rbind(d1, d2)
  }
  
  d1[seq_len(sce_pars$n), ]
}

#' @example 
# d <- generate_data(make_pars(), 1)
# head(d)
# str(d)


