logit <- function(x)
  log(x/(1-x))

inv.logit <- function(x)
  1/(1+exp(-x))

as.my.data.frame <- function(l) {
  df <- as.data.frame(t(data.frame(l )))
  df$measure <- row.names(df)
  df
}

mgrep <- function(patterns, x, ignore.case = TRUE) 
  unlist(mapply(grep, patterns, MoreArgs = list(x = x, ignore.case = ignore.case)))

compare_avg_to_ref <- function(x) {
  
  out <- list()
  
  for (f in sort(unique(x$family))) {
    fa <- x[x$family == f, ]
    
    for (s in sort(unique(fa$sce))) {
      sc <- fa[fa$sce == s, ] 
      
      for (m in sort(unique(sc$measure))) {
        me <- sc[sc$measure == m, ] 
        
        me$mean <- me$mean - me[me$method == "ref", "mean"][[1]]
        me$SD   <- sqrt(me$SD^2 + me[me$method == "ref", "SD"][[1]]^2)
        me$SE   <- sqrt(me$SE^2 + me[me$method == "ref", "SE"][[1]]^2)
        me$estd_SE <- NULL
        me[me$method == "ref", c("SD", "SE")] <- 0
        
        out[[length(out) + 1]] <- me
      }
      
    }
  }
  Reduce(rbind, out)
}

compare_to_ref <- function(object) {
  object <- as.data.frame(object)
  
  object$`ci.2.5 %`  <- NULL
  object$`ci.97.5 %` <- NULL
  
  n_methods <- length(unique(object$method))
  
  out <- list()
  
  for (f in sort(unique(object$family))) {
    fa <- object[object$family == f, ]
    
    for (s in sort(unique(fa$sce))) {
      sc <- fa[fa$sce == s, ] 
      
      for (m in sort(unique(sc$measure))) {
        me <- sc[sc$measure == m, ] 
        me <<- me
        me$error    <- me$est - rep(me[me$method == "ref", "est"], each = n_methods)
        me$sq.error <- me$error^2
        me$se       <- sqrt(me$se^2 + rep((me[me$method == "ref", "se"])^2, each = n_methods))
        me[me$method == "ref", "se"] <- 0
        
        out[[length(out) + 1]] <- me
      }
    }
  }
  Reduce(rbind, out)
}

summarise_results <- function(x)
  x %>%
  group_by(sce, measure, method, family, n_dev) %>%
  summarise(mean = mean(est),
            SD = sd(est),
            mean_SE = SD/sqrt(n()),
            
            median = median(est),
            median_lb = bootstrap(est, median, quantile, .025),
            median_ub = bootstrap(est, median, quantile, .975),
            
            bias = mean(error),
            bias_lb = bootstrap(error, mean, quantile, .025),
            bias_ub = bootstrap(error, mean, quantile, .975),
            
            median_bias = median(error),
            median_bias_lb = bootstrap(error, median, quantile, .025),
            median_bias_ub = bootstrap(error, median, quantile, .975),
            
            RMSE = sqrt(mean(sq.error)),
            RMSE_lb = bootstrap(sq.error, function(z) sqrt(mean(z)), quantile, .025),
            RMSE_ub = bootstrap(sq.error, function(z) sqrt(mean(z)), quantile, .975),
            
            median_RMSE = sqrt(median(sq.error)),
            median_RMSE_lb = bootstrap(sq.error, function(z) sqrt(median(z)), quantile, .025),
            median_RMSE_ub = bootstrap(sq.error, function(z) sqrt(median(z)), quantile, .975),
            
            estd_SE = mean(se),
            n = n())

plot_results <- function(x, 
                         measure, 
                         stat = c("RMSE", "bias", "mean", 
                                  "median RMSE", "median bias", "median")[1], 
                         family = "binomial",
                         level = .95,
                         ref = FALSE,
                         sce_levels = unique(x$sce),
                         ...) {
  library(ggplot2)
  library(gridExtra)
  library(patchwork)
  
  x <- as.data.frame(x)
  
  if (!ref) {x <- x[!x$method == "ref", ]}
  
  df <- x[x$measure == measure, ]
  df <- df[df$family  == family, ]
  df$sce <- as.factor(df$sce)
  
  colnames(df)[colnames(df) == "sce"]    <- "Scenario"
  colnames(df)[colnames(df) == "method"] <- "Method"
  colnames(df)[colnames(df) == "bias"]   <- "Bias"
  colnames(df)[colnames(df) == "mean"]   <- "Mean"
  
  df$Scenario <- factor(df$Scenario, levels = sce_levels)
  
  if (stat == "mean" || stat == "Mean") {
    stat <- "Mean"
    df$lower <- df$mean_lb
    df$upper <- df$mean_ub
  }
  
  if (stat == "bias" || stat == "Bias") {
    stat <- "Bias"
    df$lower <- df$bias_lb
    df$upper <- df$bias_ub
  }
  
  if (stat == "rmse" || stat == "RMSE") {
    stat <- "RMSE"
    df$lower <- df$RMSE_lb
    df$upper <- df$RMSE_ub
  }
  
  if (stat == "median" || stat == "Median") {
    stat <- "Median"
    df$lower <- df$median_lb
    df$upper <- df$median_ub
  }
  
  if(stat == "median_bias"|| stat == "median bias") {
    stat <- "median_bias"
    df$lower <- df$median_bias_lb
    df$upper <- df$median_bias_ub
  }
  
  if (stat == "median RMSE" || stat == "median rmse" | stat == "median_RMSE" || stat == "median_rmse") {
    stat <- "median_RMSE"
    df$lower <- df$median_RMSE_lb
    df$upper <- df$median_RMSE_ub
  }
  
  if (measure == "citl")
    mlab <- "Calibration in the large"
  
  if (measure == "slope")
    mlab <- "Calibration slope"
  
  if (measure == "wAUC")
    mlab <- "Concordance statistic"
  
  df$Method[df$Method == "dev"] <- "Development data"
  df$Method[df$Method == "ps_ignore"] <- "Unweighted / naive"
  df$Method[df$Method == "ps_odds_lin"] <- "Propensity"
  df$Method[df$Method == "ps_odds_lin_w"] <- "Propensity, with extra covariate"
  df$Method[df$Method == "ps_odds_lin_spl"] <- "Propensity, splines"
  df$Method[df$Method == "ps_odds_lin_spl_w"] <- "Propensity, splines, with extra covariate"

  possible_levels <- c(
    "Propensity",
    "Propensity, splines",
    "Propensity, with extra covariate",
    "Propensity, splines, with extra covariate",
    "Unweighted / naive",
    "Development data")

  df$Method <- factor(df$Method, levels = possible_levels)
  legend_rows <- ceiling(length(unique(df$Method))/3)
  
  # See for fixing color scheme: 
    # https://stackoverflow.com/questions/19068432/ggplot2-how-to-use-same-colors-in-different-plots-for-same-factor
  # I changed it to color blind friendly colors:
  # https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000

  selection <- possible_levels %in% unique(df$Method)
  cols <- c("#648FFF", "#DC267F", "#FE6100", "#785EF0", "#FFB000", "gray31")[selection]
  names(cols) <- possible_levels[selection]
  
  # And keep shapes constant across plots
  shapes <- 1:6
  names(shapes) <- possible_levels
  shapes <- shapes[selection]

  dodge <- 1/2
  
  labels <- c("500" = "n development = 500", "2000" = "n development = 2000")
  
  ggplot(data = df,
         aes_string(x = "Scenario", y = stat, group = "Method",  ...)) +
    geom_point(aes(color = Method, shape = Method),
               position = position_dodge(dodge)) +
    geom_errorbar(width = 1,
                  aes(ymin = lower, ymax = upper, color = Method),
                  position = position_dodge(dodge)) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = shapes) +
    facet_grid( . ~ n_dev, labeller=labeller(n_dev = labels)) +
    guides(colour = guide_legend(nrow = 2))
}

#' @param x variable to draw bootstrap samples from
#' @param fun function to apply on each bootstrap sample. E.g. mean, median, rmse
#' @param bootstrap_fun function to estimate statistic to return. E.g. SD for SE, or quantile for percentile
#' @param ... args passed to bootstrap_fun
#' @param n number of bootstrap samples
#' @param seed seed passed to set.seed, so that multiple calls have the same seed (dplyr::summarise will
#' then have the same samples for each method). Set to NULL for random seed.

bootstrap <- function(x, fun, bootstrap_fun = sd, ..., n = 5000, seed = 12345) {
  fun <- match.fun(fun)
  bootstrap_fun <- match.fun(bootstrap_fun)
  out <- rep(NA, n)
  
  if (!is.null(seed)) {set.seed(seed)}
  
  for (i in seq_len(n)) {
    z <- sample(x, size = length(x), replace = TRUE)
    out[i] <- fun(z)
  }
  bootstrap_fun(out, ...)
}

#' Get estimates from the (weighted) concordance / AUC
#' 
#' @param lp linear predictor
#' @param newdata validation data
#' @param w weights
#' @param outcome_name name of the outcome, as generated in generate_data

# Obtain estimates from (weighted) validation model/fit
get_estimates <- function(object, ...)
  UseMethod("get_estimates")

get_estimates.wAUC <- function(object, ...)
  c(est = object$estimate,
    se  = object$statistics$se,
    ci  = c("2.5 %" = object$ci.lb, "97.5 %" = object$ci.ub))
