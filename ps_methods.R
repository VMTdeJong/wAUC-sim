
############ Estimation methods ############ 
# Background functions.
ps_model <- function(data_dev, 
                     data_val, 
                     formula = make_ps_formula(data = data_dev), 
                     ps_estimation_method = "glm", 
                     family = binomial, ...) {
  data_comb <- rbind(data_dev, data_val)
  
  fitter <- match.fun(ps_estimation_method)
  fitter(formula = formula, data = data_comb, family = family, ...)
}

# data <- data.frame(y = 1, x1 = 3, x2 = 4, cluster = 1)
make_ps_formula <- function(data, splines = FALSE) {
  rhs <- grep("x", colnames(data), ignore.case = TRUE)
  
  predictors <- colnames(data)[rhs]
  if (any(splines)) {
    library(splines)
    predictors <- c(paste0("ns(", predictors[splines], ", df = 3)"), predictors[-splines])
  }
  out <- paste("set ~ ", paste(paste0(predictors, collapse = " + ")))

  formula(out)
}


# These have the same estimation functions, but have their own prediction functions below.
ps_p <- function(...) {
  out <- ps_model(...)
  class(out) <- c("ps_p", class(out))
  out
}

ps_odds <- function(...) {
  out <- ps_p(...)
  class(out) <- c("ps_odds", class(out))
  out
}

#### These methods build on the ones above. Currently the linear ones are redundant, they do not change the default setting.
ps_p_lin <- function(use.mean = FALSE, ...)
  ps_p(formula = make_ps_formula(data = list(...)$data_dev), ...)

ps_p_lin_spl <- function(use.mean = FALSE, ...)
  ps_p(formula = make_ps_formula(data = list(...)$data_dev, splines = 1L), ...)

ps_odds_lin <- function(...)
  ps_odds(formula = make_ps_formula(data = list(...)$data_dev), ...)

ps_odds_lin_spl <- function(...)
  ps_odds(formula = make_ps_formula(data = list(...)$data_dev, splines = 1L), ...)

### Method for ignoring propensity
ps_ignore <- function(...) {
  out <- list()
  class(out) <- "ps_ignore"
  out
}

############ Prediction methods ############ 
predict.ps_p <- function(object, newdata, ...)
  predict.glm(object, newdata = newdata, type = "response", ...)

predict.ps_odds <- function(object, newdata, ...) {
  p <- predict.ps_p(object, newdata = newdata, ...)
  p / (1 - p)
}

predict.ps_ignore <- function(object, newdata, ...)
  rep(1, nrow(newdata))
  
