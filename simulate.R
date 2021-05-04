#' Simulate one repetition of one scenario
#' 
#' @param pars_dev lists of scenario parameters. See sce_pars.R
#' @param pars_val ""                         
#' @param pars_ref "" 
#' @param ps_methods List of character names of propensity score memthods. See ps_methods.R
#' @param val_measures List of character names of validation measures. See val_measures.R. 
#'                     Must accept predicted value p, outcome y, propensity weight w and newdata.
#' @param save_extra Logical to indicate whether unncessary objects should be saved (for testing only).
#' @param I Number of simulation iterations
#' @param family family for prediction model AND propensity model
#' @param file File for saving

simulate_multiple <- function(pars_dev, 
                              pars_val, 
                              pars_ref,
                              ps_methods = list("ps_ignore",
                                                "ps_odds_lin",
                                                "ps_odds_lin_spl"
                                                ), 
                              val_measures = list("slope",
                                                  "citl",
                                                  "wAUC"), 
                              save_extra = F, 
                              I = 1, 
                              family = "binomial",
                              path = "test",
                              file = "multiple_runs.RData",
                              ...) {
  seed_file  <- paste0(path, "/seed_" , file)
  res_file   <- paste0(path, "/res_"  , file)
  extra_file <- paste0(path, "/extra_" , file)
  
  error_handler <- function(e) {
    print(e)
    # return(data.frame(est = rep(NA, length(ps_methods) * length(val_measures)), 
    #                   se = NA, 
    #                   'ci.2.5 %' = NA, 
    #                   'ci.97.5 %' = NA,
    #                   measure = rep(unlist(val_measures), length(ps_methods)), 
    #                   method = rep(unlist(ps_methods), length(val_measures)), 
    #                   i = i,
    #                   sce = pars_val$sce, 
    #                   family = family)
    # )
    return(NaN)
  }

  res <- list()
  seeds <- list()
  for (i in seq_len(I)) {
    if (!i %% 10) print(paste0("Starting i = ", i))
    seeds[[i]] <- .Random.seed
    save(seeds, file = seed_file)
    
    sim <- tryCatch(simulate_once(pars_dev = pars_dev,
                         pars_val = pars_val,
                         pars_ref = pars_ref,
                         ps_methods = ps_methods,
                         val_measures = val_measures, 
                         family = family,
                         save_extra = save_extra || i == 1,
                         ...),
                    error = error_handler)
    if (is.list(sim)) {
      res[[i]] <- cbind(sim$ests, 
                        i = i, 
                        sce = pars_val$sce, 
                        family = family,
                        n_dev = pars_dev$n,
                        n_val = pars_val$n,
                        n_ref = pars_ref$n)
      save(res, file = res_file)
      
      if (save_extra || i == 1)
        save(sim, file = extra_file)
    }
  }
  res <- Reduce(rbind, res)
  save(res, file = res_file)
  
  res
}

simulate_once <- function(pars_dev, 
                          pars_val, 
                          pars_ref, 
                          ps_methods, 
                          val_measures, 
                          family = binomial,
                          save_extra = F,
                          eps = 1e-8,
                          ...) {
  out <- list()
  out$call <- match.call()
  out$seed <- .Random.seed

  # Sample data
  data_dev <- generate_data(pars_dev, set = 1)  # Development
  data_val <- generate_data(pars_val, set = 0)  # Validation
  data_ref <- generate_data(pars_ref, set = NA) # For reference values for validation.
  
  # Develop model
  # pm = prediction model
  pm_fit <- glm(make_outcome_formula(data_dev), family = family, data = data_dev) 
  
  # Make predictions
  pm_p <- predict(pm_fit, newdata = data_val, type = "response")
  
  # To prevent extreme predictions when the gaussian fam is used.
  if (identical(match.fun(family), gaussian)) {
    pm_p[pm_p < eps] <- eps
    pm_p[pm_p > (1 - eps)] <- 1 - eps
  }
  
  ## Propensity score validation
  ## Loop 1: ps_methods
  # Propensity score = ps
  ests <- ps_values <- list()
  
  for (psm in seq_along(ps_methods)) {
    ps_method <- match.fun(ps_methods[[psm]])
    ps_fit <- ps_method(data_dev = data_dev, data_val = data_val, family = family, ...)  
    ps_w <- predict(ps_fit, newdata = data_val)  

    if (identical(match.fun(family), gaussian)) {
      ps_w[ps_w < eps] <- eps
      ps_w[ps_w > (1 - eps)] <- 1 - eps
    }
    
    if (save_extra) 
      ps_values[[ps_methods[[psm]]]] <- ps_w
    
    # Loop 2: val_measures
    val_ests <- list()
    for (vm in seq_along(val_measures)) {
      val_measure <- match.fun(val_measures[[vm]])
      val_fit <- val_measure(p = pm_p, lp = logit(pm_p), y = data_val$y, 
                             w = ps_w, newdata = data_val, family = family, ...)    
      val_est <- get_estimates(val_fit)                                                 
      val_ests[[val_measures[[vm]]]] <- val_est
    }
    ests[[ps_methods[[psm]]]] <- as.my.data.frame(val_ests) # simplify data, pt 1
  }
  
  # simplify data, pt 2
  for (i in seq_along(ests))
    ests[[i]]$method <- names(ests)[i]
  out$ests <- Reduce(rbind, ests)
  
  ## Validation in a reference data set
  pm_p_ref <- predict(pm_fit, newdata = data_ref, type = "response")
  
  if (identical(match.fun(family), gaussian)) {
    pm_p_ref[pm_p_ref < eps] <- eps
    pm_p_ref[pm_p_ref > (1 - eps)] <- 1 - eps
  }
  
  ref_ests <- list()
  for (vm in seq_along(val_measures)) {
    ref_measure <- match.fun(val_measures[[vm]])
    ref_fit <- ref_measure(p = pm_p_ref, lp = logit(pm_p_ref), y = data_ref$y, 
                           w = rep(1, pars_ref$n), newdata = data_ref, family = family, ...)    
    ref_est <- get_estimates(ref_fit)                                                 
    ref_ests[[val_measures[[vm]]]] <- ref_est
  }
  
  # simplify data, pt 3
  ref_ests_df <- as.my.data.frame(ref_ests)
  ref_ests_df$method <- "ref"
  out$ests <- rbind(out$ests, ref_ests_df)
  
  ###############
  
  ## 'Validation' in the development data set
  pm_p_dev <- predict(pm_fit, newdata = data_dev, type = "response")
  
  if (identical(match.fun(family), gaussian)) {
    pm_p_dev[pm_p_dev < eps] <- eps
    pm_p_dev[pm_p_dev > (1 - eps)] <- 1 - eps
  }
  
  dev_ests <- list()
  for (vm in seq_along(val_measures)) {
    dev_measure <- match.fun(val_measures[[vm]])
    dev_fit <- dev_measure(p = pm_p_dev, lp = logit(pm_p_dev), y = data_dev$y, 
                           w = rep(1, pars_dev$n), newdata = data_dev, family = family, ...)    
    dev_est <- get_estimates(dev_fit)                                                 
    dev_ests[[val_measures[[vm]]]] <- dev_est
  }
  
  # simplify data, pt 4
  dev_ests_df <- as.my.data.frame(dev_ests)
  dev_ests_df$method <- "dev"
  out$ests <- rbind(out$ests, dev_ests_df)
  
  ###########
  
  if (save_extra) {
    out$data_dev <- data_dev
    out$data_val <- data_val
    out$data_ref <- data_ref
    out$pm_fit <- pm_fit
    out$ps_fit <- ps_fit
    out$val_fit <- val_fit
    out$ref_fit <- ref_fit
    out$ps_values <- ps_values
    out$pm_p <- pm_p
  }
    
  out
}

