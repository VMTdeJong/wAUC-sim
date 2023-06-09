c_args <- commandArgs(trailingOnly = T) # These are the arguments given through LINUX 

if (length(c_args) == 0) {
  # These args are for testing only.
  sce <- 9
  I   <- 10
  rep <- 0
  family <- "binomial"
  n_dev <- 2000
  
  library(wAUC)
} else {
  # These args are for the real sim.
  sce <- as.numeric(c_args[1]) - 1
  I   <- as.numeric(c_args[2])
  rep <- as.numeric(c_args[3])
  family <- as.character(c_args[4])
  n_dev <- as.numeric(c_args[5])
  
  setwd("/hpc/julius_bs/shared/VMTdeJong/PropStd/")
  source("1_auc.R")
  source("2_weighted_auc.R")
}

print(paste0("Starting a run with scenario = ", sce, ", I = ", I, ", rep = ", rep, 
             ", n dev = ", n_dev, ", family = ", family, "."))

set.seed(rep * 100 + sce)

source("utils.R")
source("simulate.R")
source("sce_pars.R")
source("ps_methods.R")
source("generate_data.R")

# Parameters for development 
pars_dev <- switch(as.character(sce),
                   "8" = make_pars(trunc_c = c(-1, 1)),
                   "9" = make_pars(beta_w = 0),
                   "10" = make_pars(beta_w = 0,
                                    dir_w  = "x ~ w"),
                   "11" = make_pars(beta_w = 0,
                                    dir_w  = "w ~ x"),
                   "12" = make_pars(beta_w = 1),
                   "13" = make_pars(beta_w = 1,
                                    dir_w  = "x ~ w"),
                   "14" = make_pars(beta_w = 1,
                                    dir_w  = "w ~ x"),
                   make_pars()
                   )

# Parameters for validation data sets.
pars_val <- switch(as.character(sce),
                   "0" = make_pars_val(pars_dev, # = 1 on hpc
                                       sce = sce),
                   "1" = make_pars_val(pars_dev, # = 2 on hpc, etc   
                                       sd_c = .8,
                                       sce = sce),
                   "2" = make_pars_val(pars_dev,    
                                       sd_c = 1.2,
                                       sce = sce),
                   "3" = make_pars_val(pars_dev,    
                                       prev_b = 0.1,
                                       sce = sce),
                   "4" = make_pars_val(pars_dev,    
                                       prev_b = 0.4,
                                       sce = sce),
                   "5" = make_pars_val(pars_dev, 
                                       mean_c = -.4,
                                       sce = sce),
                   "6" = make_pars_val(pars_dev,
                                       mean_c = .4,
                                       sce = sce),
                   "7" = make_pars_val(pars_dev,
                                       trunc_c = c(-1, 1),
                                       sce = sce),
                   "8" = make_pars_val(pars_dev,
                                        trunc_c = c(-Inf, Inf), # i.e. no truncation, dev is truncated instead.
                                        sce = sce),
                   "9" = make_pars_val(pars_dev, 
                                       sd_c = 1.2,
                                       sce = sce),
                   "10" = make_pars_val(pars_dev, 
                                        sd_c = 1.2,
                                        sce = sce),
                   "11" = make_pars_val(pars_dev, 
                                        sd_c = 1.2,
                                        sce = sce),
                   "12" = make_pars_val(pars_dev, 
                                        sd_c = 1.2,
                                        sce = sce),
                   "13" = make_pars_val(pars_dev, 
                                        sd_c = 1.2,
                                        sce = sce),
                   "14" = make_pars_val(pars_dev, 
                                        sd_c = 1.2,
                                        sce = sce),
                   stop("invalid scenario number 'sce'"))

# Parameters for "reference" validation data set,
    # case-mix is similar to development data set
pars_ref <- make_pars_val(pars_dev, 
                          n = 1e4) 

# The sample size needs to be the same for each scenario, so better to set it here:
pars_dev$n <- n_dev
pars_val$n <- 2000

if (as.integer(sce) <= 8) {
  ps_methods <- list("ps_ignore",
                    "ps_odds_lin",
                    "ps_odds_lin_spl")
} else {
  ps_methods <- list("ps_ignore",
                    "ps_odds_lin_spl",
                    "ps_odds_lin_spl_w"
  )
}

path <- "Run_2022_11_15"
file <- paste0("sce_", sce, "_family_", family, "_ndev_", pars_dev$n, "_nval_", pars_val$n, "_r_", rep ,".RData")

print(system.time(out <- simulate_multiple(
  pars_dev,
  pars_val,
  pars_ref,
  ps_methods = ps_methods,
  I = I,  
  path = path,
  file = file,
  family = family)))




