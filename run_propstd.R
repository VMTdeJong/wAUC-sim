c_args <- commandArgs(trailingOnly = T) # These are the arguments given through LINUX 

if (length(c_args) == 0) {
  # These args are for testing only.
  sce <- 1
  I   <- 1
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
  
  setwd("/hpc/shared/julius_bs/VMTdeJong/PropStd")
  source("1_auc.R")
  source("2_weighted_auc.R")
}

print(paste0("Starting a run with scenario = ", sce, ", I = ", I, ", rep = ", rep, 
       ", family = ", family, "."))

set.seed(rep * 100 + sce)

source("utils.R")
source("simulate.R")
source("sce_pars.R")
source("ps_methods.R")
source("generate_data.R")

# Parameters for development 
pars_dev <- make_pars()

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
                   "7" = make_pars_val(pars_dev, 
                                       mean_c = -.4,
                                       sce = sce),
                   "8" = make_pars_val(pars_dev,
                                       mean_c = .4,
                                       sce = sce),
                   stop("invalid scenario number 'sce'"))

# Parameters for "reference" validation data set,
    # case-mix is similar to development data set
    # predictor-outcome associations are similar to the validation data set.
pars_ref <- make_pars_val(pars_val, 
                          n = 1e4,
                          prev_b = pars_dev$prev_b,
                          mean_c = pars_dev$mean_c,
                          sd_c   = pars_dev$sd_c) 

pars_dev <- make_pars(n = n_dev)

path <- "Run 2020 11 12"
file <- paste0("sce_", sce, "_family_", family, "_ndev_", pars_dev$n, "_nval_", pars_val$n, "_r_", rep ,".RData")

out <- simulate_multiple(pars_dev,
                  pars_val,
                  pars_ref,
                  I = I,  
                  path = path,
                  file = file,
                  family = family)




