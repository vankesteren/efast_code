# Data simulation for EFAST models
# Erik-Jan van Kesteren
# Last edited: 20200428

library(efast)
library(tidyverse)
library(pbapply)
library(parallel)


fit_models <- function(dat, M = 4, verbose = TRUE) {
  list(
    efast = efast_hemi(data = dat, M = M, lh_idx = 1:17, rh_idx = 18:34,
                       verbose = verbose, optim.force.converged = TRUE),
    efa = efast_efa(data = dat, M = M, verbose = verbose, 
                    optim.force.converged = TRUE)
  )
}

# This is how to run one iteration of the simulation
dat <- simulate_efast()
fit <- fit_models(dat)

# Create a condition grid
grid <- expand.grid(
  lambda  = c(0.5, 0.7),    # bil traits get this, lat traits this*0.85
  psi_cov = c(0, 0.5),      # covariance among latent variables
  res_cor = c(0, 0.2, 0.4)
)

# Create an output format 
result <- vector("list", length = nrow(grid))
nrep <- 120L
tbl <- 
  grid %>% 
  as_tibble() %>% 
  mutate(dataset = list(lapply(1:nrep, function(i) NULL))) %>% 
  mutate(efa_fit   = dataset,
         efast_fit = dataset)

# Create a cluster
clus <- makeCluster(detectCores() - 2)
clusterExport(clus, "fit_models")
pkg <- clusterEvalQ(clus, {
  library(efast)
  library(tidyverse)
})


# Run simulation
for (i in 1:nrow(grid)) {
  cat("\n\nCONDITION: ", i, "\n\n")
  arg <- grid[i,]
  clusterExport(clus, "arg")
  output <- pblapply(1:nrep, function(r) {
    dat <- simulate_efast(
      lam_lat  = arg$lambda*0.85, 
      lam_bil  = arg$lambda,
      psi_cov  = arg$psi_cov,
      cor_uniq = arg$res_cor
    )
    fit <- fit_models(dat, verbose = FALSE)
    list(dataset = dat, result = fit)
  }, cl = clus)
  tbl$dataset[[i]]   <- lapply(output, function(o) o$dataset)
  tbl$efa_fit[[i]]   <- lapply(output, function(o) o$result[["efa"]])
  tbl$efast_fit[[i]] <- lapply(output, function(o) o$result[["efast"]])
}

# Save the tbl
write_rds(tbl, path = "simulations/output/sim_run_efast_complete.rds")

# stop the cluster
stopCluster(clus)
