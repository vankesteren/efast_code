# Data simulation for EFAST models
# Erik-Jan van Kesteren
# Last edited: 20190614

library(efast)
library(tidyverse)
library(pbapply)
library(parallel)


fit_models <- function(dat, M = 4, verbose = TRUE, ...) {
  list(
    efast = efast_hemi(data = dat, M = M, lh_idx = 1:17, rh_idx = 18:34,
                       verbose = verbose, optim.force.converged = TRUE, ...),
    efa = efast_efa(data = dat, M = M, verbose = verbose, 
                    optim.force.converged = TRUE, ...)
  )
}


# create 100 datasets in a tibble
dat_list <- lapply(1:100, function(i, ...) simulate_efast(...), 
                   N = 650L, lam_lat = .5*.85, lam_bil = .5, 
                   psi_cov = .5, cor_uniq = .12)

grid <- expand.grid(
  M = 2:10,
  dataset = 1:100
)

# Create a cluster
clus <- makeCluster(detectCores())
clusterExport(clus, c("fit_models", "dat_list", "grid"))
pkg <- clusterEvalQ(clus, {
  library(efast)
  library(tidyverse)
})

# Run simulation
result <- pblapply(1:nrow(grid), function(i) {
  args <- grid[i,]
  dat  <- dat_list[[args$dataset]]
  fit  <- fit_models(dat, M = args$M, verbose = FALSE, rotation = 'none', control = list(iter.max = 4000))
  return(list(efa_fit = fit$efa, efast_fit = fit$efast))
}, cl = clus)

# stop the cluster
stopCluster(clus)


tbl <- 
  as_tibble(grid) %>% 
  mutate(
    efa_fit   = map(result, function(o) o$efa_fit),
    efast_fit = map(result, function(o) o$efast_fit)
  )


# Save the tbl
write_rds(tbl, path = "simulations/output/overextraction.rds")
