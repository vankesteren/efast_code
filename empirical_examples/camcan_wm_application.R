# CamCAN white matter tracts empirical analysis with EFAST
# Determining number of factors and performing likelihood ratio tests
# Then results of final model
# Erik-Jan van Kesteren
# Last edited: 20200428
library(tidyverse)
library(firatheme)
library(corrplot)
library(efast)
library(pbapply)
library(parallel)
library(xtable)
library(qgraph)
library(ggseg)

# load the data
ccwm <- read_rds("empirical_examples/data/camcan_wm.rds")


# plot the covmat
pdf("empirical_examples/plots/camcan_wm_cplot.pdf", 4, 4)
cplot(cor(ccwm), cl.pos = "r")
dev.off()

# extract area names
area_names <- 
  colnames(ccwm) %>% 
  str_extract("(?<=[LR]H\\_).*") %>% 
  unique()


#####################################
### Determining number of factors ###
#####################################

# create a grid of efast / efa parameters
results <- 
  expand.grid(method = c("EFA", "EFAST"), M = 2:11) %>% 
  as_tibble()

# create a model fit function
fit_model <- function(method, M) {
  if (method == "EFAST") {
    fit <- try(efast_hemi(data = ccwm, M = M, lh_idx = 1:21, rh_idx = 22:42, 
                          roi_names = area_names))
  } else {
    fit <- try(efast_efa(data = ccwm, M = M))
  }
}

# create a parallel cluster
clus <- makeCluster(12)
clusterExport(clus, c("ccwm", "area_names", "results", "fit_model"))
clusterEvalQ(clus, library(efast))

# fit the models
fits <- pblapply(1:nrow(results), function(i) {
  fit_model(results[i,]$method, results[i,]$M)
}, cl = clus)

stopCluster(clus)

# add them to the results
results$fit <- fits

# and save the results
write_rds(results, "empirical_examples/output/camcan_wm_fit.rds")

# optionally, load the checkpoint
results <- read_rds("empirical_examples/output/camcan_wm_fit.rds")

# extract fit measures
results_metrics <-
  results %>% 
  mutate(
    AIC    = sapply(fit, fitmeasures, fit.measures = "aic"),
    BIC    = sapply(fit, fitmeasures, fit.measures = "bic"),
    SSABIC = sapply(fit, fitmeasures, fit.measures = "bic2"),
    valid  = sapply(fit, lavaan:::lav_object_post_check)
  ) %>% 
  gather("Metric", "Value", -method, -M, -fit, -valid)

metrics_plot <- 
  results_metrics %>% 
  mutate(Convergence = case_when( valid ~ "Converged", 
                                  !valid ~ "Not converged")) %>% 
  ggplot(aes(x = M, y = Value, colour = method, alpha = Convergence, 
             size = Convergence)) + 
  geom_line(alpha = 0.3, size = 0.5) + 
  geom_point(stroke = 0) +
  theme_fira() +
  scale_colour_fira() + 
  scale_size_manual(values = c(2.5, 1.5)) +
  scale_alpha_manual(values = c(1, 0.3)) +
  labs(colour = "Model type", y = "", x = "Number of EFA factors", 
       main = "CamCAN white matter factor analysis metrics") +
  facet_wrap(~Metric) +
  theme(legend.title.align = 0) +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  guides(colour = guide_legend(order = 1), 
         alpha  = guide_legend(order = 2),
         size   = guide_legend(order = 2))

firaSave("empirical_examples/plots/camcan_wm_metrics.pdf", plot = metrics_plot, 
         width = 8, height = 4)

########################################
### Fitting and comparing best model ###
########################################

# load the checkpoint
results <- read_rds("empirical_examples/output/camcan_wm_fit.rds")

# best model is the EFAST model with M = 6
efast_fit <- results %>% filter(method == "EFAST", M == 6) %>% .$fit %>% .[[1]]
efa_fit   <- results %>% filter(method == "EFA",   M == 6) %>% .$fit %>% .[[1]]

# Network plots of latent covariance
pdf("empirical_examples/plots/camcan_wm_network_latents.pdf", 
    width = 12, height = 6)
par(mfrow = c(1,2))
Psi_efast <- crossprod(efast_fit@Model@H[[1]])
Psi_efa   <- crossprod(efa_fit@Model@H[[1]])
qgraph(Psi_efa - diag(6),   layout = "spring", posCol = firaCols[1], 
       title = 'A')
qgraph(Psi_efast - diag(6), layout = "spring", posCol = firaCols[1], 
       title = 'B')
dev.off()


# compare efast and efa using AIC, BIC, SSABIC, LR-test
diftest <- anova(efast_fit, efa_fit)
rownames(diftest) <- c("EFAST" , "EFA")
diftest$SSABIC <- c(fitmeasures(efast_fit, "bic2"), 
                    fitmeasures(efa_fit,   "bic2"))
diftest$CFI    <- c(fitmeasures(efast_fit, "cfi"), 
                    fitmeasures(efa_fit,   "cfi"))
diftest$RMSEA  <- c(fitmeasures(efast_fit, "rmsea"), 
                    fitmeasures(efa_fit,   "rmsea"))
diftest$SRMR   <- c(fitmeasures(efast_fit, "srmr"), 
                    fitmeasures(efa_fit,   "srmr"))
diftest

# latex table
xtable(diftest[, c(9:11, 1, 4:7)], digits = 3)

# plot the observed-implied covariance for both models
res_efa   <- 
  cov2cor(lavTech(efa_fit, "sampstat")[[1]]$cov) - 
  cov2cor(unclass(fitted(efa_fit)$cov))
res_efast <- 
  cov2cor(lavTech(efast_fit, "sampstat")[[1]]$cov) - 
  cov2cor(unclass(fitted(efast_fit)$cov))
srmr_efa   <- sqrt(res_efa^2)
srmr_efast <- sqrt(res_efast^2)

pdf("empirical_examples/plots/camcan_wm_obs_imp.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
cplot(srmr_efa, mar = c(0, 0, 1, 0), main = "EFA", 
      cl.lim = range(c(srmr_efa, srmr_efast)), cl.pos = "r", cex = 0.5)
cplot(srmr_efast, mar = c(0, 0, 1, 0), main = "EFAST", 
      cl.lim = range(c(srmr_efa, srmr_efast)), cl.pos = "r", cex = 0.5)
dev.off()


# plot the extracted components from the efast model
pdf("empirical_examples/plots/camcan_wm_components.pdf", width = 9, height = 3)
par(mfrow = c(1, 3))
cplot(decomposition(efast_fit)[-1], is.corr = FALSE)
dev.off()


# lateralization index
li <- lateralization(efast_fit)
li %>% 
  as_tibble %>% 
  ggplot(aes(x = str_replace_all(region, "_", " "), y = est, 
             ymax = ci.upper, ymin = ci.lower)) + 
  geom_pointrange() + 
  theme_fira() + 
  coord_flip() +
  labs(x = "", y = "Lateralization index")

firaSave("empirical_examples/plots/camcan_wm_lateralization.pdf", 
         width = 8, height = 5)

