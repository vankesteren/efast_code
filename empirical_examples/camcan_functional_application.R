# CamCAN functional data empirical analysis with EFAST
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

# load the data: list of functional covariance matrices for 5 participants
# we use the first participant at first
clist <- read_rds("empirical_examples/data/camcan_fun.rds")
ssize <- 261 # sample size
ccfun <- clist[[1]]

# plot the covmat
pdf("empirical_examples/plots/camcan_fun_cplot.pdf", 7, 7)
cplot(cor(ccfun), cl.pos = "r")
dev.off()

# set area names
area_names <- paste0("ROI", 1:45)


#####################################
### Determining number of factors ###
#####################################

# create a grid of efast / efa parameters
results <- 
  expand.grid(method = c("EFA", "EFAST"), M = 2:16) %>% 
  as_tibble()

# create a model fit function
fit_model <- function(method, M) {
  if (method == "EFAST") {
    fit <- try(efast_hemi(data = ccfun, sample.nobs = ssize, M = M, 
                          lh_idx = 1:45, rh_idx = 46:90, 
                          roi_names = area_names))
  } else {
    fit <- try(efast_efa(data = ccfun, M = M, sample.nobs = ssize))
  }
  fit
}

# create a parallel cluster
clus <- makeCluster(8)
clusterExport(clus, c("ccfun", "area_names", "results", "fit_model", "ssize"))
clusterEvalQ(clus, library(efast))

# fit the models
fits <- pblapply(1:nrow(results), function(i) {
  fit_model(results[i,]$method, results[i,]$M)
}, cl = clus)

stopCluster(clus)

# add them to the results
results$fit <- fits

# and save the results
write_rds(results, "empirical_examples/output/camcan_fun_fit.rds")

# optionally, load the checkpoint
results <- read_rds("empirical_examples/output/camcan_fun_fit.rds")

# extract fit measures
results_metrics <-
  results %>% 
  mutate(
    AIC    = sapply(fit, fitmeasures, fit.measures = "aic"),
    BIC    = sapply(fit, fitmeasures, fit.measures = "bic"),
    SSABIC = sapply(fit, fitmeasures, fit.measures = "bic2"),
    valid  = sapply(fit, lavaan:::lav_object_post_check)
  ) %>% 
  select(-fit) %>% 
  gather("Metric", "Value", -method, -M, -valid)

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
       main = "CamCAN functional network factor analysis metrics") +
  facet_wrap(~Metric) +
  theme(legend.title.align = 0) +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  guides(colour = guide_legend(order = 1), 
         alpha  = guide_legend(order = 2),
         size   = guide_legend(order = 2))

firaSave("R/empirical_examples/plots/camcan_fun_metrics.pdf", 
         plot = metrics_plot, width = 8, height = 4)

########################################
### Fitting and comparing best model ###
########################################

# load the checkpoint
results <- read_rds("empirical_examples/output/camcan_fun_fit.rds")

# best model is the EFAST model with M = 13
efast_fit <- results %>% filter(method == "EFAST", M == 13) %>% .$fit %>% .[[1]]
efa_fit   <- results %>% filter(method == "EFA",   M == 13) %>% .$fit %>% .[[1]]

# Network plots of latent covariance
pdf("empirical_examples/plots/camcan_fun_network_latents.pdf", 
    width = 12, height = 6)
par(mfrow = c(1,2))
Psi_efast <- crossprod(efast_fit@Model@H[[1]])
Psi_efa   <- crossprod(efa_fit@Model@H[[1]])
qgraph(Psi_efa - diag(13),   layout = "spring", posCol = firaCols[1], 
       title = 'A')
qgraph(Psi_efast - diag(13), layout = "spring", posCol = firaCols[1], 
       title = 'B')
dev.off()


# compare efast and efa using AIC, BIC, SSABIC, LR-test
diftest <- anova(efast_fit, efa_fit)
rownames(diftest) <- c("EFAST", "EFA")
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
res_efa <- 
  cov2cor(lavTech(efa_fit, "sampstat")[[1]]$cov) - 
  cov2cor(unclass(fitted(efa_fit)$cov))
res_efast <- 
  cov2cor(lavTech(efast_fit, "sampstat")[[1]]$cov) - 
  cov2cor(unclass(fitted(efast_fit)$cov))
srmr_efa   <- sqrt(res_efa^2)
srmr_efast <- sqrt(res_efast^2)

pdf("empirical_examples/plots/camcan_fun_obs_imp.pdf", width = 6, height = 3)
par(mfrow = c(1, 2))
cplot(srmr_efa, mar = c(0, 0, 1, 0), main = "EFA", 
      cl.lim = range(c(srmr_efa, srmr_efast)))
cplot(srmr_efast, mar = c(0, 0, 1, 0), main = "EFAST", 
      cl.lim = range(c(srmr_efa, srmr_efast)))
dev.off()


# plot the extracted components from the efast model
pdf("empirical_examples/plots/camcan_fun_components.pdf", width = 9, 
    height = 3)
par(mfrow = c(1, 3))
cplot(decomposition(efast_fit)[-1], is.corr = FALSE)
dev.off()


# lateralization index
li <- lateralization(efast_fit)
li %>% 
  as_tibble %>% 
  mutate(
    roi = region %>% str_replace_all("_", "") %>% str_pad(2, "left", "0")
  ) %>% 
  ggplot(aes(x = as_factor(roi), y = est, ymax = ci.upper, ymin = ci.lower)) + 
  geom_pointrange() + 
  theme_fira() + 
  coord_flip() +
  labs(x = "", y = "Lateralization index")
firaSave("empirical_examples/plots/camcan_fun_lateralization.pdf", width = 6, 
         height = 10)

###############################################
### Comparing the model across participants ###
###############################################

# create a group fit with 5 participants (NB: takes a long time!)
hemi_syntax <- efast_fit@external$syntax
group_fit <- lavaan(
  model          = hemi_syntax,
  sample.cov     = clist, 
  sample.nobs    = rep(261, 5),
  auto.fix.first = FALSE, 
  auto.var       = TRUE, 
  auto.efa       = TRUE, 
  information    = "observed",
  verbose        = TRUE,
  se             = "none"
)

write_rds(group_fit, "R/empirical_examples/output/fun_group_fit.rds")
group_fit <- read_rds("R/empirical_examples/output/fun_group_fit.rds")
pe <- parameterestimates(group_fit)


pe %>% 
  filter(str_detect(lhs, "F[123]$"), str_ends(rhs, "ROI[0-9]+")) %>% 
  mutate(hemi = as_factor(substr(rhs, 1, 2))) %>% 
  group_by(group, lhs) %>% 
  mutate(ROI = row_number()) %>% 
  ggplot(aes(x = ROI, y = est, colour = as_factor(group), shape = hemi)) +
  geom_point(position =  position_dodge(width = 0.1)) +
  geom_line() +
  facet_grid(cols = vars(lhs), rows = vars(group)) +
  theme_fira() +
  scale_x_continuous(breaks = c(1, 45, 90)) +
  scale_colour_fira(guide = FALSE) +
  scale_shape_manual(values = c(lh = "square", rh = "circle"), guide = FALSE) +
  labs(x = "Region of interest (lh, rh)", y = "Factor loading") +
  theme(panel.background = element_rect(colour = "#343434"))
  
firaSave("empirical_examples/plots/camcan_fun_comparison.pdf", width = 9, 
         height = 7)
