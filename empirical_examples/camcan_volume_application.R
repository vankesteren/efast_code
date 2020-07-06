# EFAST model on empirical camcan data with final model
# Determining number of factors and performing likelihood ratio tests
# Erik-Jan van Kesteren
# Last edited: 20200430
library(efast)
library(parallel)
library(pbapply)
library(tidyverse)
library(firatheme)
library(qgraph)
library(gtools)
library(ggseg)

# load the dataset
volume <- 
  read_rds("empirical_examples/data/camcan_volume.rds") %>% 
  select(-contains("unknown")) %>% # remove the unkown vars
  mutate_all(scale) %>% # mean 0, sd 1 variables
  set_names(str_extract_all(colnames(.), "(?<=ctx-).*")) %>% 
  set_names(str_replace_all(colnames(.), "-", "_")) # preprocess names
roi_names <- c(unique(str_extract_all(colnames(volume), "(?<=_).*", TRUE)))

# plot the covmat
pdf("empirical_examples/plots/camcan_volume_cplot.pdf", 7, 7)
cplot(cor(volume), cl.pos = "r")
dev.off()


#####################################
### Determining number of factors ###
#####################################
# create a grid of efast / efa parameters
grid <- 
  expand.grid(method = c("EFA", "EFAST"), M = 2:20) %>% 
  as_tibble()

# create a model fit function
fit_model <- function(method, M) {
  if (method == "EFAST") {
    fit <- try(efast_hemi(data = volume, M = M, lh_idx = 1:34, rh_idx = 35:68, 
                          roi_names = roi_names, optim.force.converged = TRUE, 
                          se = "none", rotation = "none", 
                          control = list(iter.max = 300)))
  } else {
    fit <- try(efast_efa(data = volume, M = M, optim.force.converged = TRUE, 
                         se = "none", rotation = "none",
                         control = list(iter.max = 300)))
  }
}

# create a parallel cluster
clus <- makeCluster(12)
clusterExport(clus, c("volume", "roi_names", "grid", "fit_model"))
clusterEvalQ(clus, library(efast))

# fit the models
fits <- pblapply(1:nrow(grid), function(i) {
  fit_model(grid[i,]$method, grid[i,]$M)
}, cl = clus)

stopCluster(clus)

grid$fit <- fits

# model comparison plot
grid$valid <- sapply(grid$fit, inspect, "post.check")
grid <- 
  grid %>% 
  mutate(
    BIC     = sapply(fit, BIC),
    AIC     = sapply(fit, AIC),
    SSABIC  = sapply(fit, fitmeasures, fit.measures = "bic2")
  )

grid %>% 
  gather(key = "Metric", value = "Value", BIC, AIC, SSABIC) %>%
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
  labs(colour = "Model type", y = "", x = "Number of EFA factors") +
  facet_wrap(~Metric) +
  theme(legend.title.align = 0) +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  guides(colour = guide_legend(order = 1), 
         alpha  = guide_legend(order = 2),
         size   = guide_legend(order = 2))

firaSave("empirical_examples/plots/camcan_volume_fit.pdf", width = 8, 
         height = 4)

########################################
### Fitting and comparing best model ###
########################################

# Fit the efast and efa models with 6 factors
fit_efast <- efast_hemi(volume, 6, 1:34, 35:68, roi_names, verbose = TRUE)
fit_efa   <- efast_efa(volume, M = 6, verbose = TRUE)

# Extract lambda
lambda_efast <- efast_loadings(fit_efast)
lambda_efa   <- efast_loadings(fit_efa)

# reorder lambda to be equal across models
perms <- permutations(6, 6)
mses  <- apply(perms, 1, function(p) sum((lambda_efast - lambda_efa[,p])^2))
# bonus plot
plot(mses, pch = 21, bg = "#45454588", bty = "n", axes = F, 
     xlab = "Permutation", ylab = "MSE", cex = min(mses)/mses)
lambda_efa <- lambda_efa[,perms[which.min(mses),]]


# save it all
save(lambda_efa, lambda_efast, fit_efa, fit_efast, 
     file = "empirical_examples/output/camcan_volume.Rdata")

# create pretty print table
efast_tab <- round(lambda_efast, 2)
efa_tab   <- round(lambda_efa, 2)
efast_tab[abs(lambda_efast) < 0.1] <- ""
efa_tab[abs(lambda_efa) < 0.1]     <- ""

xtable::xtable(cbind(efast_tab, efa_tab))

# Network plots of latent covariance
pdf("R/empirical_examples/plots/camcan_volume_network_latents.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
Psi_efast <- crossprod(fit_efast@Model@H[[1]][1:6, 1:6])
Psi_efa   <- crossprod(fit_efa@Model@H[[1]][perms[which.min(mses),], 
                                            perms[which.min(mses),]])

qgraph(Psi_efa - diag(6),   layout = "spring", posCol = firaCols[1], title = 'A')
qgraph(Psi_efast - diag(6), layout = "spring", posCol = firaCols[1], title = 'B')
dev.off()

# compare efast and efa using AIC, BIC, SSABIC, LR-test
diftest <- anova(fit_efa, fit_efast)
rownames(diftest) <- c("EFAST", "EFA")
diftest$SSABIC <- c(fitmeasures(fit_efast, "bic2"), 
                    fitmeasures(fit_efa,   "bic2"))
diftest$CFI    <- c(fitmeasures(fit_efast, "cfi"), 
                    fitmeasures(fit_efa,   "cfi"))
diftest$RMSEA  <- c(fitmeasures(fit_efast, "rmsea"), 
                    fitmeasures(fit_efa,   "rmsea"))
diftest$SRMR   <- c(fitmeasures(fit_efast, "srmr"), 
                    fitmeasures(fit_efa,   "srmr"))
diftest

# this goes in the paper
xtable::xtable(diftest[, c(9:11, 1, 4:7)], digits = 3)



# plot the extracted components from the efast model
pdf("empirical_examples/plots/camcan_volume_components.pdf", width = 9, height = 3)
par(mfrow = c(1, 3))
lapply(decomposition(fit_efast)[-1], corrplot, method = "shade", tl.pos = "n", 
       cl.pos = "n")
dev.off()


# lateralization index
li <- lateralization(fit_efast)
li <- li %>% as_tibble
dkt_li <- dkt %>% 
  filter(hemi == "left") %>% 
  mutate(region = str_extract(label, "(?<=_).*")) %>% 
  left_join(li, by = "region")
  

ggseg(dkt_li, mapping = aes(fill = est), 
      hemisphere = "left", colour = "white", adapt_scales = TRUE, position =  "stacked") +
  scale_fill_viridis_c() +
  theme_fira() +
  theme(line = element_blank(),
        axis.text.y = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank()) +
  labs(x = "", y = "",
       fill  = "Lateralization")

firaSave("empirical_examples/plots/camcan_volume_lateralization.pdf", width = 8, height = 3)

