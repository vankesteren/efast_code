# Process the results of the simulation
# Erik-Jan van Kesteren
# Last edited: 20200630
library(efast)
library(tidyverse)
library(firatheme)

tbl <- read_rds("simulations/output/sim_run_efast_complete.rds")

## Gather EFA and EFAST
tbl <- 
  tbl %>% 
  gather(key = "model", value = "object", efa_fit, efast_fit) %>% 
  mutate(model = toupper(str_extract(model, ".*?(?=_)")))


# get covariance of factors from an object
get_fac_cov <- function(obj) {
  faccov <- crossprod(obj@Model@H[[1]][1:4, 1:4])
  abscov <- abs(faccov[lower.tri(faccov)])
  mean(abscov)
}


## compute convergence, aic, bic, estimated absolute factor covariance
tbl <- 
  tbl %>% 
  mutate(
    convergence = map(object, function(li) {
      map_lgl(li, function(obj) {
        inspect(obj, what = "post.check")
      })
    }),
    prop_convergence = map_dbl(convergence, function(conv) {
      mean(conv)
    }),
    lb_convergence = map_dbl(convergence, function(conv) {
      qbeta(0.025, sum(conv) + 0.5, length(conv) - sum(conv) + 0.5)
    }),
    ub_convergence = map_dbl(convergence, function(conv) {
      qbeta(0.975, sum(conv) + 0.5, length(conv) - sum(conv) + 0.5)
    }),
    mean_BIC = map_dbl(object, function(li) {
      mean(map_dbl(li, function(obj) {
        BIC(obj)
      }))
    }),
    se_BIC = map_dbl(object, function(li) {
      sd(map_dbl(li, function(obj) {
        BIC(obj)
      }))
    }),
    mean_AIC = map_dbl(object, function(li) {
      mean(map_dbl(li, function(obj) {
        AIC(obj)
      }))
    }),
    se_AIC = map_dbl(object, function(li) {
      sd(map_dbl(li, function(obj) {
        AIC(obj)
      }))
    }),
    mean_est_cov = map_dbl(object, function(li) {
      mean(map_dbl(li, get_fac_cov))
    }),
    se_est_cov = map_dbl(object, function(li) {
      sd(map_dbl(li, get_fac_cov))
    }),
    psi_cov_error = mean_est_cov - psi_cov
  )

# compute average mae over lambda
get_lam <- function(lambda) {
  lam_lat <- lambda*0.85
  lam_bil <- lambda
  
  # first, generate true lambda
  lambda <- matrix(0, nrow = 34, ncol = 4)
  
  lambda[1:8,   1] <- lam_lat
  
  lambda[9:11,  2] <- lam_bil
  lambda[12:14, 3] <- lam_bil
  lambda[15:17, 4] <- lam_bil
  
  lambda[18:20, 2] <- lam_bil
  lambda[21:23, 3] <- lam_bil
  lambda[24:25, 4] <- lam_bil
  
  lambda[26:28, 2] <- lam_bil
  lambda[29:31, 3] <- lam_bil
 
  lambda[32:34, 4] <- lam_bil
  colnames(lambda) <- c("LH", "F1", "F2", "F3")
  lambda
}

factor_mae <- function(obj, lam, factor = 1) {
  est_lam <- obj@Model@GLIST$lambda
  mae <- function(a, b) as.numeric(mean(abs(a - b)))
  min(apply(est_lam, 2, mae, b = lam[, factor]))
}

tbl_facs <- tbl
for (i in 1:nrow(tbl)) {
  lam <- get_lam(tbl$lambda[i])
  lat_names <- c("LH", "F1", "F2", "F3")
  
  for (name in lat_names) {
    maes <- numeric(length(tbl$object[[i]]))
    for (o in seq_along(tbl$object[[i]])) {
      maes[o] <- factor_mae(tbl$object[[i]][[o]], lam = lam, factor = name)
    }
    tbl_facs[[paste0("maes_", name)]][i] <- list(maes)
  }
}

# tidy data
tbl_facs <- 
  tbl_facs %>% 
  gather(key = "factor", value = "maes", contains("maes_")) %>%
  separate(factor, c("remove", "factor")) %>% 
  select(-remove) %>% 
  mutate(mean_mae = sapply(maes, function(mae_list) mean(unlist(mae_list))),
         se_mae   = sapply(maes, function(mae_list) sd(unlist(mae_list))))


## PLOTS
## -----

pd <- position_dodge(width = 0.05)

# convergence
conv_plot <- 
  tbl %>% 
  filter(lambda == 0.5, N == 650) %>% 
  mutate(se_prop = sqrt(prop_convergence * (1 - prop_convergence)) / 120,
         upper = ub_convergence, lower = lb_convergence) %>% 
  ggplot(aes(x = res_cor, y = prop_convergence, colour = model, ymin = lower, 
             ymax = upper)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(width = 0.05, size = 1, position = pd) +
  theme_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  scale_colour_fira() +
  ylim(0, 1) +
  labs(
    title = "Model convergence",
    x = "Contralateral homology correlation", 
    y = "Convergence proportion", 
    colour = "Model type"
  ) +
  facet_wrap(~as_factor(psi_cov) %>% 
               fct_recode("Latent correlation = 0"   = "0", 
                          "Latent correlation = 0.5" = "0.5"))

firaSave("simulations/plots/conv_plot.pdf", plot = conv_plot, width = 9, 
         height = 6)

# convergence by N
conv_plot_n <- 
  tbl %>% 
  filter(psi_cov == 0.5) %>% 
  mutate(se_prop = sqrt(prop_convergence * (1 - prop_convergence)) / 120,
         upper = ub_convergence, lower = lb_convergence) %>% 
  ggplot(aes(x = res_cor, y = prop_convergence, colour = model, ymin = lower, 
             ymax = upper)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(width = 0.05, size = 1, position = pd) +
  theme_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  scale_colour_fira() +
  ylim(0, 1) +
  labs(
    title = "Model convergence",
    subtitle = "Latent covariance = 0.5",
    x = "Contralateral homology correlation", 
    y = "Convergence proportion", 
    colour = "Model type"
  ) +
  facet_grid(
    as_factor(lambda) %>% 
      fct_recode("Loadings = 0.5" = "0.5", "Loadings = 0.7" = "0.7") 
    ~ 
    as_factor(N) %>% 
      fct_recode("N = 650" = "650", "N = 130" = "130", "N = 65" = "65") %>% 
      fct_rev()
  )

firaSave("simulations/plots/conv_plot_n.pdf", plot = conv_plot_n, width = 9, 
         height = 6)

# bic
BIC_plot <- 
  tbl %>% 
  filter(lambda == 0.5, N == 650) %>% 
  mutate(upper = mean_BIC + 1.96*se_BIC,
         lower = mean_BIC - 1.96*se_BIC) %>% 
  ggplot(aes(x = res_cor, y = mean_BIC, colour = model, ymin = lower, 
             ymax = upper)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(width = 0.05, size = 1, position = pd) +
  theme_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  scale_colour_fira() +
  labs(
    title = "BIC",
    x = "Contralateral homology correlation", 
    y = "BIC", 
    colour = "Model type"
  ) +
  facet_wrap(~as_factor(psi_cov) %>% 
               fct_recode("Latent correlation = 0"   = "0", 
                          "Latent correlation = 0.5" = "0.5"))

firaSave("simulations/plots/BIC_plot.pdf", plot = BIC_plot, width = 9, 
         height = 6)


# aic
AIC_plot <- 
  tbl %>% 
  filter(lambda == 0.5, N == 650) %>% 
  mutate(upper = mean_AIC + 1.96*se_AIC,
         lower = mean_AIC - 1.96*se_AIC) %>% 
  ggplot(aes(x = res_cor, y = mean_AIC, colour = model, ymin = lower, 
             ymax = upper)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(width = 0.05, size = 1, position = pd) +
  theme_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  scale_colour_fira() +
  labs(
    title = "AIC",
    x = "Contralateral homology correlation", 
    y = "AIC", 
    colour = "Model type"
  ) +
  facet_wrap(~as_factor(psi_cov) %>% 
               fct_recode("Latent correlation = 0"   = "0", 
                          "Latent correlation = 0.5" = "0.5"))


firaSave("simulations/plots/AIC_plot.pdf", plot = AIC_plot, width = 9, 
         height = 6)


# covariance error
psi_cov_plot <- 
  tbl %>% 
  filter(psi_cov == 0.5, N == 650) %>% 
  mutate(lower = mean_est_cov - 1.96*se_est_cov, 
         upper = mean_est_cov + 1.96*se_est_cov) %>% 
  ggplot(aes(x = res_cor, y = mean_est_cov, colour = model, ymin = lower, 
             ymax = upper)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(width = 0.05, size = 1, position = pd) +
  theme_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  scale_colour_fira() +
  labs(
    title    = "Factor covariance estimates",
    subtitle = "True factor covariance: 0.5",
    x        = "Contralateral homology correlation", 
    y        = "Estimated factor covariance", 
    colour   = "Model type"
  ) +
  facet_wrap(~as_factor(lambda) %>% 
               fct_recode("Factor loading = 0.5" = "0.5",
                          "Factor loading = 0.7" = "0.7"))

firaSave("simulations/plots/psi_cov_plot.pdf", plot = psi_cov_plot, width = 9, 
         height = 6)

# Lambda mae
mae_plot <- 
  tbl_facs %>% 
  filter(lambda == 0.5, psi_cov == 0.5, N == 650) %>% 
  mutate(upper = mean_mae + 1.96*se_mae,
         lower = mean_mae - 1.96*se_mae) %>% 
  ggplot(aes(x = res_cor, y = mean_mae, colour = model, 
             ymin = lower, ymax = upper)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(width = 0.05, size = 1, position = pd) +
  facet_grid(~factor %>% fct_recode(
    "Bilateral Factor 1" = "F1",
    "Bilateral Factor 2" = "F2",
    "Bilateral Factor 3" = "F3",
    "Lateral Factor" = "LH"
  )) +
  theme_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  scale_colour_fira() +
  labs(
    title = "Estimation error of factor loadings",
    x = "Contralateral homology correlation", 
    y = "Mean absolute error", 
    colour = "Model type"
  )


firaSave("simulations/plots/mae_plot.pdf", plot = mae_plot, width = 9, 
         height = 6)

## SUPPLEMENTARY PLOTS
## -------------------

mae_plot_full <- 
  tbl_facs %>% 
  filter(N == 650) %>% 
  mutate(upper = mean_mae + 1.96*se_mae,
         lower = mean_mae - 1.96*se_mae) %>% 
  ggplot(aes(x = res_cor, y = mean_mae, colour = model, 
             ymin = lower, ymax = upper)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(width = 0.05, size = 1, position = pd) +
  facet_grid(paste("Loadings = ", lambda) + 
               paste("Covariance = ", psi_cov) ~
               factor %>% fct_recode(
                 "Bilateral Factor 1" = "F1",
                 "Bilateral Factor 2" = "F2",
                 "Bilateral Factor 3" = "F3",
                 "Lateral Factor"     = "LH"
                 ), scales = "free_y"
             ) +
  theme_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434"), 
        legend.position = "top") +
  scale_colour_fira() +
  labs(
    title = "Estimation error of factor loadings",
    x = "Contralateral homology correlation", 
    y = "Mean absolute error", 
    colour = "Model type"
  )

firaSave("simulations/plots/mae_plot_full.pdf", plot = mae_plot_full,
         width = 10, height = 12)


mae_plot_sample <- 
  tbl_facs %>% 
  filter(lambda == 0.7, psi_cov == 0.5) %>% 
  mutate(upper = mean_mae + 1.96*se_mae,
         lower = mean_mae - 1.96*se_mae) %>% 
  ggplot(aes(x = res_cor, y = mean_mae, colour = model, 
             ymin = lower, ymax = upper)) +
  geom_line(size = 1, position = pd) +
  geom_point(size = 2.5, position = pd) +
  geom_errorbar(width = 0.05, size = 1, position = pd) +
  facet_grid(as_factor(paste("N = ", N)) ~
               factor %>% fct_recode(
                 "Bilateral Factor 1" = "F1",
                 "Bilateral Factor 2" = "F2",
                 "Bilateral Factor 3" = "F3",
                 "Lateral Factor"     = "LH"
               )
  ) +
  theme_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434"), 
        legend.position = "top") +
  scale_colour_fira() +
  labs(
    title = "Estimation error of factor loadings",
    subtitle = "Loadings = 0.7, covariance = 0.5",
    x = "Contralateral homology correlation", 
    y = "Mean absolute error", 
    colour = "Model type"
  )

firaSave("simulations/plots/mae_plot_sample.pdf", plot = mae_plot_sample,
         width = 10, height = 9)
