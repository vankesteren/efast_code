# Process the results of the overextraction simulation
# Erik-Jan van Kesteren
# Last edited: 20200430
library(efast)
library(tidyverse)
library(firatheme)


tbl <- read_rds("simulations/output/overextraction.rds")


# look at factor extraction by AIC and BIC criterion
tbl_tidy <- 
  tbl %>% 
  gather("model", "fit", efa_fit, efast_fit) %>% 
  mutate(
    aic    = map_dbl(fit, AIC),
    bic    = map_dbl(fit, BIC),
    ssabic = map_dbl(fit, function(x) fitmeasures(x)["bic2"]),
  ) %>% 
  gather("criterion", "value", aic, bic, ssabic) %>%
  mutate(
    model     = toupper(str_extract(model, ".*(?=_)")),
    criterion = toupper(criterion)
  )

tbl_tidy %>% 
  ggplot(aes(x = M, y = value, colour = model, fill = model, alpha = dataset)) +
  geom_vline(xintercept = 4, lty = "dashed") +
  geom_point(alpha = 0.1, position = position_jitter(width = 0.1), pch = 21) +
  geom_line(data = tbl_tidy %>% group_by(M, model, criterion) %>% 
              summarise(value = mean(value)), 
            aes(alpha = NULL), size = 1) +
  facet_wrap(~criterion) + 
  theme_fira() +
  scale_colour_fira() +
  scale_fill_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434")) +
  labs(x = "Number of factors", y = "Criterion value", colour = "Model", fill = "Model")

firaSave("simulations/plots/extraction_lines.pdf", width = 9, height = 6)

tbl_tidy %>% 
  group_by(dataset, criterion, model) %>% 
  summarise(nfac = which.min(value) + 1) %>% 
  ggplot(aes(x = model, y = nfac, fill = model, colour = model)) +
  geom_hline(yintercept = 4, lty = "dashed") +
  geom_point(position = position_jitter(width = 0.15, height = 0.1), 
             alpha = 0.4) +
  geom_boxplot(outlier.alpha = 0, colour = "black", fill = "transparent", 
               width = 0.5) +
  facet_grid(~criterion) +
  scale_y_continuous(breaks = 1:11, limits = c(1, 11)) +
  theme_fira() +
  scale_fill_fira() +
  scale_colour_fira() +
  theme(panel.border = element_rect(fill = NA, colour = "#343434"),
        legend.position = "none") +
  labs(x = "", y = "Extracted number of factors")

firaSave("simulations/plots/extraction_boxes.pdf", width = 9, height = 6)


# Let's look at convergence
tbl_tidy_2 <- 
  tbl %>% 
  gather("model", "fit", efa_fit, efast_fit) %>% 
  mutate(
    model       = toupper(str_extract(model, ".*(?=_)")),
    convergence = map_lgl(fit, inspect, "post.check")
  ) %>% 
  group_by(M, model) %>% 
  summarise(conv_prop = mean(convergence), 
            lb = qbeta(0.025, sum(convergence) + .5, 
                       length(convergence) - sum(convergence) + .5),
            ub = qbeta(0.975, sum(convergence) + .5, 
                       length(convergence) - sum(convergence) + .5))
# 95% Jeffrey's interval

pd <- position_dodge(0.2)
tbl_tidy_2 %>% 
  ggplot(aes(x = M, y = conv_prop, colour = model, ymin = lb, ymax = ub)) +
  geom_line(size = 1, position = pd) +
  geom_errorbar(width = 0.2, size = 1, position = pd) +
  geom_point(position = pd, size = 2.5) +
  ylim(0, 1) +
  theme_fira() + 
  scale_colour_fira() +
  labs(y = "Convergence proportion", x = "Number of factors", colour = "Model")


firaSave("simulations/plots/extraction_converge.pdf", width = 9, height = 6)
