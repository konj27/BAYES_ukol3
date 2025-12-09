library(bayesrules)
library(bayesplot)
library(tidyverse)
library(psych)
library(ggcorrplot)
library(glue)
library(rstanarm)
library(broom.mixed)
library(tidybayes)
library(loo)
# ukol 1 Honza Petr

# ukol 2 Jenda Anička
# 11.10 ---------------------------------------------------
penguin_data <- penguins_bayes %>% 
  filter(species %in% c("Adelie", "Gentoo")) %>% na.omit()
plot(x = penguin_data$body_mass_g, y = penguin_data$flipper_length_mm)
ggplot(data = penguin_data)+
  geom_point(mapping = aes(x=body_mass_g, y=flipper_length_mm,color = species))+
  ylab("flipper length (mm)")+
  xlab("body mass (g)")

summarize_data <- penguin_data %>% 
  reframe(flipper_length_mm, body_mass_g, species)%>%
  na.omit()

summary(summarize_data[summarize_data$species == "Adelie", c("flipper_length_mm","body_mass_g")])
summary(summarize_data[summarize_data$species == "Gentoo", c("flipper_length_mm","body_mass_g")])

cor(summarize_data[summarize_data$species == "Adelie", c("flipper_length_mm","body_mass_g")])
cor(summarize_data[summarize_data$species == "Gentoo", c("flipper_length_mm","body_mass_g")])

R2_mod = lm(body_mass_g ~ flipper_length_mm, summarize_data)
R2 <- summary(R2_mod)$r.squared

glue("R^2 = {R2}")

penguin_model = stan_glm(
  body_mass_g ~ flipper_length_mm + species,
  data = summarize_data,
  family = gaussian,
  prior_intercept = normal(5000, 350),
  prior = normal(5000, 350, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

prior_summary(penguin_model)
mcmc_trace(penguin_model)
mcmc_dens_overlay(penguin_model)
mcmc_acf(penguin_model)

pp_check(penguin_model)
tidy(penguin_model, effects = c("fixed", "aux"))
df_predict <- data.frame(flipper_length_mm = 197, species = "Adelie")

predictions_1 <- posterior_predict(penguin_model, newdata = df_predict)

mcmc_hist(predictions_1) + labs(x ="Y") 
table(predictions_1)
colMeans(predictions_1)
#11.11 ------------------------------------------------------------------

penguin_model_inter = stan_glm(
  body_mass_g ~ flipper_length_mm * species,
  data = summarize_data,
  family = gaussian,
  prior_intercept = normal(5000, 350),
  prior = normal(5000, 350, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

pp_check(penguin_model_inter)
tidy(penguin_model_inter, effects = c("fixed", "aux"))

summarize_data %>%
  add_fitted_draws(penguin_model_inter, n = 50) %>%
  ggplot(aes(x = flipper_length_mm, y = body_mass_g, color = species)) +
  geom_line(aes(y = .value, group = paste(species, .draw)), alpha = .1) +
  geom_point(data = summarize_data, size = 0.5)+
  ylab("body mass (g)")+
  xlab("flipper length (mm)")
# 11.12 ------------------------------------------------------------
penguin_data_2 <- penguin_data %>% 
  reframe(body_mass_g, flipper_length_mm, bill_length_mm, bill_depth_mm) %>%
  na.omit()

penguin_model_2 = stan_glm(
  body_mass_g ~ flipper_length_mm + bill_length_mm + bill_depth_mm,
  data = penguin_data_2,
  family = gaussian,
  prior_intercept = normal(5000, 350),
  prior = normal(5000, 350, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

posterior_interval(penguin_model_2, prob = 0.95)
# 11.13 ------------------------------------------------------------
penguin_model_3 <- stan_glm(
  body_mass_g ~ flipper_length_mm,
  data = penguin_data_2,
  family = gaussian,
  prior_intercept = normal(5000, 350),
  prior = normal(5000, 350, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735
)
pp_check(penguin_model_3)

penguin_model_4 <- stan_glm(
  body_mass_g ~ species,
  data = penguin_data,
  family = gaussian,
  prior_intercept = normal(5000, 350),
  prior = normal(5000, 350, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735
)
pp_check(penguin_model_4)

pp_check(penguin_model, plotfun = "dens_overlay", nreps = 50)
pp_check(penguin_model_2, plotfun = "dens_overlay", nreps = 50)
pp_check(penguin_model_3, plotfun = "dens_overlay", nreps = 50)
pp_check(penguin_model_4, plotfun = "dens_overlay", nreps = 50)
penguins_complete <- penguins_bayes %>% 
  select(flipper_length_mm, body_mass_g, species, 
         bill_length_mm, bill_depth_mm) %>% 
  na.omit() 

cv_1 <- prediction_summary_cv(
  model = penguin_model,
  data  = penguins_complete,
  k     = 10
)

cv_2 <- prediction_summary_cv(
  model = penguin_model_2,
  data  = penguins_complete,
  k     = 10
)

cv_3 <- prediction_summary_cv(
  model = penguin_model_3,
  data  = penguins_complete,
  k     = 10
)

cv_4 <- prediction_summary_cv(
  model = penguin_model_4,
  data  = penguins_complete,
  k     = 10
)

cv_1
cv_2
cv_3
cv_4


loo_1 <- loo(penguin_model)
loo_2 <- loo(penguin_model_2)
loo_3 <- loo(penguin_model_3)
loo_4 <- loo(penguin_model_4)

loo_1$estimates["elpd_loo", ]
loo_2$estimates["elpd_loo", ]
loo_3$estimates["elpd_loo", ]
loo_4$estimates["elpd_loo", ]

# Direct comparison
loo_compare(loo_1, loo_2, loo_3, loo_4)

# 11.14 ------------------------------------------------------------

penguins_numeric <- penguins_bayes %>% 
  select(flipper_length_mm, body_mass_g, 
         bill_length_mm, bill_depth_mm) %>% 
  na.omit() 
cor(penguins_numeric)
ggcorrplot(cor(penguins_numeric))

ggplot(data = penguin_data)+
  geom_point(mapping = aes(x=bill_length_mm, y=body_mass_g,color = sex))+
  ylab("body mass (g)")+
  xlab("bill length (mm)")

ggplot(data = penguin_data)+
  geom_point(mapping = aes(x=bill_length_mm, y=body_mass_g,color = island))+
  ylab("body mass (g)")+
  xlab("bill length (mm)")

ggplot(data = penguin_data)+
  geom_point(mapping = aes(x=bill_length_mm, y=body_mass_g,color = species))+
  ylab("body mass (g)")+
  xlab("bill length (mm)")

penguin_our_model_1 = stan_glm(
  bill_length_mm ~ flipper_length_mm + body_mass_g+ island,
  data = penguin_data,
  family = gaussian,
  prior_intercept = normal(40, 15),
  prior = normal(0, 2.5, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

prior_summary(penguin_our_model_1)
mcmc_trace(penguin_our_model_1)
mcmc_dens_overlay(penguin_our_model_1)
mcmc_acf(penguin_our_model_1)
pp_check(penguin_our_model_1)

penguin_our_model_2 = stan_glm(
  bill_length_mm ~ flipper_length_mm + body_mass_g + sex,
  data = penguin_data,
  family = gaussian,
  prior_intercept = normal(40, 15),
  prior = normal(0, 2.5, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

pp_check(penguin_our_model_2)

penguin_our_model_3 = stan_glm(
  bill_length_mm ~ flipper_length_mm + body_mass_g * island,
  data = penguin_data,
  family = gaussian,
  prior_intercept = normal(40, 15),
  prior = normal(0, 2.5, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

pp_check(penguin_our_model_3)

# ukol 3 Honza Petr

# ukol 4 Jenda Anička