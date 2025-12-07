library(bayesrules)
library(bayesplot)
library(tidyverse)
library(psych)
library(ggcorrplot)
library(glue)
library(rstanarm)
library(broom.mixed)
library(tidybayes)

# ukol 1 Honza Petr

# ukol 2 Jenda Anička
# 11.10 ---------------------------------------------------
penguin_data <- penguins_bayes %>% 
  filter(species %in% c("Adelie", "Gentoo"))
plot(x = penguin_data$body_mass_g, y = penguin_data$flipper_length_mm)
ggplot(data = penguin_data)+
  geom_point(mapping = aes(x=body_mass_g, y=flipper_length_mm,color = species))

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

# ukol 3 Honza Petr

# ukol 4 Jenda Anička