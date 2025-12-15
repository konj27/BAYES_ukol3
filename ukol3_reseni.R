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
library(e1071)
# ukol 1 Honza Petr

# Načtení knihoven
library(bayesrules)
library(rstanarm)
library(tidyverse)
library(janitor)

# Načtení dat a omit pro chyběhící hodnoty
data(penguins_bayes)
df <- penguins_bayes %>% na.omit()

# 9.16 a)
prior_model <- stan_glm(
  flipper_length_mm ~ bill_length_mm,
  data = df,
  family = gaussian,
  prior_intercept = normal(200, 25),
  prior = normal(0, 10),             
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 10000, seed = 84735,
  prior_PD = TRUE # vzorkovani z prioru
)

# 9.16 b)
prior_summary(prior_model)

# 9.16 c)
df %>%
  add_fitted_draws(prior_model, n = 100) %>%
  ggplot(aes(x = bill_length_mm, y = flipper_length_mm)) +
  geom_line(aes(y = .value, group = .draw), alpha = 0.2) +
  labs(x = "Délka zobáku (mm)", y = "Délka ploutve (mm)")

# 9.17 a)
ggplot(df, aes(x = bill_length_mm, y = flipper_length_mm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Délka zobáku (mm)", y = "Délka ploutve (mm)")

# 9.18 a)
post_model <- update(prior_model, prior_PD = FALSE)

# 9.18 b)
df %>%
  add_fitted_draws(post_model, n = 100) %>%
  ggplot(aes(x = bill_length_mm, y = flipper_length_mm)) +
  geom_line(aes(y = .value, group = .draw), alpha = 0.2, color = "blue") +
  geom_point(data = df)

# 9.18 c)
tidy(post_model, conf.int = TRUE, conf.level = 0.90)

# Získání řetězců parametrů
chains <- as.data.frame(post_model)

# 9.19 a) Výpočet predikcí manuálně
pablo_predictions <- chains %>%
  mutate(
    mu = `(Intercept)` + `bill_length_mm` * 51,
    y_new = rnorm(n(), mean = mu, sd = sigma)
  )

# 9.19 b)
ggplot(pablo_predictions) +
  geom_density(aes(x = mu), fill = "blue", alpha = 0.5) +
  geom_density(aes(x = y_new), fill = "orange", alpha = 0.5) +
  labs(x = "Délka ploutve (mm)")

# 9.19 c)
quantile(pablo_predictions$y_new, probs = c(0.1, 0.9))

# 9.19 e)
predict_check <- posterior_predict(post_model, newdata = data.frame(bill_length_mm = 51))
quantile(predict_check, probs = c(0.1, 0.9))

# 9.20 b) Scatterplot
ggplot(df, aes(x = body_mass_g, y = flipper_length_mm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Hmotnost (g)", y = "Délka ploutve (mm)")

# 9.20 d)
model_mass <- stan_glm(
  flipper_length_mm ~ body_mass_g,
  data = df,
  family = gaussian,
  prior = normal(0.01, 0.002),
  prior_intercept = normal(0, 100, autoscale = TRUE),
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 10000, seed = 84735
)

# 9.20 e)
model_mass %>%
  as.data.frame() %>%
  ggplot(aes(x = body_mass_g)) +
  geom_density(fill = "purple", alpha = 0.6) +
  stat_function(fun = dnorm, args = list(mean = 0.01, sd = 0.002),
                color = "red", linetype = "dashed") +
  labs(x = "Beta 1 (vliv hmotnosti)") +
  xlim(0.005, 0.017)

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

penguin_our_model_data <- penguins_bayes %>% na.omit()

penguins_numeric <- penguin_our_model_data %>% 
  select(flipper_length_mm, body_mass_g, 
         bill_length_mm, bill_depth_mm)

ggcorrplot(cor(penguins_numeric)) #nahled na korelace pro numericke promenne

#par pruzkumnych grafu
ggplot(data = penguin_our_model_data)+
  geom_point(mapping = aes(x=flipper_length_mm, y=bill_length_mm,color = species))+
  xlab("flipper length (mm)")+
  ylab("bill length (mm)")

ggplot(data = penguin_our_model_data)+
  geom_point(mapping = aes(x=flipper_length_mm, y=bill_length_mm,color = island))+
  xlab("flipper length (mm)")+
  ylab("bill length (mm)")

ggplot(data = penguin_our_model_data)+
  geom_point(mapping = aes(x=flipper_length_mm, y=bill_length_mm,color = sex))+
  xlab("flipper length (mm)")+
  ylab("bill length (mm)")

#model 1 - pouze kvantitativni promenne
penguin_our_model_1 = stan_glm(
  bill_length_mm ~ flipper_length_mm + body_mass_g,
  data = penguin_our_model_data,
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

#model 2 - pridani jedne dulezite kategoricke promenne
penguin_our_model_2 = stan_glm(
  bill_length_mm ~ flipper_length_mm + body_mass_g + island,
  data = penguin_our_model_data,
  family = gaussian,
  prior_intercept = normal(40, 15),
  prior = normal(0, 2.5, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

prior_summary(penguin_our_model_2)
mcmc_trace(penguin_our_model_2)
mcmc_dens_overlay(penguin_our_model_2)
mcmc_acf(penguin_our_model_2)
pp_check(penguin_our_model_2)

#model 3 - pridani interakcniho clenu
penguin_our_model_3 = stan_glm(
  bill_length_mm ~ body_mass_g + flipper_length_mm * species,
  data = penguin_our_model_data,
  family = gaussian,
  prior_intercept = normal(40, 15),
  prior = normal(0, 2.5, autoscale = TRUE), 
  prior_aux = exponential(1, autoscale = TRUE),
  chains = 4, iter = 5000*2, seed = 84735)

prior_summary(penguin_our_model_3)
mcmc_trace(penguin_our_model_3)
mcmc_dens_overlay(penguin_our_model_3)
mcmc_acf(penguin_our_model_3)
pp_check(penguin_our_model_3)

tidy(penguin_our_model_1, effects = c("fixed", "aux"), conf.int = TRUE, conf.level = 0.80)
tidy(penguin_our_model_2, effects = c("fixed", "aux"), conf.int = TRUE, conf.level = 0.80)
tidy(penguin_our_model_3, effects = c("fixed", "aux"), conf.int = TRUE, conf.level = 0.80)


loo_bill1 <- loo(penguin_our_model_1)
loo_bill2 <- loo(penguin_our_model_2)
loo_bill3 <- loo(penguin_our_model_3)

loo_bill1$estimates[1]
loo_bill2$estimates[1]
loo_bill3$estimates[1]

# porovnani modelu - serazeni od nejlepsiho
loo_compare(loo_bill1, loo_bill2, loo_bill3)

# ukol 3 Honza Petr
#AirBnB

if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesrules, tidyverse, rstanarm, bayesplot, broom.mixed, patchwork)

set.seed(2025)

data("airbnb_small")

#Kontrola průměru a rozptylu (Variance >> Mean indikuje overdispersion)
stats <- airbnb_small %>% 
  summarize(prumer_reviews = mean(reviews), 
            rozptyl_reviews = var(reviews))
print(stats)

#Vizualizace vztahu 
plot_boxplot <- ggplot(airbnb_small, aes(x = room_type, y = reviews, fill = room_type)) +
  geom_boxplot(alpha = 0.6) +
  labs(title = "Rozdělení recenzí podle typu pokoje",
       y = "Počet recenzí", x = "Typ pokoje") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("eda_boxplot.png", plot = plot_boxplot, width = 6, height = 4)

#Model 1: Poisson, Mean = Variance => neplatí
model_poisson <- stan_glm(
  reviews ~ rating + district + room_type + accommodates,
  data = airbnb_small,
  family = poisson(link = "log"),
  prior = normal(0, 2.5, autoscale = TRUE),
  prior_intercept = normal(0, 2.5, autoscale = TRUE),
  chains = 4, iter = 2000, refresh = 0, seed = 2025
)

#Model 2: Negativně Binomický, řeší overdispersion
model_negbin <- stan_glm(
  reviews ~ rating + district + room_type + accommodates,
  data = airbnb_small,
  family = neg_binomial_2(link = "log"),
  prior = normal(0, 2.5, autoscale = TRUE),
  prior_intercept = normal(0, 2.5, autoscale = TRUE),
  chains = 4, iter = 2000, refresh = 0, seed = 2025
)

#Posterior Predictive Check 
ppc_pois <- pp_check(model_poisson) + labs(title = "Poisson ", subtitle = "Model nepokrývá variabilitu")
ppc_nb   <- pp_check(model_negbin) + labs(title = "Negativně Binomický", subtitle = "Model sedí na data")

#Spojení grafů
final_ppc <- ppc_pois / ppc_nb
ggsave("airbnb_ppc_check.png", plot = final_ppc, width = 8, height = 8)

#Porovnání pomocí Leave-One-Out Cross-Validation
loo_poisson <- loo(model_poisson)
loo_negbin  <- loo(model_negbin)

#Vypíšeme srovnání 
comp <- loo_compare(loo_poisson, loo_negbin)
print(comp)

#Tabulka koeficientů pro vítězný modelß
final_results <- tidy(model_negbin, exponentiate = TRUE, conf.int = TRUE)
print(final_results)


# ukol 4 Jenda Anička

#priprava dat
vaccines_data <- pulse_of_the_nation %>%
  select(vaccines_are_safe, age, education, party, climate_change,science_is_honest,income) %>%
  drop_na()%>%
  mutate(
    vaccines_are_safe = factor(vaccines_are_safe,levels = c(
      "Strongly Agree",
      "Somewhat Agree",
      "Neither Agree nor Disagree",
      "Somewhat Disagree",
      "Strongly Disagree")),
    education = factor(education, levels = c(
      "Graduate degree",
      "College degree",
      "Some college",
      "High school",
      "Other")),
    science_is_honest = factor(science_is_honest,levels = c(
      "Strongly Agree",
      "Somewhat Agree",
      "Neither Agree nor Disagree",
      "Somewhat Disagree",
      "Strongly Disagree")),
    party = factor(party),
    climate_change = factor(climate_change),
    age = as.numeric(age)
)

head(vaccines_data)


#pruzkum dat - vizualizace
ggplot(vaccines_data, 
       aes(fill = education, x = vaccines_are_safe)) + 
  geom_bar(position = "fill")

ggplot(vaccines_data, 
       aes(fill = party, x = vaccines_are_safe)) + 
  geom_bar(position = "fill")

ggplot(vaccines_data, 
       aes(fill = climate_change, x = vaccines_are_safe)) + 
  geom_bar(position = "fill")

ggplot(vaccines_data, 
       aes(fill = science_is_honest, x = vaccines_are_safe)) + 
  geom_bar(position = "fill")

vaccines_data %>%
  ggplot(aes(x = age)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ vaccines_are_safe)

vaccines_data %>%
  ggplot(aes(x = income)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ vaccines_are_safe)

#naivni modely
#model1 - pouziti science_is_honest
naive_model_1 <- naiveBayes(vaccines_are_safe ~ age + education + party + income, data = vaccines_data)
#model2 - pouziti climate_change
naive_model_2 <- naiveBayes(vaccines_are_safe ~ age + climate_change + science_is_honest, data = vaccines_data)

#predikce
vaccines_pred1 <- predict(naive_model_1, newdata = vaccines_data, type = "class")
vaccines_pred2 <- predict(naive_model_2, newdata = vaccines_data, type = "class")

set.seed(84735)

# CV pro model 1
cv_vaccines_1 <- naive_classification_summary_cv(
  model = naive_model_1,
  data  = vaccines_data,
  y     = "vaccines_are_safe",
  k     = 10
)
cv_vaccines_1$cv
cv_vaccines_1$folds


# CV pro model 2
cv_vaccines_2 <- naive_classification_summary_cv(
  model = naive_model_2,
  data  = vaccines_data,
  y     = "vaccines_are_safe",
  k     = 10
)
cv_vaccines_2$cv
cv_vaccines_2$folds


#presnost obou modelu
acc1 <- mean(vaccines_data$vaccines_are_safe == vaccines_pred1) 
acc2 <- mean(vaccines_data$vaccines_are_safe == vaccines_pred2)

acc1
acc2

#konfuzni matice model 1
vaccines_data %>% 
  mutate(pred1 = vaccines_pred1) %>%
  tabyl(vaccines_are_safe,pred1) %>%
  adorn_percentages("row")%>%
  adorn_pct_formatting(digits = 2)%>%
  adorn_ns()

#konfuzni matice model 2
vaccines_data %>% 
  mutate(pred2 = vaccines_pred2) %>%
  tabyl(vaccines_are_safe,pred2) %>%
  adorn_percentages("row")%>%
  adorn_pct_formatting(digits = 2)%>%
  adorn_ns()



