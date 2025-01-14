library(ggplot2)
library(posterior)
library(bayesplot)
library(brms)
library(ggpubr)
library(tidyr)
library(dplyr)
library(ggcorrplot)
library(gridExtra)
library(dplyr)
library(loo)
library(invgamma)
library(AER)
library(viridis)
library(corrplot)
library(tidyverse)

options(mc.cores = 4)

data("Guns", package = "AER")

Guns$year <- as.numeric(as.character(Guns$year))
guns_yes <- Guns %>% filter(law == "yes")
guns_no <- Guns %>% filter(law == "no")

guns_yes_mean <- guns_yes %>%
  group_by(year) %>%
  summarise(total_violent = mean(violent))

guns_no_mean <- guns_no %>%
  group_by(year) %>%
  summarise(total_violent = mean(violent))

# Average crime incidents through year with Yes law
plot_yes <- ggplot(guns_yes_mean, aes(x = year, y = total_violent)) +
  geom_line(color = "blue") +
  labs(title = "Crime Rate Over Time (Law in effect)", x = "Year", y = "Average Violent Crime") +
  theme_minimal()
plot_yes

# Average crime incidents through year with No law
plot_no <- ggplot(guns_no_mean, aes(x = year, y = total_violent)) +
  geom_line(color = "red") +
  labs(title = "Crime Rate Over Time (No Law in effect)", x = "Year", y = "Average Violent Crime") +
  theme_minimal()
plot_no

guns_yes_mean$law <- "Law in Effect"
guns_no_mean$law <- "No Law in Effect"

combined_data <- rbind(guns_yes_mean, guns_no_mean)

plot_mix <- ggplot(combined_data, aes(x = year, y = total_violent, color = law)) +
  geom_line(size = 1) +
  labs(
    title = "Crime Rate Over Time (Law in Effect vs. No Law in Effect)",
    x = "Year",
    y = "Average Violent Crime"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal()
plot_mix

# Boxplot of Crime Incidents over Yes/No Shall-carry Laws
ggplot(Guns) +
  aes(x = law, y = violent, fill = law) +
  geom_boxplot() +
  labs(x = "Concealed-carry laws", y = "Crime Rate(per100,000)", fill = "Concealed-carry laws")

Guns$law_numeric <- ifelse(Guns$law == "yes", 1, 0)

# Normalized relevant variables
Guns_normalized <- Guns %>%
  select(-state, -year, -law) %>%
  mutate(across(where(is.factor), as.numeric)) %>%
  scale() %>% 
  as.data.frame()

# Correlation matrix
correlation_matrix <- cor(Guns_normalized)

# Correlation matrix
corrplot(correlation_matrix, 
         method = "color",       
         type = "full",          
         col = colorRampPalette(c("lightskyblue", "white", "coral"))(200), 
         addCoef.col = "black",  
         tl.col = "black",       
         tl.srt = 45,         
         number.cex = 0.7)

# Normalized dataset
Guns_normalized <- Guns_normalized %>%
  mutate(
    state = Guns$state,
    year  = Guns$year,
    law   = Guns$law
)

# Pooled model 
gun_pooled_formula <- bf(
  violent ~ 1 + prisoners #+ afam + cauc + male + population + income + density
            + law,
  family = "gaussian",
  center = FALSE
)

(gun_pooled_priors <- c(
  prior(normal(0, 1), class = "b", coef = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "prisoners"),
  # prior(normal(0, 1), class = "b", coef = "afam"),
  # prior(normal(0, 1), class = "b", coef = "cauc"),
  # prior(normal(0, 1), class = "b", coef = "male"),
  # prior(normal(0, 1), class = "b", coef = "population"),
  # prior(normal(0, 1), class = "b", coef = "income"),
  # prior(normal(0, 1), class = "b", coef = "density"),
  prior(normal(0, 1), class = "b", coef = "lawyes"),
  prior(normal(0, 2), class = "sigma", lb= 0)
))

gun_pooled_fit <- brm(
  formula = gun_pooled_formula,
  prior = gun_pooled_priors,
  data = Guns_normalized,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  sample_prior = "yes"
)

gun_pooled_fit <- add_criterion(
  gun_pooled_fit,
  criterion = "loo",
  moment_match = TRUE,
  overwrite = TRUE
)


# Law model
gun_law_formula <- bf(
  violent ~ 1 + law + (1 + law | state) + (1 | year),
  family = "gaussian",
  center = FALSE
)

(gun_law_priors <- c(
  prior(normal(0, 1), class = "b", coef = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "lawyes"),
  prior(normal(0, 2), class = "sigma", lb= 0),
  prior(normal(0, 5), class = "sd", group = "state", coef = "Intercept"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "lawyes"),
  prior(normal(0, 5), class = "sd", group = "year", coef = "Intercept")
))

gun_law_fit <- brm(
  formula = gun_law_formula,
  prior = gun_law_priors,
  data = Guns_normalized,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

gun_law_fit <- add_criterion(
  gun_law_fit,
  criterion = "loo",
  moment_match = TRUE,
  reloo = TRUE,
  overwrite = TRUE
)

# All model
gun_all_formula <- bf(
  violent ~ 1 + prisoners + law + # + afam + cauc + male + population + income + density
    (1 + prisoners + law | state) + (1 | year), # + afam + cauc + male + population + income + density
  family = "gaussian",
  center = FALSE
)

(gun_all_priors <- c(
  prior(normal(0, 1), class = "b", coef = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "prisoners"),
  # prior(normal(0, 1), class = "b", coef = "afam"),
  # prior(normal(0, 1), class = "b", coef = "cauc"),
  # prior(normal(0, 1), class = "b", coef = "male"),
  # prior(normal(0, 1), class = "b", coef = "population"),
  # prior(normal(0, 1), class = "b", coef = "income"),
  # prior(normal(0, 1), class = "b", coef = "density"),
  prior(normal(0, 1), class = "b", coef = "lawyes"),
  prior(normal(0, 2), class = "sigma", lb = 0),
  
  prior(normal(0, 5), class = "sd", group = "year", coef = "Intercept"),
  
  prior(normal(0, 5), class = "sd", group = "state", coef = "Intercept"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "prisoners"),
  # prior(normal(0, 5), class = "sd", group = "state", coef = "afam"),
  # prior(normal(0, 5), class = "sd", group = "state", coef = "cauc"),
  # prior(normal(0, 5), class = "sd", group = "state", coef = "male"),
  # prior(normal(0, 5), class = "sd", group = "state", coef = "population"),
  # prior(normal(0, 5), class = "sd", group = "state", coef = "income"),
  # prior(normal(0, 5), class = "sd", group = "state", coef = "density"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "lawyes")
))

gun_all_fit <- brm(
  formula = gun_all_formula,
  prior = gun_all_priors,
  data = Guns_normalized,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  sample_prior = "yes",
  save_pars = save_pars(all = TRUE)
)

gun_all_fit <- add_criterion(
  gun_all_fit,
  criterion = "loo",
  moment_match = TRUE,
  reloo = TRUE,
  overwrite = TRUE
)


# Load models
gun_law_fit <- readRDS("notebooks/bda2024/gun_law_model.rds")
gun_all_fit <- readRDS("notebooks/bda2024/gun_all_model.rds")


# Convergence diagnostics
mcmc_rhat_hist(rhat(gun_pooled_fit)) + labs(title = "Pooled Model") + theme(plot.title = element_text(hjust = 0.5))
mcmc_rhat_hist(rhat(gun_law_fit)) + labs(title = "Law Model") + theme(plot.title = element_text(hjust = 0.5))
mcmc_rhat_hist(rhat(gun_all_fit)) + labs(title = "All Model") + theme(plot.title = element_text(hjust = 0.5))

mcmc_neff_hist(neff_ratio(gun_pooled_fit)) + labs(title = "Pooled Model") + theme(plot.title = element_text(hjust = 0.5))
mcmc_neff_hist(neff_ratio(gun_law_fit)) + labs(title = "Law Model") + theme(plot.title = element_text(hjust = 0.5))
mcmc_neff_hist(neff_ratio(gun_all_fit)) + labs(title = "All Model") + theme(plot.title = element_text(hjust = 0.5))

rstan::check_hmc_diagnostics(gun_pooled_fit$fit)
rstan::check_hmc_diagnostics(gun_law_fit$fit)
rstan::check_hmc_diagnostics(gun_all_fit$fit)


# Posterior predictive check
pp_check(gun_pooled_fit) + labs(title = "Pooled Model") + theme(plot.title = element_text(hjust = 0.5))

pp_check(gun_law_fit) + labs(title = "Law Model") + theme(plot.title = element_text(hjust = 0.5))

pp_check(gun_all_fit) + labs(title = "All Model") + theme(plot.title = element_text(hjust = 0.5))


# Model comparison
pooled_loo <- loo(gun_pooled_fit)
law_loo <- loo(gun_law_fit)
all_loo <- loo(gun_all_fit)
pooled_loo
law_loo
all_loo
loo_compare(pooled_loo, law_loo, all_loo)

# k diagnostics
plot(
  pooled_loo,
  diagnostic = c("k"),
  label_points = FALSE,
  main = "Pooled-Model diagnostic plot"
)

plot(
  law_loo,
  diagnostic = c("k"),
  label_points = FALSE,
  main = "Law-Model diagnostic plot"
)

plot(
  all_loo,
  diagnostic = c("k"),
  label_points = FALSE,
  main = "All-Model diagnostic plot"
)


# Prior sensitivity analysis
# Pooled model change
gun_pooled_c_formula <- bf(
  violent ~ 1 + prisoners + afam + cauc + male + 
    population + income + density + law,
  family = "gaussian",
  center = FALSE
)

(gun_pooled_c_priors <- c(
  prior(normal(3, 5), class = "b", coef = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "prisoners"),
  prior(normal(0, 1), class = "b", coef = "afam"),
  prior(normal(0, 1), class = "b", coef = "cauc"),
  prior(normal(0, 1), class = "b", coef = "male"),
  prior(normal(0, 1), class = "b", coef = "population"),
  prior(normal(0, 1), class = "b", coef = "income"),
  prior(normal(0, 1), class = "b", coef = "density"),
  prior(normal(0, 1), class = "b", coef = "lawyes"),
  prior(normal(0, 2), class = "sigma", lb= 0)
))

gun_pooled_c_fit <- brm(
  formula = gun_pooled_c_formula,
  prior = gun_pooled_c_priors,
  data = Guns,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  sample_prior = "yes"
)

m_before <- tibble(value = posterior_samples(gun_pooled_fit)$prior_b_Intercept)
m_after <- tibble(value = posterior_samples(gun_pooled_c_fit)$prior_b_Intercept)
m_before$type <- "Original"
m_after$type <- "Modified"
ggplot(rbind(m_before, m_after), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "Pooled Model",
    x = "Prior Mean",
    y = "Density"
  ) +
  theme_minimal()

m_before <- tibble(value = posterior_samples(gun_pooled_fit)$b_Intercept)
m_after <- tibble(value = posterior_samples(gun_pooled_c_fit)$b_Intercept)
m_before$type <- "Original"
m_after$type <- "Modified"
ggplot(rbind(m_before, m_after), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "Pooled Model",
    x = "Posterior",
    y = "Density"
  ) +
  theme_minimal()




# Law model change
gun_law_c_formula <- bf(
  violent ~ 1 + law + (1 + law | state) + (1 | year),
  family = "gaussian",
  center = FALSE
)

get_prior(gun_law_formula, data = Guns)

(gun_law_c_priors <- c(
  prior(normal(3, 5), class = "b", coef = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "lawyes"),
  prior(normal(0, 2), class = "sigma", lb= 0),
  prior(normal(2, 3), class = "sd", group = "state", coef = "Intercept"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "lawyes"),
  prior(normal(0, 5), class = "sd", group = "year", coef = "Intercept")
))

gun_law_c_fit <- brm(
  formula = gun_law_c_formula,
  prior = gun_law_c_priors,
  data = Guns,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  sample_prior = "yes"
)

gun_law_fit_prior <- brm(
  formula = gun_law_formula,
  prior = gun_law_priors,
  data = Guns,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  sample_prior = "only"
)

m_before <- tibble(value = posterior_samples(gun_law_fit_prior)$b_Intercept)
m_after <- tibble(value = posterior_samples(gun_law_c_fit)$prior_b_Intercept)
m_before$type <- "Original"
m_after$type <- "Modified"
ggplot(rbind(m_before, m_after), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "Law Model",
    x = "Prior Mean",
    y = "Density"
  ) +
  theme_minimal()

m_before <- tibble(value = posterior_samples(gun_law_fit_prior)$sd_state__Intercept)
m_after <- tibble(value = posterior_samples(gun_law_c_fit)$prior_sd_state__Intercept)
m_before$type <- "Original"
m_after$type <- "Modified"
ggplot(rbind(m_before, m_after), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "Law Model",
    x = "Prior Sigma",
    y = "Density"
  ) +
  theme_minimal()

m_before <- tibble(value = posterior_samples(gun_law_fit)$sd_state__Intercept)
m_after <- tibble(value = posterior_samples(gun_law_c_fit)$sd_state__Intercept)
m_before$type <- "Original"
m_after$type <- "Modified"
ggplot(rbind(m_before, m_after), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "Law Model",
    x = "Posterior",
    y = "Density"
  ) +
  theme_minimal()




# All model change
gun_all_c_formula <- bf(
  violent ~ 1 + prisoners + afam + cauc + male + population + income + density + law + 
    (1 + prisoners + afam + cauc + male + population + income + density + law | state) + (1 | year),
  family = "gaussian",
  center = FALSE
)

get_prior(gun_all_c_formula, data = Guns)

(gun_all_c_priors <- c(
  prior(normal(3, 5), class = "b", coef = "Intercept"),
  prior(normal(0, 1), class = "b", coef = "prisoners"),
  prior(normal(0, 1), class = "b", coef = "afam"),
  prior(normal(0, 1), class = "b", coef = "cauc"),
  prior(normal(0, 1), class = "b", coef = "male"),
  prior(normal(0, 1), class = "b", coef = "population"),
  prior(normal(0, 1), class = "b", coef = "income"),
  prior(normal(0, 1), class = "b", coef = "density"),
  prior(normal(0, 1), class = "b", coef = "lawyes"),
  prior(normal(0, 2), class = "sigma", lb = 0),
  
  prior(normal(0, 5), class = "sd", group = "year", coef = "Intercept"),
  
  prior(normal(0, 1), class = "sd", group = "state", coef = "Intercept"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "prisoners"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "afam"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "cauc"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "male"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "population"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "income"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "density"),
  prior(normal(0, 5), class = "sd", group = "state", coef = "lawyes")
))

gun_all_c_fit <- brm(
  formula = gun_all_c_formula,
  prior = gun_all_c_priors,
  data = Guns,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  sample_prior = "yes"
)

gun_all_fit_prior <- brm(
  formula = gun_all_formula,
  prior = gun_all_priors,
  data = Guns,
  iter = 2000,
  warmup = 1000,
  chains = 4,
  sample_prior = "only"
)

m_before <- tibble(value = posterior_samples(gun_all_fit_prior)$b_Intercept)
m_after <- tibble(value = posterior_samples(gun_all_c_fit)$prior_b_Intercept)
m_before$type <- "Original"
m_after$type <- "Modified"
ggplot(rbind(m_before, m_after), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "All Model",
    x = "Prior Mean",
    y = "Density"
  ) +
  theme_minimal()

m_before <- tibble(value = posterior_samples(gun_all_fit_prior)$sd_state__Intercept)
m_after <- tibble(value = posterior_samples(gun_all_c_fit)$prior_sd_state__Intercept)
m_before$type <- "Original"
m_after$type <- "Modified"
ggplot(rbind(m_before, m_after), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "All Model",
    x = "Prior Sigma",
    y = "Density"
  ) +
  theme_minimal()

m_before <- tibble(value = posterior_samples(gun_all_fit)$sd_state__Intercept)
m_after <- tibble(value = posterior_samples(gun_all_c_fit)$sd_state__Intercept)
m_before$type <- "Original"
m_after$type <- "Modified"
ggplot(rbind(m_before, m_after), aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  labs(
    title = "All Model",
    x = "Posterior",
    y = "Density"
  ) +
  theme_minimal()


# Save models
saveRDS(gun_law_fit, file = "gun_law_model.rds")
saveRDS(gun_all_fit, file = "gun_all_model.rds")



#---------------------------------------------
# E X P E R I M E N T    S P A C E
#---------------------------------------------

gun_law_test_fit <- brm(
  formula = gun_law_formula,
  prior = gun_law_priors,
  data = Guns,
  iter = 6000,
  warmup = 1000,
  chains = 4,
  control = list(adapt_delta = 0.95)
)



















