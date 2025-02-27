here::i_am("gc2/victor.R")

library(conflicted)
library(here)
library(tidyverse)
conflicts_prefer(dplyr::filter, dplyr::lag)

data <- read_csv(here("data/annual_bat_counts.csv"))

data %>% count(AcousticSite, NightOf) %>% filter(n > 1)

data %>%
    count(AcousticSite, name = "num_nights") %>%
    count(num_nights, name = "num_sites")

n_distinct(data$Year)

data %>%
    group_by(AcousticSite) %>%
    summarize(years = str_c(Year, collapse = ",")) %>%
    count(years, sort = TRUE, name = "num_sites")

data %>% count(Year, sort = TRUE, name = "num_sites")

data2 <- data %>%
    select(
        AcousticSite, Year, dist_to_water,
        ALL, ends_with("_count"),
        starts_with("nights_")
    ) %>%
    pivot_longer(
        starts_with("nights_"),
        names_to = "num_nights_year", names_prefix = "nights_",
        values_to = "num_nights"
    ) %>%
    filter(Year == num_nights_year) %>%
    select(!num_nights_year)

mod1 <- glm(ALL ~ dist_to_water, poisson, data2)
summary(mod1)

tibble(fitted_val = fitted(mod1), squared_resid = resid(mod1)^2) %>%
    ggplot(aes(fitted_val, squared_resid)) +
    geom_point(alpha = 0.3) +
    geom_abline(linetype = "dashed", slope = 1, intercept = 0) +
    labs(x = "Fitted Value", y = "Squared Residual") +
    theme_bw()

tibble(sqrt_fitted_val = sqrt(fitted(mod1)), abs_resid = abs(resid(mod1))) %>%
    ggplot(aes(sqrt_fitted_val, abs_resid)) +
    geom_point(alpha = 0.3) +
    geom_abline(linetype = "dashed", slope = 1, intercept = 0) +
    labs(x = expression(sqrt("Fitted Value")), y = "Absolute Residual") +
    theme_bw()

mod2 <- glm(ALL ~ dist_to_water, quasipoisson, data2)
summary(mod2)

tibble(fitted_val = fitted(mod2), squared_resid = resid(mod2)^2) %>%
    ggplot(aes(fitted_val, squared_resid)) +
    geom_point(alpha = 0.3) +
    geom_abline(linetype = "dashed", slope = 1, intercept = 0) +
    labs(x = "Fitted Value", y = "Squared Residual") +
    theme_bw()

tibble(sqrt_fitted_val = sqrt(fitted(mod2)), abs_resid = abs(resid(mod2))) %>%
    ggplot(aes(sqrt_fitted_val, abs_resid)) +
    geom_point(alpha = 0.3) +
    geom_abline(linetype = "dashed", slope = 1, intercept = 0) +
    labs(x = expression(sqrt("Fitted Value")), y = "Absolute Residual") +
    theme_bw()
