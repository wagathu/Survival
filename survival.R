
# Loading packages --------------------------------------------------------

library(pacman)
p_load(
  tidyverse,
  survival,
  mice,
  GGally,
  ggfortify,
  survminer,
  highcharter,
  ggsurvfit,
  tidycmprsk,
  gtsummary,
  lubridate,
  condSURV
)

# Importing data ----------------------------------------------------------

setwd("~/GitHub")
data <- haven::read_dta("dhs_data2/KEBR8AFL_2022.DTA")

# Data with required columns only -----------------------------------------

# Function for extracting the required columns only
clean_data <- function(data) {
  data <- data |>
    dplyr::select(
      caseid,
      # unique identifier
      v005,
      # Women's individual sample weight
      v006,
      # Month of Interview
      v007,
      # Year of Interview
      v008,
      # Date of interview (CMC)
      v008a,
      # Date of Interview (CDC)
      v012,
      # Respondent's current age
      v013,
      # Age in 5-year groups
      v101,
      # De facto region of residence. Region in which the respondent was interviewed.
      v130,
      # Religion
      v131,
      # Ethnicity
      v139,
      # De jure region of residence
      b0,
      # Child is twin
      b1,
      # Month of birth of child
      b2,
      # Year of birth of child
      b3,
      # Century month code for the date of birth of the child
      b4,
      # Sex of child
      b5,
      # Whether child was alive or dead at the time of interview.
      b6,
      # Age at death
      b7,
      # Age at death (months, imputed)
      b8,
      # Current age of child (if not dead)
      b9,
      # Child lives with whom,
      b11 ,
      # birth interval
      b20,
      # Duration of pregnancy in months
      m15,
      # Place of delivery
      m17,
      # Delivery by cesarean section
      b18,
      #Century Day Code of birth (CDC)
      m19,
      # Birth weight in kilograms (3 decimals)
      m19a,
      # Weight at birth/recall
      v106,
      # mother's education
      m4,
      # The duration of breastfeeding of the child in months.
      m5,
      # duration child breastfed
      m14,
      # Number of antenatal visits during pregnancy
      m70,
      # Whether child's health was checked after discharge
      m15,
      # Place of delivery of the child.
      m18        # Size of child as reported subjectively by the respondent
    )
  
  return(data)
}

df <- clean_data(data)
#haven::as_factor()

# 5 years before the survey -----------------------------------------------

df <- df |>
  mutate(duration = v007 - b2) |>
  filter(duration <= 5)


# Calculating Survival time in days ---------------------------------------

cdc_function <- function(cdc) {
  # Convert CDC to a Date object
  cdc_date <- as.Date("1899-12-31") + cdc
  
  # Calculate the year
  year <- year(cdc_date)
  
  # Step 1: Calculate the number of days since 1900 to yyyy-1
  cumulative_days_1 <- sum(sapply(1900:(year - 1), function(y) 365 + as.integer((y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0)))
  
  # Step 2: Calculate the number of days from yyyy to mm/dd
  cumulative_days_2 <- 0
  xdays <- 0
  
  for (i in 1:(month(cdc_date) - 1)) {
    if (i %in% c(1, 3, 5, 7, 8, 10)) {
      xdays <- 31
    } else if (i == 2) {
      xdays <- 28 + as.integer((year %% 4 == 0 & year %% 100 != 0) | year %% 400 == 0)
    } else if (i %in% c(4, 6, 9, 11)) {
      xdays <- 30
    }
    
    cumulative_days_2 <- cumulative_days_2 + xdays
  }
  
  # Calculate CDC
  cdc_result <- cumulative_days_1 + cumulative_days_2 + (day(cdc_date) - 1)
  
  # Convert the CDC result to a Date object
  result_date <- as.Date("1899-12-31") + cdc_result
  
  return(result_date)
}

# Test the function with the provided example
# cdc_example <- 43302  # CDC for 21 of July, 2018
# result_date <- cdc_function(cdc_example)
# print(result_date)


# Selecting Infants -------------------------------------------------------

df <- df |>
  mutate(surv_time2  = v008a - b18) |> 
  mutate(
    surv_time = case_when(
      substr(b6, 1, 1) == "1" ~ substr(b6, nchar(b6) - 1, nchar(b6)),
      as.numeric(b7) == 1 ~ "30",
      is.na(b7) ~ "40",
      TRUE ~ as.character(as.numeric(b7) * 30)
    ),
    surv_time = as.numeric(surv_time)
  ) |>
  mutate(status = ifelse(b5 == 1, 0, 1)) # 1-dead, 0-censored

# Note: the Surv() function in the {survival} package
# accepts by default TRUE/FALSE, where TRUE is event
# and FALSE is censored; 1/0 where 1 is event and
# 0 is censored; or 2/1 where 2 is event and 1 is censored.
# Please take care to ensure the event indicator is properly formatted.

# Survival Time --------------------------------------------------------------------

df |> 
    ggplot() +
    geom_density(aes(x = surv_time))

  df_time <- df |>
  #filter(surv_time <= 30) |>
  as_factor()

## Density Plot
df_time |>
  filter(status == 1) |>
  #dplyr::filter(!is.na(b7)) |>
  ggplot2::ggplot(aes(x = as.numeric(surv_time2))) +
  ggplot2::geom_density(linewidth = .78) +
  theme_bw() +
  xlab("time")

## Box-Plot
df_time |>
  ggplot(aes(y = surv_time)) +
  geom_boxplot()


# Hazard Rate -------------------------------------------------------------
# 1. Kaplan Mier
#surv_obj <- Surv(df$surv_time, df$status)
km_fit <- survfit(Surv(surv_time, status == 1) ~ 1, df_time)
summary(km_fit)

survfit2(Surv(surv_time, status == 1) ~ 1, data = df_time) %>%
  ggsurvfit() +
  labs(x = "Days",
       y = "Overall survival probability")

autoplot(km_fit) +
  theme_classic()

cumhaz = km_fit$cumhaz
ggsurvplot(km_fit)

# Cumulative hazard rate --------------------------------------------------

CR <- data.frame(
  time = km_fit$time,
  n.risk = km_fit$n.risk,
  n.event = km_fit$n.event,
  survival = km_fit$surv,
  cumalative.hazard = km_fit$cumhaz,
  upper = km_fit$upper,
  lower = km_fit$lower
) |>
  mutate(harz = c(0, diff(cumalative.hazard)))

CR |>
  ggplot(aes(x = factor(time), y = survival)) +
  geom_col(aes(x = factor(time), y = survival), fill = "steelblue") +
  labs(title = "Survival probabilities", y = "S(X)", x = "Time in months") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme_classic() +
  coord_flip()



