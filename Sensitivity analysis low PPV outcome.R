library(ggplot2)
library(magrittr)
library(purrr)
library(survival)
set.seed(108)

#write a function to simulate data in each trial with 13500 patients
simulate_trial <- function(death_rate_HCTZ, PPV) {
  #HR of 0.825
  death_rate_CTD <- death_rate_HCTZ * 0.825
  #block randomize with 13500 patients and block sizes of 6
  treatment <- blockrand::blockrand(n = 13500, 
                       num.levels = 2,
                       levels = c("CTD", "HCTZ"),
                       block.sizes = 3)$treatment
  #Give hazard rate for each treatment
  death_rate <- ifelse(treatment == "HCTZ", death_rate_HCTZ, death_rate_CTD)
  #simulate true event times for every patient
  true_event_time <- purrr::map_dbl(death_rate, ~rexp(1, .))
  #simulate false positive event times for every patient
  false_event_rate <- death_rate_HCTZ / PPV - death_rate_HCTZ
  false_event_time <- rexp(13500, false_event_rate)
  #determine which patients have the first 1055 events with the trial ending at 1055 events
  first_event_time <- pmin(false_event_time, true_event_time, na.rm = TRUE)
  trial_end_date <- sort(first_event_time)[1055]
  follow_up_time <- pmin(first_event_time, trial_end_date)
  event <- first_event_time <= follow_up_time
  #collect data needed for PH model
  trial_data <- cbind.data.frame(treatment, follow_up_time, event)
  trial_data
}

#write a function that gives the percentage of 1000 trials that show a difference in the primary outcome (aka the power)
simulate_trial_power <- function(death_rate_HCTZ, PPV) {
  pvalues <- 1:1000 %>%
    map(~simulate_trial(death_rate_HCTZ, PPV)) %>%
    map(~survdiff(Surv(follow_up_time, event) ~ treatment, .)) %>%
    map("pvalue") %>%
    map_lgl(~ 0.05 > .)
  mean(pvalues)
}

#make a vectors of PPV and death rate for each level to investigate
death_rate_HCTZ <- 0.03 / 365
PPV <- seq(0.5, 1, 0.01)

#determine the power at each death rate and PPV
trial_power <- map_dbl(PPV, ~simulate_trial_power(death_rate_HCTZ, .))
graph_data <- cbind.data.frame(PPV, trial_power)

#graph the results
ggplot(graph_data, aes(x = PPV, y = trial_power)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_reverse() +
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
  ylab("Power") +
  xlab("Positive predictive value in HCTZ group") +
  theme_minimal()


