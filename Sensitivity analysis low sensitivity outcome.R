library(ggplot2)
library(magrittr)
library(purrr)
library(survival)
set.seed(108)

#write a function to simulate data in each trial with 13500 patients
simulate_trial <- function(event_rate_HCTZ, sensitivity) {
  #HR of 0.825
  event_rate_CTD <- event_rate_HCTZ * 0.825
  #block randomize with 13500 patients and block sizes of 6
  treatment <- blockrand::blockrand(n = 13500, 
                       num.levels = 2,
                       levels = c("CTD", "HCTZ"),
                       block.sizes = 3)$treatment
  #Give hazard rate for each treatment
  event_rate <- ifelse(treatment == "HCTZ", event_rate_HCTZ, event_rate_CTD)
  #simulate event times for every patient
  event_time <- purrr::map_dbl(event_rate, ~rexp(1, .))
  #determine if the patient's event was counted
  event_detected <- as.logical(rbinom(13500, 1, sensitivity))
  counted_event_time <- purrr::map2_dbl(event_detected, event_time, ~ifelse(.x, .y, Inf))
  #determine which patients have the first 1055 events with the trial ending at 1055 events
  last_counted_event_time <- sort(counted_event_time)[1055]
  follow_up_time <- pmin(counted_event_time, last_counted_event_time)
  event <- counted_event_time <= last_counted_event_time
  #keep the data needed to make the PH model
  trial_data <- cbind.data.frame(treatment, follow_up_time, event)
  trial_data
}

#write a function that gives the percentage of 1000 trials that show a difference in the primary outcome (aka the power)
simulate_trial_power <- function(event_rate_HCTZ, sensitvity) {
  pvalues <- 1:1000 %>%
    map(~simulate_trial(event_rate_HCTZ, sensitivity)) %>%
    map(~survdiff(Surv(follow_up_time, event) ~ treatment, .)) %>%
    map("pvalue") %>%
    map_lgl(~ 0.05 > .)
  mean(pvalues)
}

#make a vectors of sensitivity and death rate for each level to investigate
event_rate_HCTZ <- 0.03 / 365
sensitivity <- seq(0.5, 1, 0.01)

#determine the power at each death rate and sensitivity
trial_power <- map_dbl(sensitivity, ~simulate_trial_power(event_rate_HCTZ, .))
graph_data <- cbind.data.frame(senstivitity, trial_power)

#graph the results
ggplot(graph_data, aes(x = sensitivity, y = trial_power)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_reverse() +
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
  ylab("Power") +
  xlab("Sensitivity") +
  theme_minimal()

