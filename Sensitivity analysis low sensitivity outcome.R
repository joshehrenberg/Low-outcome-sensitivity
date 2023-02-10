library(ggplot2)
library(magrittr)
set.seed(108)

#write a function to simulate data in each trial with 13500 patients
simulate_trial <- function(death_rate_HCTZ, sensitivity) {
  #HR of 0.825
  death_rate_CTD <- death_rate_HCTZ * 0.825
  #block randomize with 6000 patients and block sizes of 6
  treatment <- blockrand::blockrand(n = 13500, 
                       num.levels = 2,
                       levels = c("CTD", "HCTZ"),
                       block.sizes = 3)$treatment
  #Give hazard rate for each treatment
  death_rate <- vector("numeric", length = 13500)
  for(i in seq_len(13500)) {
    death_rate[i] <- if(treatment[i] == "HCTZ") {
      death_rate[i] <- death_rate_HCTZ
    }
    else {
      death_rate[i] <- death_rate_CTD
    }
  }
  #simulate true event times for every patient
  true_event_time <- purrr::map_dbl(death_rate, ~rexp(1, .))
  #simulate false positive event times for every patient
  false_event_rate <- death_rate_HCTZ / sensitivity - death_rate_HCTZ
  false_event_time <- rexp(13500, false_event_rate)
  #determine which patients have the first 1055 events with the trial ending at 1055 events
  first_event_time <- pmin(false_event_time, true_event_time, na.rm = TRUE)
  trial_end_date <- sort(first_event_time)[1055]
  follow_up_time <- pmin(first_event_time, trial_end_date)
  event <- first_event_time <= follow_up_time
  #determine if the trial detected an effect
  trial_data <- cbind.data.frame(treatment, follow_up_time, event)
  trial_data
}

#write a function that gives the percentage of 1000 trials that show a difference in the primary outcome (aka the power)
simulate_trials <- function(death_rate_HCTZ, sensitivity) {
  death_rates_HCTZ <- rep(death_rate_HCTZ, times = 1000)
  sensitivities <- rep(sensitivity, times = 1000)
  trials <- purrr::map2(death_rates_HCTZ, sensitivities, simulate_trial)
  trials
}

#make a data frame of all sensitivities and death_rates for the analysis
death_rate_HCTZ <- 0.03 / 365
sensitivity <- seq(0.5, 1, 0.01)
power_data <- cbind.data.frame(death_rate_HCTZ, sensitivity)


#determine the power at each death rate and sensitivity
trial_power <- vector("list", length(sensitivity))
for(i in seq_len(length(sensitivity))) {
  trial_power[[i]] <- simulate_trials(death_rate_HCTZ, sensitivity[i]) %>%
    purrr::map(~survival::survdiff(Surv(follow_up_time, event) ~ treatment, .)) %>%
    purrr::map("pvalue") %>%
    purrr::map_lgl(~ 0.05 > .) %>%
    mean()
}
trial_power <- unlist(trial_power)
graph_data <- cbind.data.frame(sensitivity, trial_power)

#graph the results
ggplot(graph_data, aes(x = sensitivity, y = trial_power)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_reverse() +
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1)) +
  ylab("Power") +
  xlab("Sensitivity in HCTZ group") +
  theme_minimal()


