calibrate_intercept <- function(lp_without_intercept, target) {
  uniroot(function(b) mean(plogis(b + lp_without_intercept)) - target,
          interval = c(-12, 12))$root
}

simulate_rwe_cohort <- function(n = 10000, treatment_prevalence = .40,
                                age_selection = .35, comorbidity_selection = .55,
                                severity_selection = .80, unmeasured_confounding = 0,
                                positivity_stress = 1, true_odds_ratio = .80,
                                seed = 2026) {
  set.seed(seed)
  d <- data.frame(
    id = seq_len(n), age = pmin(pmax(rnorm(n, 67, 9), 40), 90),
    female = rbinom(n, 1, .52), comorbidity = rpois(n, 1.8),
    severity = rnorm(n), prior_event = rbinom(n, 1, .22),
    healthcare_use = rpois(n, 4), u = rnorm(n)
  )
  z_age <- as.numeric(scale(d$age)); z_comorb <- as.numeric(scale(d$comorbidity))
  measured_lp <- positivity_stress * (age_selection*z_age +
    comorbidity_selection*z_comorb + severity_selection*d$severity +
    .25*d$female + .45*d$prior_event + .15*as.numeric(scale(d$healthcare_use)))
  treatment_lp <- measured_lp + unmeasured_confounding*d$u
  intercept <- calibrate_intercept(treatment_lp, treatment_prevalence)
  d$ps_true <- plogis(intercept + treatment_lp)
  d$treatment <- rbinom(n, 1, d$ps_true)
  outcome_base <- -2.35 + .025*(d$age-67) + .16*d$comorbidity +
    .48*d$severity + .55*d$prior_event + .04*d$healthcare_use +
    unmeasured_confounding*d$u
  log_or <- log(true_odds_ratio)
  d$p0 <- plogis(outcome_base)
  d$p1 <- plogis(outcome_base + log_or)
  d$outcome <- rbinom(n, 1, ifelse(d$treatment==1,d$p1,d$p0))
  attr(d,"true_ate") <- mean(d$p1-d$p0)
  d
}
