# example_run.R
# End-to-end demonstration of the TTE simulation toolkit on the real microbiome
# count matrix.  Run from the TTE_sims/ project root, e.g.:
#   Rscript examples/example_run.R

# ---- source all modules -------------------------------------------------- #
setwd("~/Library/CloudStorage/OneDrive-UniversityofUtah/TTE_sims")
for (f in list.files("R", full.names = TRUE, pattern = "\\.R$")) source(f)

# ---- load real count matrix ---------------------------------------------- #
load("example_data.Rdata")
# example_data: 144 subjects x 1731 taxa (data frame of raw counts)
cat(sprintf("Count matrix: %d subjects x %d taxa\n",
            nrow(example_data), ncol(example_data)))

# ============================================================
# SCENARIO 1 — 1-stage Weibull-PH
# ============================================================
# Use the first 50 CLR-transformed taxa as predictors.
# Betas are all positive (+0.15) so higher CLR abundance → higher hazard →
# shorter survival under PH.
cat("\n=== Scenario 1: 1-stage Weibull PH ===\n")
set.seed(42)
betas_1 <- rep(0.15, 50)

res1 <- simulate_tte(
  counts          = example_data,
  design          = "one_stage",
  transform_args  = list(method = "clr", pseudocount = 0.5),
  predictor_args  = list(index = 1:50, betas = betas_1),
  tte_args        = list(model         = "ph",
                         dist          = "weibull",
                         params        = list(shape = 1.5),
                         target_median = 5),
  censor_args     = list(type = "administrative"),
  target_event_prop = 0.60,
  seed            = 42L
)

print(res1)

# Direction-of-effect check: subjects with eta above the median should have
# shorter observed times under PH (higher hazard → shorter survival).
eta1   <- res1$eta
hi_eta <- eta1 > median(eta1)
med_hi <- median(res1$outcome$time[hi_eta])
med_lo <- median(res1$outcome$time[!hi_eta])
cat(sprintf("  PH direction: median time (high eta) = %.3f\n", med_hi))
cat(sprintf("                median time (low  eta) = %.3f\n", med_lo))
stopifnot("PH: high-eta subjects should have shorter times" = med_hi < med_lo)
cat("  [OK] PH direction confirmed.\n")

# Quick Cox PH sanity check
if (requireNamespace("survival", quietly = TRUE)) {
  df1 <- data.frame(time   = res1$outcome$time,
                    status = res1$outcome$status,
                    eta    = res1$eta)
  fit1 <- survival::coxph(survival::Surv(time, status) ~ eta, data = df1)
  coef1 <- coef(fit1)[["eta"]]
  cat(sprintf("  CoxPH coefficient on eta: %.4f (expected > 0)\n", coef1))
  stopifnot("CoxPH coefficient should be positive under PH" = coef1 > 0)
  cat("  [OK] CoxPH coefficient positive.\n")
}

# ============================================================
# SCENARIO 2 — 2-stage AFT (mixture-cure)
# ============================================================
# Occurrence predictor: taxa 51-70 (drives cure/noncure membership).
# Timing predictor:    taxa 71-90 with NEGATIVE betas so that under AFT
#   (positive eta → longer survival) negative betas → shorter survival
#   among the susceptible (illustrates the sign-flip documentation).
cat("\n=== Scenario 2: 2-stage Weibull AFT (mixture-cure) ===\n")
set.seed(7)
betas_occ    <- rnorm(20, mean = 0, sd = 0.8)
betas_timing <- rep(-0.2, 20)   # negative → shorter AFT survival

res2 <- simulate_tte(
  counts               = example_data,
  design               = "two_stage",
  transform_args       = list(method = "clr", pseudocount = 0.5),
  predictor_occurrence = list(index = 51:70,  betas = betas_occ),
  predictor_timing     = list(index = 71:90,  betas = betas_timing),
  tte_args             = list(model         = "aft",
                              dist          = "weibull",
                              params        = list(shape = 1.2),
                              target_median = 6),
  censor_args          = list(type                = "administrative_plus_random",
                              administrative_time = 15,
                              random              = list(dist = "exponential",
                                                         rate = 0.08)),
  target_event_prop    = 0.55,
  seed                 = 7L
)

print(res2)

# Susceptibility check: π_i should vary and mean ≈ target_event_prop = 0.55
pi_mean <- mean(res2$tte$pi)
cat(sprintf("  Mean susceptibility pi: %.3f (target 0.55)\n", pi_mean))
stopifnot("mean(pi) should be close to 0.55" = abs(pi_mean - 0.55) < 1e-6)

# Cured subjects must have latent_time == Inf and status == 0
n_cured        <- sum(!res2$tte$event_eligible)
latent_cured   <- res2$tte$latent_time[!res2$tte$event_eligible]
status_cured   <- res2$outcome$status[!res2$tte$event_eligible]
stopifnot("Cured latent times must be Inf"  = all(is.infinite(latent_cured)))
stopifnot("Cured subjects must be censored" = all(status_cured == 0L))
cat(sprintf("  Cured subjects: %d (all have Inf latent time and status=0) [OK]\n",
            n_cured))

# AFT direction-of-effect among eligible subjects only.
# AFT convention (built into sim_tte): higher eta_timing → longer survival.
# Our negative betas mean subjects with large CLR values get more negative eta,
# but among whichever subjects happen to have higher eta the survival is longer.
elig      <- res2$tte$event_eligible
eta2_t    <- res2$eta$timing[elig]
time_elig <- res2$outcome$time[elig]
hi2     <- eta2_t > median(eta2_t)
med2_hi <- median(time_elig[hi2])
med2_lo <- median(time_elig[!hi2])
cat(sprintf("  AFT: median time (high eta) = %.3f\n", med2_hi))
cat(sprintf("       median time (low  eta) = %.3f\n", med2_lo))
# AFT fundamental property: higher eta → longer survival time
stopifnot("AFT: high-eta eligible subjects should have longer times" =
            med2_hi > med2_lo)
cat("  [OK] AFT direction confirmed.\n")

# Optional survreg sanity check (among eligible subjects)
if (requireNamespace("survival", quietly = TRUE)) {
  df2 <- data.frame(time   = time_elig,
                    status = res2$outcome$status[elig],
                    eta    = eta2_t)
  fit2 <- survival::survreg(survival::Surv(time, status) ~ eta,
                             data = df2, dist = "weibull")
  coef2 <- coef(fit2)[["eta"]]
  cat(sprintf("  survreg AFT coefficient on eta: %.4f (expected > 0: higher eta → longer time)\n",
              coef2))
  stopifnot("survreg coef should be positive (AFT: higher eta → longer time)" =
              coef2 > 0)
  cat("  [OK] survreg coefficient positive.\n")
}

# ============================================================
# SCENARIO 3 — 1-stage with custom fn (non-linear predictor)
# ============================================================
cat("\n=== Scenario 3: 1-stage exponential PH with custom fn ===\n")
# Use the product of the first and second CLR taxon as the predictor.
res3 <- simulate_tte(
  counts          = example_data,
  design          = "one_stage",
  transform_args  = list(method = "clr"),
  predictor_args  = list(index = 1:2,
                         fn    = function(a, b) 0.3 * a - 0.3 * b),
  tte_args        = list(model         = "ph",
                         dist          = "exponential",
                         params        = list(shape = 1),
                         target_median = 4),
  censor_args     = list(type = "administrative"),
  target_event_prop = 0.65,
  seed            = 99L
)

print(res3)
stopifnot("custom fn result is numeric eta" = is.numeric(res3$eta))
cat("  [OK] Custom fn scenario completed.\n")

cat("\nAll scenarios completed without error.\n")
