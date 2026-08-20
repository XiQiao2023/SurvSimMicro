source(file.path("..", "..", "R", "utils.R"))
source(file.path("..", "..", "R", "tte.R"))

# ---- helpers --------------------------------------------------------------- #

# Stochastic dominance: P(X <= t) >= P(Y <= t) for all t  <=>  X <=_st Y.
# We use a simple check: median(group_high_eta) < median(group_low_eta)
# after splitting on a binary covariate.

# ---- 1-stage: baseline median calibration --------------------------------- #

test_that("1-stage exponential: sample median ≈ target_median at η=0", {
  set.seed(101)
  n   <- 5000
  out <- sim_tte(eta_timing = rep(0, n), design = "1stage",
                 model = "ph", dist = "exponential",
                 params = list(shape = 1), target_median = 8)

  expect_equal(median(out$latent_time), 8, tolerance = 0.05)   # within 5 %
})

test_that("1-stage Weibull shape=1.5: sample median ≈ target_median at η=0", {
  set.seed(202)
  n   <- 5000
  out <- sim_tte(eta_timing = rep(0, n), design = "1stage",
                 model = "ph", dist = "weibull",
                 params = list(shape = 1.5), target_median = 3)

  expect_equal(median(out$latent_time), 3, tolerance = 0.05)
})

# ---- 1-stage structure ----------------------------------------------------- #

test_that("1-stage returns correct list structure", {
  set.seed(1)
  out <- sim_tte(eta_timing = rnorm(50), design = "1stage",
                 model = "ph", dist = "weibull",
                 params = list(shape = 1), target_median = 5)

  expect_true(is.list(out))
  expect_named(out, c("id", "latent_time", "event_eligible",
                      "pi", "intercept", "lambda", "shape"))
  expect_equal(length(out$latent_time), 50)
  expect_true(all(out$event_eligible))
  expect_null(out$pi)
  expect_null(out$intercept)
})

test_that("1-stage lambda is analytically correct", {
  target_median <- 6
  shape         <- 1.5
  expected_lambda <- log(2) / target_median^shape

  out <- sim_tte(eta_timing = rep(0, 10), design = "1stage",
                 model = "ph", dist = "weibull",
                 params = list(shape = shape), target_median = target_median)

  expect_equal(out$lambda, expected_lambda, tolerance = 1e-12)
  expect_equal(out$shape,  shape,           tolerance = 1e-12)
})

# ---- PH sign convention: higher η ⇒ shorter survival ---------------------- #

test_that("PH: higher eta_timing produces stochastically shorter survival", {
  set.seed(303)
  n    <- 2000
  eta  <- c(rep(-1, n / 2), rep(1, n / 2))
  out  <- sim_tte(eta_timing = eta, design = "1stage",
                  model = "ph", dist = "weibull",
                  params = list(shape = 1.5), target_median = 5)

  med_lo <- median(out$latent_time[eta == -1])
  med_hi <- median(out$latent_time[eta ==  1])
  expect_lt(med_hi, med_lo)   # higher η ⇒ shorter time
})

# ---- AFT sign convention: higher η ⇒ longer survival ---------------------- #

test_that("AFT: higher eta_timing produces stochastically longer survival", {
  set.seed(404)
  n    <- 2000
  eta  <- c(rep(-1, n / 2), rep(1, n / 2))
  out  <- sim_tte(eta_timing = eta, design = "1stage",
                  model = "aft", dist = "weibull",
                  params = list(shape = 1.5), target_median = 5)

  med_lo <- median(out$latent_time[eta == -1])
  med_hi <- median(out$latent_time[eta ==  1])
  expect_gt(med_hi, med_lo)   # higher η ⇒ longer time
})

test_that("AFT baseline median unchanged at η=0", {
  set.seed(505)
  n   <- 5000
  out <- sim_tte(eta_timing = rep(0, n), design = "1stage",
                 model = "aft", dist = "weibull",
                 params = list(shape = 1.5), target_median = 4)

  expect_equal(median(out$latent_time), 4, tolerance = 0.05)
})

# ---- 2-stage: event proportion ------------------------------------------- #

test_that("2-stage: fraction event_eligible ≈ target_event_prop", {
  set.seed(606)
  n   <- 5000
  out <- sim_tte(eta_timing      = rnorm(n, sd = 0.5),
                 design          = "2stage",
                 model           = "ph",
                 dist            = "weibull",
                 params          = list(shape = 1),
                 target_median   = 5,
                 eta_occurrence  = rnorm(n),
                 target_event_prop = 0.65)

  expect_equal(mean(out$event_eligible), 0.65, tolerance = 0.03)
})

test_that("2-stage: π_i varies across subjects (not constant)", {
  set.seed(707)
  n   <- 500
  out <- sim_tte(eta_timing      = rep(0, n),
                 design          = "2stage",
                 model           = "ph",
                 dist            = "exponential",
                 params          = list(shape = 1),
                 target_median   = 5,
                 eta_occurrence  = rnorm(n),
                 target_event_prop = 0.5)

  # π should vary (not all identical)
  expect_gt(sd(out$pi), 0.01)
})

test_that("2-stage: cured subjects have latent_time == Inf", {
  set.seed(808)
  n   <- 500
  out <- sim_tte(eta_timing      = rep(0, n),
                 design          = "2stage",
                 model           = "ph",
                 dist            = "exponential",
                 params          = list(shape = 1),
                 target_median   = 5,
                 eta_occurrence  = rnorm(n),
                 target_event_prop = 0.5)

  cured_times     <- out$latent_time[!out$event_eligible]
  eligible_times  <- out$latent_time[ out$event_eligible]

  expect_true(all(is.infinite(cured_times)))
  expect_true(all(is.finite(eligible_times)))
})

test_that("2-stage: intercept shifts mean π to target_event_prop exactly", {
  set.seed(909)
  n   <- 2000
  eta_occ <- rnorm(n, sd = 2)
  out <- sim_tte(eta_timing      = rep(0, n),
                 design          = "2stage",
                 model           = "ph",
                 dist            = "exponential",
                 params          = list(shape = 1),
                 target_median   = 5,
                 eta_occurrence  = eta_occ,
                 target_event_prop = 0.3)

  # The solved intercept should make mean(π_i) very close to 0.3
  expect_equal(mean(out$pi), 0.3, tolerance = 1e-6)
})

# ---- error handling -------------------------------------------------------- #

test_that("unknown dist stops with informative message", {
  expect_error(
    sim_tte(rep(0, 10), design = "1stage", model = "ph",
            dist = "gompertz", params = list(shape = 1), target_median = 5),
    regexp = "not yet implemented"
  )
})

test_that("unknown model stops with informative message", {
  expect_error(
    sim_tte(rep(0, 10), design = "1stage", model = "loglogistic",
            dist = "weibull", params = list(shape = 1), target_median = 5),
    regexp = "not yet implemented"
  )
})

test_that("unknown design stops with informative message", {
  expect_error(
    sim_tte(rep(0, 10), design = "3stage", model = "ph",
            dist = "weibull", params = list(shape = 1), target_median = 5),
    regexp = "not yet implemented"
  )
})

test_that("2-stage without eta_occurrence errors", {
  expect_error(
    sim_tte(rep(0, 10), design = "2stage", model = "ph",
            dist = "exponential", params = list(shape = 1),
            target_median = 5, target_event_prop = 0.5),
    regexp = "eta_occurrence"
  )
})

test_that("2-stage without target_event_prop errors", {
  expect_error(
    sim_tte(rep(0, 10), design = "2stage", model = "ph",
            dist = "exponential", params = list(shape = 1),
            target_median = 5, eta_occurrence = rep(0, 10)),
    regexp = "target_event_prop"
  )
})

test_that("non-positive target_median errors", {
  expect_error(
    sim_tte(rep(0, 10), design = "1stage", model = "ph",
            dist = "weibull", params = list(shape = 1), target_median = -1),
    regexp = "positive"
  )
})

test_that("non-positive shape errors", {
  expect_error(
    sim_tte(rep(0, 10), design = "1stage", model = "ph",
            dist = "weibull", params = list(shape = -1), target_median = 5),
    regexp = "positive"
  )
})
