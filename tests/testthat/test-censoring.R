source(file.path("..", "..", "R", "utils.R"))
source(file.path("..", "..", "R", "tte.R"))
source(file.path("..", "..", "R", "censoring.R"))

# Shared fixture: 500 finite latent times drawn from Exp(rate=0.2), median≈3.5
make_latent <- function(n = 500, seed = 1) {
  set.seed(seed)
  stats::rexp(n, rate = 0.2)
}

# ---- return structure ----------------------------------------------------- #

test_that("returns a data.frame with the required columns", {
  T_lat <- make_latent()
  out   <- sim_censor(T_lat, type = "administrative", administrative_time = 10)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("id", "time", "status", "latent_time", "censor_time"))
  expect_equal(nrow(out), length(T_lat))
  expect_equal(out$id, seq_len(length(T_lat)))
})

# ---- type = "administrative" ---------------------------------------------- #

test_that("administrative: everyone past cutoff is censored at cutoff", {
  T_lat  <- make_latent()
  cutoff <- 3
  out    <- sim_censor(T_lat, type = "administrative",
                       administrative_time = cutoff)

  # All censor times equal the cutoff
  expect_true(all(out$censor_time == cutoff))
  # Subjects with latent_time > cutoff must be censored and have time == cutoff
  beyond <- T_lat > cutoff
  expect_true(all(out$status[beyond]  == 0L))
  expect_true(all(out$time[beyond]    == cutoff))
  # Subjects with latent_time <= cutoff must be events
  within <- T_lat <= cutoff
  expect_true(all(out$status[within]  == 1L))
  expect_equal(out$time[within], T_lat[within])
})

test_that("administrative: time == pmin(latent_time, admin_time)", {
  T_lat <- make_latent()
  cutoff <- 5
  out   <- sim_censor(T_lat, type = "administrative",
                      administrative_time = cutoff)

  expect_equal(out$time, pmin(T_lat, cutoff))
})

# ---- type = "none" -------------------------------------------------------- #

test_that("none: all finite latent times are events", {
  T_lat <- make_latent()
  out   <- sim_censor(T_lat, type = "none")

  expect_true(all(out$status == 1L))
})

test_that("none: censor time is max(finite) + 10", {
  T_lat    <- make_latent()
  expected <- max(T_lat) + 10
  out      <- sim_censor(T_lat, type = "none")

  expect_true(all(out$censor_time == expected))
})

# ---- type = "random" ------------------------------------------------------ #

test_that("random exponential: censor times are positive", {
  set.seed(42)
  T_lat <- make_latent()
  out   <- sim_censor(T_lat, type = "random",
                      random = list(dist = "exponential", rate = 0.3))

  expect_true(all(out$censor_time > 0))
})

test_that("random exponential with admin cap: censor times <= admin_time", {
  set.seed(42)
  T_lat  <- make_latent()
  cap    <- 6
  out    <- sim_censor(T_lat, type = "random",
                       administrative_time = cap,
                       random = list(dist = "exponential", rate = 0.3))

  expect_true(all(out$censor_time <= cap))
})

test_that("random exponential without cap: times can exceed any fixed value", {
  set.seed(7)
  T_lat <- stats::rexp(2000, rate = 0.01)   # very long latent times
  out   <- sim_censor(T_lat, type = "random",
                      random = list(dist = "exponential", rate = 0.001))

  # With a very slow censoring rate at least some censor times should be large
  expect_true(max(out$censor_time) > 100)
})

# ---- type = "administrative_plus_random" ---------------------------------- #

test_that("admin_plus_random: censor times <= administrative_time", {
  set.seed(11)
  T_lat <- make_latent()
  cap   <- 4
  out   <- sim_censor(T_lat, type = "administrative_plus_random",
                      administrative_time = cap,
                      random = list(dist = "exponential", rate = 0.5))

  expect_true(all(out$censor_time <= cap))
})

# ---- 2-stage handoff: Inf latent times always censored ------------------- #

test_that("Inf latent times always produce status = 0", {
  set.seed(20)
  n     <- 200
  T_lat <- stats::rexp(n, rate = 0.2)
  # Force half of subjects to be 'cured'
  T_lat[sample(n, n %/% 2)] <- Inf

  out <- sim_censor(T_lat, type = "administrative", administrative_time = 20)

  inf_rows <- is.infinite(T_lat)
  expect_true(all(out$status[inf_rows]  == 0L))
  expect_true(all(out$time[inf_rows]    == 20))
})

test_that("Inf latent times: time equals censor_time", {
  set.seed(21)
  n     <- 100
  T_lat <- c(stats::rexp(50, rate = 0.2), rep(Inf, 50))
  out   <- sim_censor(T_lat, type = "administrative", administrative_time = 10)

  inf_rows <- is.infinite(T_lat)
  expect_equal(out$time[inf_rows], out$censor_time[inf_rows])
})

# ---- target_event_prop calibration --------------------------------------- #

test_that("target_event_prop calibration: realised event rate ≈ target", {
  set.seed(30)
  n     <- 5000
  T_lat <- stats::rexp(n, rate = 0.2)
  target <- 0.55
  out    <- sim_censor(T_lat, type = "administrative",
                       target_event_prop = target)

  realized <- mean(out$status)
  expect_equal(realized, target, tolerance = 0.02)
})

test_that("target_event_prop overrides supplied administrative_time", {
  T_lat  <- make_latent()
  target <- 0.40
  # Supply a deliberately wrong admin_time; calibration should override it
  out    <- sim_censor(T_lat, type = "administrative",
                       administrative_time = 999,
                       target_event_prop   = target)

  realized <- mean(out$status)
  expect_equal(realized, target, tolerance = 0.02)
})

test_that("target_event_prop with type='random' is ignored (no admin cutoff to set)", {
  # Calibration via quantile only works for type="administrative" where all
  # censoring is at a deterministic cutoff.  For type="random", passing
  # target_event_prop has no effect (no admin cutoff to solve); the call must
  # at minimum succeed and return a valid data frame.
  set.seed(31)
  T_lat <- stats::rexp(500, rate = 0.2)
  out   <- sim_censor(T_lat, type = "random",
                      random            = list(dist = "exponential", rate = 0.3),
                      target_event_prop = 0.70)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 500L)
})

# ---- error handling ------------------------------------------------------ #

test_that("unknown type stops with informative message", {
  expect_error(
    sim_censor(make_latent(), type = "interval"),
    regexp = "not yet implemented"
  )
})

test_that("unknown random$dist stops with informative message", {
  expect_error(
    sim_censor(make_latent(), type = "random",
               random = list(dist = "weibull", rate = 0.1)),
    regexp = "not yet implemented"
  )
})

test_that("administrative without administrative_time errors", {
  expect_error(
    sim_censor(make_latent(), type = "administrative"),
    regexp = "administrative_time"
  )
})

test_that("administrative_plus_random without administrative_time errors", {
  expect_error(
    sim_censor(make_latent(), type = "administrative_plus_random",
               random = list(dist = "exponential", rate = 0.2)),
    regexp = "administrative_time"
  )
})

test_that("random without rate errors", {
  expect_error(
    sim_censor(make_latent(), type = "random",
               random = list(dist = "exponential")),
    regexp = "rate"
  )
})
