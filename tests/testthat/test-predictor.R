source(file.path("..", "..", "R", "utils.R"))
source(file.path("..", "..", "R", "transform.R"))
source(file.path("..", "..", "R", "predictor.R"))

# Shared fixture: 10 subjects x 6 taxa CLR matrix with named dims.
make_clr <- function() {
  set.seed(42)
  counts <- matrix(rpois(60, lambda = 10) + 1L, nrow = 10,
                   dimnames = list(paste0("S", 1:10), paste0("T", 1:6)))
  sim_transform(counts, method = "clr")
}

# ---- linear path --------------------------------------------------------- #

test_that("linear path matches manual matrix multiply", {
  X     <- make_clr()
  idx   <- c(1L, 3L, 5L)
  betas <- c(0.5, -1.2, 2.0)

  eta_fn  <- sim_predictor(X, index = idx, betas = betas)
  eta_man <- as.numeric(X[, idx, drop = FALSE] %*% betas)

  expect_equal(eta_fn, eta_man)
})

test_that("linear path returns a plain numeric vector of length n", {
  X   <- make_clr()
  eta <- sim_predictor(X, index = 1:3, betas = c(1, 0, -1))

  expect_true(is.numeric(eta))
  expect_equal(length(eta), nrow(X))
  expect_null(names(eta))
})

test_that("betas length mismatch errors", {
  X <- make_clr()
  expect_error(
    sim_predictor(X, index = 1:3, betas = c(1, 2)),
    regexp = "length\\(betas\\)"
  )
})

# ---- mutual-exclusion checks -------------------------------------------- #

test_that("supplying both betas and fn errors", {
  X <- make_clr()
  expect_error(
    sim_predictor(X, index = 1:2, betas = c(1, 1),
                  fn = function(a, b) a + b),
    regexp = "exactly one"
  )
})

test_that("supplying neither betas nor fn errors", {
  X <- make_clr()
  expect_error(
    sim_predictor(X, index = 1:2),
    regexp = "exactly one"
  )
})

# ---- custom fn path ------------------------------------------------------- #

test_that("custom fn path reproduces manual calculation", {
  X   <- make_clr()
  i1  <- 2L
  i2  <- 4L

  eta_fn  <- sim_predictor(X, index = c(i1, i2),
                            fn = function(a, b) a + 2 * b)
  eta_man <- X[, i1] + 2 * X[, i2]

  expect_equal(eta_fn, as.numeric(eta_man))
})

test_that("custom fn path returns plain numeric vector of length n", {
  X   <- make_clr()
  eta <- sim_predictor(X, index = 1:3,
                       fn = function(a, b, c) a * b + c^2)

  expect_true(is.numeric(eta))
  expect_equal(length(eta), nrow(X))
})

test_that("fn returning wrong length errors", {
  X <- make_clr()
  expect_error(
    sim_predictor(X, index = 1:2, fn = function(a, b) sum(a + b)),
    regexp = "length"
  )
})

test_that("fn returning non-numeric errors", {
  X <- make_clr()
  expect_error(
    sim_predictor(X, index = 1:2, fn = function(a, b) as.character(a + b)),
    regexp = "numeric"
  )
})

# ---- index normalisation -------------------------------------------------- #

test_that("character-name indexing gives same result as integer indexing", {
  X      <- make_clr()
  betas  <- c(1.0, -0.5)

  eta_int  <- sim_predictor(X, index = c(2L, 5L),            betas = betas)
  eta_chr  <- sim_predictor(X, index = c("T2", "T5"),        betas = betas)

  expect_equal(eta_int, eta_chr)
})

test_that("out-of-range integer index errors", {
  X <- make_clr()
  expect_error(sim_predictor(X, index = 99L, betas = 1),
               regexp = "out of range")
})

test_that("unknown column name errors", {
  X <- make_clr()
  expect_error(sim_predictor(X, index = "NOPE", betas = 1),
               regexp = "Unknown column name")
})
