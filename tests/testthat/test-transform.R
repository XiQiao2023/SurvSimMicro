source(file.path("..", "..", "R", "utils.R"))
source(file.path("..", "..", "R", "transform.R"))

# Shared small fixture: 4 subjects x 5 taxa, mix of zeros and positive counts.
make_counts <- function() {
  matrix(
    c(0, 1, 2, 3, 4,
      5, 0, 0, 10, 2,
      0, 0, 0,  0, 0,   # all-zero row
      1, 1, 1,  1, 1),  # uniform row
    nrow = 4, byrow = TRUE,
    dimnames = list(
      paste0("S", 1:4),
      paste0("T", 1:5)
    )
  )
}

test_that("output dimensions and dimnames match input", {
  m   <- make_counts()
  out <- sim_transform(m, method = "clr")

  expect_equal(dim(out), dim(m))
  expect_equal(rownames(out), rownames(m))
  expect_equal(colnames(out), colnames(m))
})

test_that("every CLR row sums to ~0", {
  m    <- make_counts()
  out  <- sim_transform(m, method = "clr")
  rsums <- rowSums(out)

  expect_true(all(abs(rsums) < 1e-8),
              info = paste("Row sums:", paste(rsums, collapse = ", ")))
})

test_that("all-zero row produces all-zero CLR row", {
  m   <- make_counts()
  out <- sim_transform(m, method = "clr")

  # Row 3 is all zeros; after pseudocount it is uniform, so CLR = 0 everywhere.
  expect_true(all(out[3, ] == 0),
              info = paste("All-zero row CLR:", paste(out[3, ], collapse = ", ")))
})

test_that("CLR matches known analytical values", {
  # For a 1x3 matrix with counts [1, 2, 4] and pseudocount 0.5:
  #   adjusted: [1.5, 2.5, 4.5]
  #   log:      log(c(1.5, 2.5, 4.5))
  #   CLR:      log(c) - mean(log(c))
  counts <- matrix(c(1, 2, 4), nrow = 1,
                   dimnames = list("S1", paste0("T", 1:3)))
  adj   <- c(1.5, 2.5, 4.5)
  expected <- log(adj) - mean(log(adj))

  out <- sim_transform(counts, method = "clr", pseudocount = 0.5)
  expect_equal(as.numeric(out[1, ]), expected, tolerance = 1e-12)
})

test_that("CLR is scale-invariant for strictly positive data with negligible pseudocount", {
  # Scale invariance holds exactly when pseudocount = 0, i.e., for strictly
  # positive compositions: CLR(c*x) = CLR(x).  With a pseudocount the property
  # is only approximate; we verify it becomes negligible when counts >> pseudocount.
  set.seed(1)
  counts <- matrix(rpois(20, lambda = 1000) + 1L, nrow = 4)
  counts_scaled        <- counts
  counts_scaled[1, ]   <- counts[1, ] * 7L

  out        <- sim_transform(counts,        method = "clr", pseudocount = 0.5)
  out_scaled <- sim_transform(counts_scaled, method = "clr", pseudocount = 0.5)

  # With counts ~ 1000 and pseudocount 0.5, relative error < 0.1 %
  expect_equal(out[1, ], out_scaled[1, ], tolerance = 1e-3)
})

test_that("CLR with pseudocount = 0 errors on zero counts", {
  m <- make_counts()   # contains zeros
  expect_error(sim_transform(m, method = "clr", pseudocount = 0),
               regexp = "positive")
})

test_that("unimplemented method stops with informative error", {
  m <- make_counts()
  expect_error(sim_transform(m, method = "alr"),
               regexp = "not yet implemented")
  expect_error(sim_transform(m, method = "ilr"),
               regexp = "not yet implemented")
  expect_error(sim_transform(m, method = "magic"),
               regexp = "not yet implemented")
})

test_that("data frame input is accepted and coerced", {
  df  <- as.data.frame(make_counts())
  out <- sim_transform(df, method = "clr")
  expect_true(is.matrix(out))
  expect_equal(nrow(out), nrow(df))
})

test_that("non-numeric input errors", {
  m        <- make_counts()
  m_chr    <- matrix(as.character(m), nrow = nrow(m))
  expect_error(sim_transform(m_chr, method = "clr"),
               regexp = "numeric")
})

test_that("negative counts error", {
  m      <- make_counts()
  m[1,1] <- -1
  expect_error(sim_transform(m, method = "clr"),
               regexp = "negative")
})
