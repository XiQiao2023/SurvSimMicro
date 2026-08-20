#' Transform a raw count matrix for survival modelling
#'
#' Applies a compositional or count-based transformation to an
#' \eqn{n \times p} matrix of raw microbiome counts, producing a
#' real-valued matrix of the same dimensions suitable for use as covariates
#' in downstream simulation modules.
#'
#' Currently supported transformations:
#' \describe{
#'   \item{\code{"clr"}}{Centred log-ratio (CLR): for each subject (row) add a
#'     pseudocount of \code{pseudocount} to avoid \eqn{\log(0)}, then compute
#'     \eqn{\text{clr}(x_i)_j = \log(x_{ij} + \delta) - \frac{1}{p}\sum_k \log(x_{ik} + \delta)}.
#'     The result has row sums of zero by construction.}
#' }
#' Additional transforms (ALR, ILR, …) will be added in future releases;
#' passing an unrecognised \code{method} value stops with a clear error.
#'
#' @param counts Numeric matrix (\eqn{n \times p}) of raw microbiome counts
#'   (rows = subjects, columns = taxa).  Coerced and validated by
#'   \code{.check_count_matrix}: must be non-negative with at least one row
#'   and column.
#' @param method Character scalar naming the transformation.  Currently
#'   \code{"clr"} only.  Case-sensitive.
#' @param pseudocount Numeric scalar \eqn{> 0} added to every count before
#'   taking logarithms.  Defaults to \code{0.5}.  Ignored by transforms that
#'   do not use logarithms.
#' @param ... Additional arguments passed to method-specific helpers (reserved
#'   for future transforms such as ALR's \code{ref} taxon or ILR's basis
#'   matrix).
#'
#' @return A numeric matrix of the same dimensions as \code{counts} containing
#'   the transformed values.  Row and column names are preserved.
#'
#' @seealso \code{\link{simulate_tte}}, \code{\link{.check_count_matrix}}
#'
#' @examples
#' \dontrun{
#' load("example_data.Rdata")
#' clr_mat <- sim_transform(example_data, method = "clr")
#' }
#'
#' @export
sim_transform <- function(counts,
                          method      = "clr",
                          pseudocount = 0.5,
                          ...) {
  m <- .check_count_matrix(counts)

  out <- switch(
    method,
    clr = .clr(m, pseudocount = pseudocount),
    alr = stop("method 'alr' not yet implemented"),
    ilr = stop("method 'ilr' not yet implemented"),
    stop("method '", method, "' not yet implemented")
  )

  out
}


# Apply pseudocount zero-handling: add delta to every element.
# Isolated here so alternative strategies (e.g. multiplicative replacement,
# GBM imputation) can be swapped in without touching .clr().
.add_pseudocount <- function(m, pseudocount) {
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L || pseudocount <= 0)
    stop("'pseudocount' must be a single positive number.")
  m + pseudocount
}


# Centred log-ratio transform, computed row-wise.
# Row sums of the result are zero by construction:
#   sum_j [log(x_j + d) - mean_k log(x_k + d)] = 0.
.clr <- function(m, pseudocount) {
  m_adj  <- .add_pseudocount(m, pseudocount)
  log_m  <- log(m_adj)

  if (any(!is.finite(log_m)))
    stop("Non-finite values produced during CLR log step; ",
         "check for zero pseudocount or non-positive counts.")

  # Subtract row geometric mean (= row mean of logs)
  row_means <- rowMeans(log_m)
  out       <- log_m - row_means   # recycling sweeps each row

  dimnames(out) <- dimnames(m)
  out
}
