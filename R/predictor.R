#' Construct a linear predictor from a transformed count matrix
#'
#' Selects a subset of columns from a CLR-transformed microbiome matrix and
#' computes a linear predictor \eqn{\eta} of length \eqn{n} via one of two
#' mutually-exclusive paths:
#'
#' \describe{
#'   \item{Linear path (\code{betas} supplied)}{
#'     \eqn{\eta = X_{\mathcal{S}} \beta}, a standard matrix–vector product.
#'     \code{length(betas)} must equal \code{length(index)}.}
#'   \item{Custom path (\code{fn} supplied)}{
#'     Each selected column is passed as a separate positional argument to
#'     \code{fn}, so a user may write e.g.
#'     \code{fn = function(x1, x2, x3) x1 * x2 + sin(x3)}.  The function
#'     must return a numeric vector of length \eqn{n}.}
#' }
#'
#' Exactly one of \code{betas} or \code{fn} must be supplied; providing both
#' or neither is an error.
#'
#' The returned \eqn{\eta} is passed directly to \code{\link{sim_tte}}.  In a
#' \strong{2-stage mixture-cure design} call \code{sim_predictor} twice with
#' different \code{index}/\code{betas}/\code{fn} arguments — once to build
#' \eqn{\eta_{\text{occurrence}}} (logit probability of being susceptible) and
#' once to build \eqn{\eta_{\text{timing}}} (failure-time linear predictor
#' among the susceptible).
#'
#' @param X Numeric matrix (\eqn{n \times p}).  Typically the CLR-transformed
#'   output of \code{\link{sim_transform}}.  Must be numeric with at least one
#'   row and column.
#' @param index Integer vector of column positions \emph{or} character vector
#'   of column names identifying the active predictor set \eqn{\mathcal{S}}.
#'   Validated by \code{\link{.as_index}}.
#' @param betas Numeric vector of length \eqn{|\mathcal{S}|}.  Coefficients
#'   for the linear path \eqn{\eta = X_{\mathcal{S}} \beta}.  Exactly one of
#'   \code{betas} and \code{fn} must be non-\code{NULL}.
#' @param fn A function whose formal arguments correspond, in order, to the
#'   selected columns of \code{X}.  It is called as
#'   \code{do.call(fn, unname(as.list(as.data.frame(X[, index, drop=FALSE]))))}
#'   and must return a numeric vector of length \eqn{n}.  Exactly one of
#'   \code{betas} and \code{fn} must be non-\code{NULL}.
#' @param ... Additional arguments reserved for future use (currently ignored).
#'
#' @return A numeric vector of length \eqn{n} containing the linear predictor
#'   \eqn{\eta_i} for each subject, on the model's link scale.  No names are
#'   attached; the vector is ready to pass to \code{\link{sim_tte}}.
#'
#' @seealso \code{\link{sim_transform}}, \code{\link{sim_tte}},
#'   \code{\link{simulate_tte}}, \code{\link{.as_index}}
#'
#' @examples
#' \dontrun{
#' load("example_data.Rdata")
#' clr <- sim_transform(example_data, method = "clr")
#'
#' # Linear path
#' eta <- sim_predictor(clr, index = 1:10, betas = rep(0.1, 10))
#'
#' # Custom (non-linear) path
#' eta2 <- sim_predictor(clr, index = c(1, 2),
#'                       fn = function(a, b) a + 2 * b)
#'
#' # 2-stage: call twice with different index/betas
#' eta_occ    <- sim_predictor(clr, index = 1:5,   betas = rnorm(5))
#' eta_timing <- sim_predictor(clr, index = 6:15,  betas = rnorm(10))
#' }
#'
#' @export
sim_predictor <- function(X, index, betas = NULL, fn = NULL, ...) {

  # ---- input validation -------------------------------------------------- #
  if (!is.matrix(X) || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
  if (nrow(X) == 0L || ncol(X) == 0L)
    stop("'X' must have at least one row and one column.")

  idx <- .as_index(index, p = ncol(X), colnames = colnames(X))

  has_betas <- !is.null(betas)
  has_fn    <- !is.null(fn)

  if (has_betas && has_fn)
    stop("Supply exactly one of 'betas' or 'fn', not both.")
  if (!has_betas && !has_fn)
    stop("Supply exactly one of 'betas' or 'fn'; neither was provided.")

  n  <- nrow(X)
  Xs <- X[, idx, drop = FALSE]

  # ---- linear path ------------------------------------------------------- #
  if (has_betas) {
    betas <- as.numeric(betas)
    if (length(betas) != length(idx))
      stop("length(betas) (", length(betas), ") must equal length(index) (",
           length(idx), ").")
    eta <- as.numeric(Xs %*% betas)
    return(eta)
  }

  # ---- custom fn path ---------------------------------------------------- #
  if (!is.function(fn))
    stop("'fn' must be a function.")

  eta <- do.call(fn, unname(as.list(as.data.frame(Xs))))

  if (!is.numeric(eta))
    stop("'fn' must return a numeric vector; got '", class(eta)[1L], "'.")
  if (length(eta) != n)
    stop("'fn' returned a vector of length ", length(eta),
         " but X has ", n, " rows.")

  as.numeric(eta)
}
