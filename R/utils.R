#' Null-coalescing operator
#'
#' Returns \code{x} when \code{x} is not \code{NULL}; otherwise returns
#' \code{y}.  Modelled after the \code{\%||\%} idiom used throughout the
#' simulation toolkit.
#'
#' @param x Any R object.  Returned as-is unless it is \code{NULL}.
#' @param y Any R object.  Returned when \code{x} is \code{NULL}.
#'
#' @return \code{x} if \code{!is.null(x)}, else \code{y}.
#'
#' @examples
#' NULL %||% 42        # 42
#' "hello" %||% "bye" # "hello"
#'
#' @export
`%||%` <- function(x, y) if (!is.null(x)) x else y


#' Validate and coerce a raw count matrix
#'
#' Coerces \code{x} to a numeric matrix and performs three sanity checks
#' required before any downstream simulation step: the matrix must be numeric,
#' must contain no negative values, and must have at least one row and one
#' column.  On success the coerced numeric matrix is returned so callers need
#' not repeat the coercion.
#'
#' @param x An object that can be coerced to a numeric matrix via
#'   \code{as.matrix} — typically a data frame or matrix of raw microbiome
#'   counts.  Rows correspond to subjects and columns to taxa.
#'
#' @return A numeric matrix with the same dimensions and values as \code{x}
#'   (after coercion).
#'
#' @section Errors:
#' Stops with a descriptive message when:
#' \itemize{
#'   \item The resulting matrix is not numeric (e.g., character columns were
#'     present).
#'   \item Any element is strictly negative.
#'   \item The matrix has zero rows or zero columns.
#' }
#'
#' @examples
#' m <- matrix(c(0, 3, 5, 12), nrow = 2)
#' .check_count_matrix(m)  # returns m unchanged
#'
#' @keywords internal
.check_count_matrix <- function(x) {
  m <- as.matrix(x)
  if (!is.numeric(m))
    stop("Count matrix must be numeric; got mode '", mode(m), "'.")
  if (nrow(m) == 0L)
    stop("Count matrix must have at least one row.")
  if (ncol(m) == 0L)
    stop("Count matrix must have at least one column.")
  if (any(m < 0, na.rm = TRUE))
    stop("Count matrix contains negative values; raw counts must be >= 0.")
  m
}


#' Normalise a column index to integer positions
#'
#' Accepts a column selector that is either a vector of integer positions or a
#' character vector of column names and returns a validated integer position
#' vector in the range \code{[1, p]}.  This allows all module functions to
#' accept either style of column reference without duplicating the validation
#' logic.
#'
#' @param index Integer vector of 1-based column positions \emph{or} character
#'   vector of column names to look up in \code{colnames}.  Mixed integer /
#'   character vectors are not supported; the type of the first element
#'   determines dispatch.
#' @param p Integer scalar.  Total number of columns in the matrix.  Integer
#'   positions must lie in \code{[1, p]}.
#' @param colnames Character vector of length \code{p} giving the column names
#'   of the matrix.  Required when \code{index} is a character vector; may be
#'   \code{NULL} when \code{index} is integer.
#'
#' @return Integer vector of validated column positions, each in
#'   \code{[1, p]}.
#'
#' @section Errors:
#' Stops with a descriptive message when:
#' \itemize{
#'   \item \code{index} is character and any name is not found in
#'     \code{colnames}.
#'   \item \code{index} is integer and any position is outside
#'     \code{[1, p]}.
#' }
#'
#' @examples
#' .as_index(c(1L, 3L), p = 5L, colnames = letters[1:5])
#' .as_index(c("a", "c"), p = 5L, colnames = letters[1:5])
#'
#' @keywords internal
.as_index <- function(index, p, colnames) {
  if (is.character(index)) {
    if (is.null(colnames))
      stop("'colnames' must be supplied when 'index' is a character vector.")
    bad <- setdiff(index, colnames)
    if (length(bad) > 0L)
      stop("Unknown column name(s): ",
           paste(bad, collapse = ", "), ".")
    match(index, colnames)
  } else {
    idx <- as.integer(index)
    out_of_range <- idx < 1L | idx > as.integer(p)
    if (any(out_of_range))
      stop("Column position(s) out of range [1, ", p, "]: ",
           paste(idx[out_of_range], collapse = ", "), ".")
    idx
  }
}
