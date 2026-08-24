#' End-to-end simulation of time-to-event data from microbiome counts
#'
#' Top-level driver that chains all four simulation modules in order:
#'
#' \enumerate{
#'   \item \strong{Module 1 — Transform} (\code{\link{sim_transform}}): CLR-
#'     transform the raw count matrix.
#'   \item \strong{Module 2 — Predictor} (\code{\link{sim_predictor}}): build
#'     one linear predictor (1-stage) or two (2-stage, one per role).
#'   \item \strong{Module 3 — TTE} (\code{\link{sim_tte}}): draw latent
#'     failure times; in the 2-stage design also draws cure/noncure membership.
#'   \item \strong{Module 4 — Censoring} (\code{\link{sim_censor}}): apply
#'     censoring to produce the observed \code{(time, status)} outcome.
#' }
#'
#' The caller's RNG state is saved before any randomness is consumed and
#' restored on exit, exactly as \code{sim_generate()} does in \code{test.R}.
#' Sequential calls with different seeds therefore produce independent datasets
#' without side-effects on the global RNG.
#'
#' @section Design modes:
#' \describe{
#'   \item{\strong{1-stage} (\code{design = "one_stage"})}{
#'     A single linear predictor \eqn{\eta_{\text{timing}}} (built from
#'     \code{predictor_args}) drives the Weibull PH or AFT failure-time model.
#'     All subjects are eligible for the event; the observed event proportion
#'     is controlled by the censoring step via \code{target_event_prop}.}
#'   \item{\strong{2-stage / mixture-cure} (\code{design = "two_stage"})}{
#'     Two separate predictors are built: \code{predictor_occurrence} for the
#'     logit probability of being susceptible (\eqn{\eta_{\text{occ}}}), and
#'     \code{predictor_timing} for the failure-time latency
#'     (\eqn{\eta_{\text{timing}}}).  The uniroot intercept in
#'     \code{\link{sim_tte}} is solved so that
#'     \eqn{\bar{\pi} = \texttt{target\_event\_prop}} (the expected noncured
#'     fraction).  Observed censoring is then applied on top by
#'     \code{censor_args}.}
#' }
#'
#' @param counts Numeric matrix or data frame (\eqn{n \times p}) of raw
#'   microbiome counts.  Validated and CLR-transformed by
#'   \code{\link{sim_transform}}.
#' @param design Character scalar: \code{"one_stage"} or \code{"two_stage"}.
#' @param transform_args Named list forwarded to \code{\link{sim_transform}}.
#'   Recognised keys: \code{method} (default \code{"clr"}),
#'   \code{pseudocount} (default \code{0.5}).
#' @param predictor_args Named list forwarded to \code{\link{sim_predictor}}
#'   in the 1-stage design.  Must include \code{index} and exactly one of
#'   \code{betas} or \code{fn}.
#' @param predictor_occurrence Named list forwarded to
#'   \code{\link{sim_predictor}} for the occurrence (incidence) predictor in
#'   the 2-stage design.  Must include \code{index} and \code{betas} or
#'   \code{fn}.
#' @param predictor_timing Named list forwarded to \code{\link{sim_predictor}}
#'   for the latency (timing) predictor in the 2-stage design.
#' @param tte_args Named list forwarded to \code{\link{sim_tte}}.  Recognised
#'   keys: \code{model} (\code{"ph"} or \code{"aft"}), \code{dist}
#'   (\code{"weibull"} or \code{"exponential"}), \code{params}
#'   (e.g. \code{list(shape = 1.5)}), \code{target_median}.  Do not pass
#'   \code{design}, \code{eta_timing}, \code{eta_occurrence}, or
#'   \code{target_event_prop} here — those are wired automatically.
#' @param censor_args Named list forwarded to \code{\link{sim_censor}}.
#'   Recognised keys: \code{type}, \code{administrative_time}, \code{random},
#'   \code{target_event_prop}.  In the 1-stage design, if
#'   \code{target_event_prop} is not present in \code{censor_args} it is
#'   injected from the top-level \code{target_event_prop} argument.
#' @param target_event_prop Numeric scalar in \eqn{(0, 1)}.  In the 1-stage
#'   design: forwarded to \code{\link{sim_censor}} to calibrate the
#'   administrative censoring time.  In the 2-stage design: forwarded to
#'   \code{\link{sim_tte}} to set the expected susceptible fraction
#'   \eqn{\bar{\pi}}.
#' @param seed Integer or \code{NULL}.  Random seed.  When non-\code{NULL}
#'   the global RNG is set before any randomness is consumed and restored to
#'   its prior state on exit.
#'
#' @return A named list of class \code{"sim_tte_data"} with components:
#'   \describe{
#'     \item{\code{outcome}}{A \code{data.frame} with columns \code{id}
#'       (integer), \code{time} (numeric), \code{status} (integer:
#'       1 = event, 0 = censored).  This is the analysis-ready dataset.}
#'     \item{\code{clr}}{Numeric matrix (\eqn{n \times p}).  CLR-transformed
#'       count matrix from Module 1.}
#'     \item{\code{eta}}{Linear predictor(s).  A numeric vector of length
#'       \eqn{n} (1-stage) or a named list with elements \code{timing} and
#'       \code{occurrence} (2-stage).}
#'     \item{\code{tte}}{Full output list from \code{\link{sim_tte}}:
#'       \code{latent_time}, \code{event_eligible}, \code{pi},
#'       \code{intercept}, \code{lambda}, \code{shape}.}
#'     \item{\code{censor}}{Full output \code{data.frame} from
#'       \code{\link{sim_censor}}: \code{id}, \code{time}, \code{status},
#'       \code{latent_time}, \code{censor_time}.}
#'     \item{\code{call}}{A list recording all top-level arguments for
#'       reproducibility.}
#'   }
#'
#' @seealso \code{\link{sim_transform}}, \code{\link{sim_predictor}},
#'   \code{\link{sim_tte}}, \code{\link{sim_censor}}
#'
#' @examples
#' \dontrun{
#' load("example_data.Rdata")
#'
#' # 1-stage Weibull PH — first 30 taxa, uniform betas
#' result_1stage <- simulate_tte(
#'   counts          = example_data,
#'   design          = "one_stage",
#'   transform_args  = list(method = "clr"),
#'   predictor_args  = list(index = 1:30, betas = rep(0.1, 30)),
#'   tte_args        = list(model = "ph", dist = "weibull",
#'                          params = list(shape = 1.5), target_median = 5),
#'   censor_args     = list(type = "administrative"),
#'   target_event_prop = 0.6,
#'   seed            = 1L
#' )
#' head(result_1stage$outcome)
#'
#' # 2-stage mixture-cure AFT
#' result_2stage <- simulate_tte(
#'   counts               = example_data,
#'   design               = "two_stage",
#'   transform_args       = list(method = "clr"),
#'   predictor_occurrence = list(index = 1:15, betas = rnorm(15)),
#'   predictor_timing     = list(index = 16:30, betas = rnorm(15)),
#'   tte_args             = list(model = "aft", dist = "weibull",
#'                               params = list(shape = 1.2), target_median = 4),
#'   censor_args          = list(type = "administrative",
#'                               administrative_time = 12),
#'   target_event_prop    = 0.55,
#'   seed                 = 2L
#' )
#' mean(result_2stage$outcome$status)
#' }
#'
#' @export
simulate_tte <- function(counts,
                         design               = "one_stage",
                         transform_args       = list(),
                         predictor_args       = list(),
                         predictor_occurrence = list(),
                         predictor_timing     = list(),
                         tte_args             = list(),
                         censor_args          = list(),
                         target_event_prop    = 0.5,
                         seed                 = NULL) {

  if (!design %in% c("one_stage", "two_stage"))
    stop("design '", design, "' not recognised; use 'one_stage' or 'two_stage'.")

  # ---- RNG management (restore caller's state on exit) ------------------- #
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  # ---- Module 1: transform ----------------------------------------------- #
  clr_mat <- do.call(sim_transform,
                     c(list(counts = counts), transform_args))

  # ---- Module 2: predictor(s) -------------------------------------------- #
  if (design == "one_stage") {
    eta_timing <- do.call(sim_predictor,
                          c(list(X = clr_mat), predictor_args))
    eta_out <- eta_timing
  } else {
    if (length(predictor_occurrence) == 0L)
      stop("'predictor_occurrence' must be supplied for design = 'two_stage'.")
    if (length(predictor_timing) == 0L)
      stop("'predictor_timing' must be supplied for design = 'two_stage'.")
    eta_occurrence <- do.call(sim_predictor,
                              c(list(X = clr_mat), predictor_occurrence))
    eta_timing     <- do.call(sim_predictor,
                              c(list(X = clr_mat), predictor_timing))
    eta_out <- list(occurrence = eta_occurrence, timing = eta_timing)
  }

  # ---- Module 3: latent failure times ------------------------------------ #
  # Map top-level design names to the internal "1stage"/"2stage" convention
  tte_design <- switch(design, one_stage = "1stage", two_stage = "2stage")

  tte_call <- c(
    list(eta_timing = eta_timing,
         design     = tte_design),
    tte_args
  )

  if (design == "two_stage") {
    tte_call$eta_occurrence    <- eta_occurrence
    tte_call$target_event_prop <- target_event_prop
  }

  tte_result <- do.call(sim_tte, tte_call)

  # ---- Module 4: censoring ----------------------------------------------- #
  # In 1-stage, inject target_event_prop into censor_args if not already there
  cens_call <- c(list(latent_time = tte_result$latent_time), censor_args)
  if (design == "one_stage" && is.null(cens_call$target_event_prop)) {
    cens_call$target_event_prop <- target_event_prop
  }

  censor_result <- do.call(sim_censor, cens_call)

  # ---- assemble return object -------------------------------------------- #
  outcome <- data.frame(
    id     = censor_result$id,
    time   = censor_result$time,
    status = censor_result$status
  )

  out <- list(
    outcome = outcome,
    clr     = clr_mat,
    eta     = eta_out,
    tte     = tte_result,
    censor  = censor_result,
    call    = list(
      design               = design,
      transform_args       = transform_args,
      predictor_args       = predictor_args,
      predictor_occurrence = predictor_occurrence,
      predictor_timing     = predictor_timing,
      tte_args             = tte_args,
      censor_args          = censor_args,
      target_event_prop    = target_event_prop,
      seed                 = seed
    )
  )

  class(out) <- "sim_tte_data"
  out
}


#' Print method for sim_tte_data
#'
#' @param x A \code{sim_tte_data} object returned by \code{\link{simulate_tte}}.
#' @param ... Ignored.
#' @export
print.sim_tte_data <- function(x, ...) {
  n       <- nrow(x$outcome)
  p       <- ncol(x$clr)
  n_evt   <- sum(x$outcome$status)
  med_t   <- median(x$outcome$time)
  design  <- x$call$design

  cat("sim_tte_data  [n =", n, "| p =", p, "| design =", design, "]\n")
  cat("  Events    :", n_evt, "/", n,
      sprintf("(%.1f %%)\n", 100 * n_evt / n))
  cat("  Median observed time :", round(med_t, 3), "\n")
  if (!is.null(x$tte$pi)) {
    cat("  Mean susceptibility pi:", round(mean(x$tte$pi), 3), "\n")
  }
  invisible(x)
}
