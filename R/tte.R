#' Simulate latent failure times
#'
#' Draws latent (uncensored) failure times for \eqn{n} subjects under a
#' proportional-hazards (PH) or accelerated failure-time (AFT) model, using
#' \pkg{simsurv} as the survival engine.  The baseline scale parameter is
#' back-solved \emph{analytically} from \code{target_median}, so that a subject
#' with \eqn{\eta = 0} has exactly the requested median survival.
#'
#' @section Baseline calibration (analytic, not iterative):
#' The \pkg{simsurv} Weibull parameterisation uses cumulative hazard
#' \eqn{H_0(t) = \lambda t^\gamma}.  The baseline median \eqn{m_0} satisfies
#' \eqn{\lambda m_0^\gamma = \log 2}, so:
#' \deqn{\lambda = \frac{\log 2}{m_0^\gamma}.}
#' \code{dist = "exponential"} is handled as Weibull with \eqn{\gamma = 1},
#' giving \eqn{\lambda = \log(2) / m_0}.  Because the calibration is at
#' \eqn{\eta = 0}, it is unaffected by the choice of \code{model} (PH or AFT).
#'
#' @section Sign conventions (read before setting betas):
#' \describe{
#'   \item{\strong{PH} (\code{model = "ph"}):}{
#'     Person hazard \eqn{h_i(t) = h_0(t) \exp(\eta_i)}.
#'     \strong{Positive \eqn{\eta} \eqn{\Rightarrow} higher hazard
#'     \eqn{\Rightarrow} shorter survival.}
#'     Implemented via \code{betas = c(eta = 1)} in \pkg{simsurv}.}
#'   \item{\strong{AFT} (\code{model = "aft"}):}{
#'     Acceleration factor \eqn{\exp(\eta_i)}.
#'     \strong{Positive \eqn{\eta} \eqn{\Rightarrow} longer survival.}
#'     For the Weibull family (the unique PH \eqn{\cap} AFT distribution) an
#'     AFT model with predictor \eqn{\eta} is algebraically equivalent to a PH
#'     model with predictor \eqn{-\gamma \eta}.  The implementation therefore
#'     reuses the same \pkg{simsurv} PH call with
#'     \code{betas = c(eta = -shape)} (\eqn{-1} for exponential).
#'     \strong{The baseline median is unchanged at \eqn{\eta = 0}}, so the
#'     analytic calibration still holds.}
#' }
#'
#' @section 1-stage vs. 2-stage (mixture-cure) design:
#' \describe{
#'   \item{\strong{1-stage} (\code{design = "1stage"}):}{
#'     Every subject gets a latent failure time drawn from the Weibull PH/AFT
#'     model parameterised by \code{eta_timing}.  All subjects are
#'     \code{event_eligible}.  The observed event proportion is controlled
#'     entirely by the censoring step (\code{\link{sim_censor}});
#'     \code{target_event_prop} is ignored here.}
#'   \item{\strong{2-stage / mixture-cure} (\code{design = "2stage"}):}{
#'     The population is a mixture of a \emph{cured} group that never
#'     experiences the event and a \emph{noncured} (susceptible) group that
#'     eventually fails.  The population survivor function is
#'     \deqn{S(t) = (1 - \pi_i) + \pi_i S_1(t),}
#'     where \eqn{\pi_i} is the person-specific susceptibility probability and
#'     \eqn{S_1(t)} is the noncured survivor function.
#'
#'     \strong{Stage 1 — incidence (cure/noncure membership):}
#'     \eqn{\pi_i = \mathrm{plogis}(a + \eta_{\mathrm{occ},i})}.
#'     The scalar intercept \eqn{a} is solved with \code{uniroot()} so that
#'     \eqn{\bar{\pi} = \texttt{target\_event\_prop}}.  Each subject is then
#'     drawn as \eqn{\mathrm{event\_eligible}_i \sim \mathrm{Bernoulli}(\pi_i)}.
#'
#'     \strong{Stage 2 — latency (failure time among susceptible):}
#'     Latent failure times are drawn for ALL subjects via \pkg{simsurv} using
#'     \code{eta_timing}, then set to \code{Inf} for cured subjects
#'     (\code{event_eligible == 0}).  Censoring is applied in
#'     \code{\link{sim_censor}}, not here.}
#' }
#'
#' @param eta_timing Numeric vector of length \eqn{n}.  Linear predictor for
#'   the failure-time model (latency); output of \code{\link{sim_predictor}}.
#' @param design Character scalar: \code{"1stage"} or \code{"2stage"}.
#' @param model Character scalar: \code{"ph"} (proportional hazards) or
#'   \code{"aft"} (accelerated failure time).  See Sign conventions above.
#' @param dist Character scalar: \code{"weibull"} or \code{"exponential"}.
#'   \code{"exponential"} is Weibull with \eqn{\gamma = 1}.
#' @param params Named list of distribution parameters.  For Weibull supply
#'   \code{list(shape = \eqn{\gamma})} with \eqn{\gamma > 0}.  Defaults to
#'   \code{list(shape = 1)} (exponential).
#' @param target_median Positive numeric scalar.  Baseline (\eqn{\eta = 0})
#'   median survival time.  \eqn{\lambda} is derived analytically; do not
#'   supply \eqn{\lambda} directly.
#' @param eta_occurrence Numeric vector of length \eqn{n} or \code{NULL}.
#'   Logit-scale linear predictor for cure/noncure membership.  Required when
#'   \code{design = "2stage"}.  Output of a separate \code{\link{sim_predictor}}
#'   call using columns/betas chosen to model incidence.
#' @param target_event_prop Numeric scalar in \eqn{(0, 1)} or \code{NULL}.
#'   Target mean susceptibility \eqn{\bar{\pi}}.  Required when
#'   \code{design = "2stage"}.  The uniroot intercept is solved so that
#'   \eqn{\mathrm{mean}(\pi_i) = \texttt{target\_event\_prop}}.
#' @param ... Reserved for future arguments (currently ignored).
#'
#' @return A named list with components:
#'   \describe{
#'     \item{\code{id}}{Integer vector \code{1:n}.}
#'     \item{\code{latent_time}}{Numeric vector of length \eqn{n}.  Latent
#'       failure time per subject; \code{Inf} for cured subjects in the
#'       2-stage design.}
#'     \item{\code{event_eligible}}{Logical vector of length \eqn{n}.
#'       \code{TRUE} for susceptible (noncured) subjects; always \code{TRUE}
#'       in the 1-stage design.}
#'     \item{\code{pi}}{Numeric vector of length \eqn{n} or \code{NULL}.
#'       Person-specific susceptibility probabilities \eqn{\pi_i}; \code{NULL}
#'       in the 1-stage design.}
#'     \item{\code{intercept}}{Numeric scalar or \code{NULL}.  Solved logit
#'       intercept \eqn{a}; \code{NULL} in the 1-stage design.}
#'     \item{\code{lambda}}{Numeric scalar.  Calibrated Weibull rate parameter
#'       \eqn{\lambda = \log(2) / m_0^\gamma}.}
#'     \item{\code{shape}}{Numeric scalar.  Weibull shape \eqn{\gamma} used.}
#'   }
#'
#' @importFrom simsurv simsurv
#' @importFrom stats plogis rbinom uniroot
#'
#' @seealso \code{\link{sim_predictor}}, \code{\link{sim_censor}},
#'   \code{\link{simulate_tte}}, \code{\link[simsurv]{simsurv}}
#'
#' @examples
#' \dontrun{
#' # 1-stage exponential with flat (zero) predictor: median should ≈ 5
#' out <- sim_tte(eta_timing = rep(0, 500), design = "1stage",
#'                model = "ph", dist = "exponential",
#'                params = list(shape = 1), target_median = 5)
#' median(out$latent_time)   # ≈ 5
#'
#' # 2-stage mixture cure
#' set.seed(1)
#' n   <- 1000
#' eta_occ    <- rnorm(n)
#' eta_timing <- rnorm(n, sd = 0.5)
#' out2 <- sim_tte(eta_timing, design = "2stage", model = "ph",
#'                 dist = "weibull", params = list(shape = 1.5),
#'                 target_median = 5,
#'                 eta_occurrence = eta_occ, target_event_prop = 0.7)
#' mean(out2$event_eligible)  # ≈ 0.7
#' }
#'
#' @export
sim_tte <- function(eta_timing,
                    design            = "1stage",
                    model             = "ph",
                    dist              = "weibull",
                    params            = list(shape = 1),
                    target_median,
                    eta_occurrence    = NULL,
                    target_event_prop = NULL,
                    ...) {

  # ---- input validation -------------------------------------------------- #
  if (!is.numeric(eta_timing) || length(eta_timing) < 1L)
    stop("'eta_timing' must be a non-empty numeric vector.")
  n <- length(eta_timing)

  if (!is.numeric(target_median) || length(target_median) != 1L ||
      target_median <= 0)
    stop("'target_median' must be a single positive number.")

  # ---- resolve shape / dist ---------------------------------------------- #
  shape <- switch(
    dist,
    weibull     = {
      g <- params$shape %||% 1
      if (!is.numeric(g) || length(g) != 1L || g <= 0)
        stop("'params$shape' must be a single positive number.")
      g
    },
    exponential = 1,
    stop("dist '", dist, "' not yet implemented")
  )

  if (!model %in% c("ph", "aft"))
    stop("model '", model, "' not yet implemented")

  if (!design %in% c("1stage", "2stage"))
    stop("design '", design, "' not yet implemented")

  # ---- analytic scale calibration ---------------------------------------- #
  # λ·m0^γ = log(2)  ⟹  λ = log(2) / target_median^γ
  lambda <- log(2) / (target_median ^ shape)

  # ---- simsurv beta coefficient per model -------------------------------- #
  # PH : betas = c(eta = 1)         positive η ⇒ higher hazard ⇒ shorter time
  # AFT: betas = c(eta = -shape)    positive η ⇒ longer time (sign flip!)
  beta_eta <- switch(model,
    ph  = 1,
    aft = -shape
  )

  # ---- draw latent failure times (for all n subjects) -------------------- #
  x_frame <- data.frame(eta = eta_timing)
  raw <- simsurv::simsurv(dist    = "weibull",
                          lambdas = lambda,
                          gammas  = shape,
                          x       = x_frame,
                          betas   = c(eta = beta_eta))
  latent_time <- raw$eventtime   # length n

  # ---- 1-stage: done ----------------------------------------------------- #
  if (design == "1stage") {
    return(list(
      id             = seq_len(n),
      latent_time    = latent_time,
      event_eligible = rep(TRUE, n),
      pi             = NULL,
      intercept      = NULL,
      lambda         = lambda,
      shape          = shape
    ))
  }

  # ---- 2-stage: mixture-cure membership ---------------------------------- #
  if (is.null(eta_occurrence))
    stop("'eta_occurrence' must be supplied when design = '2stage'.")
  if (!is.numeric(eta_occurrence) || length(eta_occurrence) != n)
    stop("'eta_occurrence' must be a numeric vector of length n (", n, ").")
  if (is.null(target_event_prop) ||
      !is.numeric(target_event_prop) || length(target_event_prop) != 1L ||
      target_event_prop <= 0 || target_event_prop >= 1)
    stop("'target_event_prop' must be a single number in (0, 1) ",
         "when design = '2stage'.")

  # Solve intercept a so that mean(plogis(a + eta_occurrence)) == target_event_prop
  f_root <- function(a) mean(stats::plogis(a + eta_occurrence)) - target_event_prop
  # Bracket: plogis(-Inf+eta) -> 0, plogis(+Inf+eta) -> 1; widen bracket safely
  bracket_lo <- -500 - max(abs(eta_occurrence))
  bracket_hi <-  500 + max(abs(eta_occurrence))
  root <- stats::uniroot(f_root, lower = bracket_lo, upper = bracket_hi,
                         tol = .Machine$double.eps^0.5)
  intercept <- root$root

  pi_i          <- stats::plogis(intercept + eta_occurrence)
  event_eligible <- as.logical(stats::rbinom(n, size = 1L, prob = pi_i))

  # Cured subjects get latent time Inf — they never experience the event
  latent_time[!event_eligible] <- Inf

  list(
    id             = seq_len(n),
    latent_time    = latent_time,
    event_eligible = event_eligible,
    pi             = pi_i,
    intercept      = intercept,
    lambda         = lambda,
    shape          = shape
  )
}
