#' Apply censoring to latent failure times
#'
#' Takes a vector of latent failure times (output of \code{\link{sim_tte}}) and
#' produces the final observed survival outcome \code{(time, status)} by
#' generating latent censoring times and taking the minimum of each subject's
#' failure time and censoring time.
#'
#' Four censoring mechanisms are supported (via \code{type}):
#' \describe{
#'   \item{\code{"none"}}{Censoring time is \code{max(finite latent_time) + 10}
#'     for every subject — effectively no censoring.  Only \code{Inf} latent
#'     times (cured subjects from a 2-stage design) will be censored.}
#'   \item{\code{"administrative"}}{Every subject is censored at the fixed
#'     calendar time \code{administrative_time}.}
#'   \item{\code{"random"}}{Independent censoring times are drawn from the
#'     distribution specified in \code{random} (currently
#'     \code{random$dist = "exponential"} with rate \code{random$rate}).  If
#'     \code{administrative_time} is also supplied, draw times are capped at
#'     that value.}
#'   \item{\code{"administrative_plus_random"}}{Same as \code{"random"} but
#'     \code{administrative_time} is required and always applied as the cap.}
#' }
#'
#' @section \code{Inf} latent times (mixture-cure):
#' Subjects with \code{latent_time = Inf} (cured subjects from a 2-stage
#' design) are always censored regardless of type: their observed time equals
#' the censoring time and their status is 0.
#'
#' @section Optional calibration via \code{target_event_prop}:
#' When \code{target_event_prop} is supplied and the chosen \code{type}
#' involves an administrative cutoff, \code{administrative_time} is solved
#' automatically as the \code{target_event_prop}-quantile of the
#' \emph{finite} latent failure times.  This ensures that approximately
#' \code{target_event_prop} of subjects (among those with finite latent times)
#' become events at the administrative cutoff.  The argument is ignored for
#' \code{type = "none"}.  Note that in a 2-stage design, cured subjects
#' (\code{Inf} latent times) are always censored and therefore reduce the
#' realised event rate below \code{target_event_prop}.
#'
#' @param latent_time Numeric vector of length \eqn{n}.  Latent failure times
#'   as returned by the \code{latent_time} component of \code{\link{sim_tte}}.
#'   \code{Inf} values for cured subjects are handled correctly.
#' @param type Character scalar.  One of \code{"none"},
#'   \code{"administrative"}, \code{"random"}, or
#'   \code{"administrative_plus_random"}.
#' @param administrative_time Positive numeric scalar.  Fixed administrative
#'   censoring time.  Required for \code{type = "administrative"} and
#'   \code{"administrative_plus_random"}.  Optional cap for \code{"random"}.
#'   Overridden by \code{target_event_prop} calibration when supplied.
#' @param random Named list controlling random censoring.  Used fields:
#'   \code{dist} (character scalar; currently \code{"exponential"} only) and
#'   \code{rate} (positive numeric scalar; the exponential rate parameter).
#'   Required when \code{type} is \code{"random"} or
#'   \code{"administrative_plus_random"}.
#' @param target_event_prop Numeric scalar in \eqn{(0, 1)} or \code{NULL}.
#'   When non-\code{NULL} and \code{type} involves an administrative cutoff,
#'   \code{administrative_time} is set to the corresponding empirical quantile
#'   of the finite latent failure times.  See the calibration section above.
#' @param ... Reserved for future arguments (currently ignored).
#'
#' @return A \code{data.frame} with one row per subject and columns:
#'   \describe{
#'     \item{\code{id}}{Integer: subject index \code{1:n}.}
#'     \item{\code{time}}{Numeric: observed event or censoring time
#'       \eqn{= \min(T_i, C_i)}.}
#'     \item{\code{status}}{Integer: \code{1} if the event was observed
#'       (\eqn{T_i \le C_i}), \code{0} if censored.}
#'     \item{\code{latent_time}}{Numeric: the original latent failure time
#'       passed in (including any \code{Inf} values).}
#'     \item{\code{censor_time}}{Numeric: the latent censoring time applied
#'       to each subject.}
#'   }
#'
#' @importFrom stats quantile rexp
#'
#' @seealso \code{\link{sim_tte}}, \code{\link{simulate_tte}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' tte_out <- sim_tte(eta_timing = rnorm(200), design = "1stage",
#'                    model = "ph", dist = "weibull",
#'                    params = list(shape = 1.5), target_median = 5)
#'
#' # Administrative censoring at t = 8
#' obs <- sim_censor(tte_out$latent_time, type = "administrative",
#'                   administrative_time = 8)
#'
#' # Random exponential censoring, no admin cap
#' obs2 <- sim_censor(tte_out$latent_time, type = "random",
#'                    random = list(dist = "exponential", rate = 0.15))
#'
#' # Calibrate admin time to hit 60 % event rate
#' obs3 <- sim_censor(tte_out$latent_time, type = "administrative",
#'                    target_event_prop = 0.60)
#' }
#'
#' @export
sim_censor <- function(latent_time,
                       type                = "administrative",
                       administrative_time = NULL,
                       random              = list(dist = "exponential",
                                                  rate = NULL),
                       target_event_prop   = NULL,
                       ...) {

  # ---- input validation -------------------------------------------------- #
  if (!is.numeric(latent_time) || length(latent_time) < 1L)
    stop("'latent_time' must be a non-empty numeric vector.")
  n <- length(latent_time)

  finite_times <- latent_time[is.finite(latent_time)]
  if (length(finite_times) == 0L)
    stop("All 'latent_time' values are Inf; no finite failure times to censor.")

  # ---- optional calibration of administrative_time ----------------------- #
  # Applies whenever type involves an admin cutoff and target_event_prop given.
  # Calibration by quantile only works for types where ALL censoring is
  # deterministically at the admin cutoff.  For "random" and
  # "administrative_plus_random" the random draws produce additional censoring
  # below the cap, so the realised event rate would be lower than the target.
  uses_admin <- type %in% c("administrative")
  if (!is.null(target_event_prop) && uses_admin) {
    if (!is.numeric(target_event_prop) || length(target_event_prop) != 1L ||
        target_event_prop <= 0 || target_event_prop >= 1)
      stop("'target_event_prop' must be a single number in (0, 1).")
    administrative_time <- stats::quantile(finite_times,
                                           probs = target_event_prop,
                                           type  = 1,
                                           names = FALSE)
  }

  # ---- generate latent censoring times ----------------------------------- #
  censor_time <- switch(
    type,

    none = {
      rep(max(finite_times) + 10, n)
    },

    administrative = {
      if (is.null(administrative_time))
        stop("'administrative_time' must be supplied when type = 'administrative'.")
      rep(as.numeric(administrative_time), n)
    },

    random = {
      .draw_random_censor(n, random, administrative_time)
    },

    administrative_plus_random = {
      if (is.null(administrative_time))
        stop("'administrative_time' must be supplied when ",
             "type = 'administrative_plus_random'.")
      .draw_random_censor(n, random, administrative_time)
    },

    stop("type '", type, "' not yet implemented")
  )

  # ---- assemble observed outcome ----------------------------------------- #
  time   <- pmin(latent_time, censor_time)
  status <- as.integer(latent_time <= censor_time)

  data.frame(
    id          = seq_len(n),
    time        = as.numeric(time),
    status      = status,
    latent_time = as.numeric(latent_time),
    censor_time = as.numeric(censor_time)
  )
}


# Draw random censoring times from random$dist, capping at admin_time if given.
.draw_random_censor <- function(n, random, admin_time) {
  dist <- random$dist %||% stop("'random$dist' must be specified.")

  C <- switch(
    dist,
    exponential = {
      rate <- random$rate
      if (is.null(rate) || !is.numeric(rate) || length(rate) != 1L || rate <= 0)
        stop("'random$rate' must be a single positive number for ",
             "dist = 'exponential'.")
      stats::rexp(n, rate = rate)
    },
    stop("random$dist '", dist, "' not yet implemented")
  )

  if (!is.null(admin_time)) {
    C <- pmin(C, as.numeric(admin_time))
  }

  C
}
