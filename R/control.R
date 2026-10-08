# Value of a control parameter for d smoothing parameters: element d (the last
# one if the vector is shorter, so a scalar applies to all d)
resolve_by_d <- function(param, d) {
  as.integer(param[[min(as.integer(d), length(param))]])
}

# Default control parameters. N_draw, N_final and n_quad have one value per
# number d = 1, ..., 4 of smoothing parameters (see resolve_by_d()).
bets_control_defaults <- function() {
  list(
    N_iter_max = 30,
    N_draw     = c(128L, 256L, 512L, 1024L),
    N_final    = c(200L, 300L, 500L, 500L),
    nu0 = 3,
    psi0 = NULL,
    phi_min = 0.8,
    phi_max = 0.98,
    min_ess = NULL,
    is_df = 5,
    is_scale = 4,
    lr = 0.9,
    c_inflate_eta = 3,
    prior_models = NULL,
    integration = "auto",
    n_quad = c(21L, 21L, 9L, 7L),
    n_scan = 64L
  )
}

# Checks control and fills in the defaults
resolve_bets_control <- function(control = list()) {
  if (!is.list(control) || (length(control) > 0 &&
                            (is.null(names(control)) || any(names(control) == "")))) {
    stop("control must be a named list")
  }
  defaults <- bets_control_defaults()

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0) {
    stop(unknown_control_message(unknown, names(defaults)), call. = FALSE)
  }

  ctrl <- utils::modifyList(defaults, control)
  check_bets_control(ctrl)
  ctrl
}

# Checks the values of the control entries (the length of prior_models is
# checked in bets(), where the number of models is known)
check_bets_control <- function(ctrl) {
  # numeric, length in len, no NA, and ok(x) for all values
  valid <- function(x, ok, len = 1) {
    is.numeric(x) && length(x) %in% len && !anyNA(x) && isTRUE(all(ok(x)))
  }
  integer_from <- function(lo) function(x) is.finite(x) & x >= lo & x == round(x)
  positive <- function(x) is.finite(x) & x > 0
  bad <- function(name, what) stop(sprintf("control$%s must be %s", name, what), call. = FALSE)
  by_d <- ", or one value per number of smoothing parameters"

  if (!valid(ctrl$N_iter_max, integer_from(1))) bad("N_iter_max", "a positive integer")
  for (nm in c("N_draw", "N_final", "n_quad")) {
    if (!valid(ctrl[[nm]], integer_from(1), 1:4)) bad(nm, paste0("a positive integer", by_d))
  }
  if (!is.null(ctrl$min_ess) && !valid(ctrl$min_ess, positive, 1:4)) {
    bad("min_ess", paste0("NULL or a positive number", by_d))
  }
  if (!valid(ctrl$n_scan, integer_from(0))) bad("n_scan", "a non-negative integer")
  for (nm in c("is_df", "is_scale", "c_inflate_eta")) {
    if (!valid(ctrl[[nm]], positive)) bad(nm, "a positive number")
  }
  if (!valid(ctrl$lr, function(x) x > 0 & x <= 1)) bad("lr", "a number in (0, 1]")

  # nu0 > 2, so that the prior mean of \sigma^2 exists
  if (!valid(ctrl$nu0, function(x) is.finite(x) & x > 2)) bad("nu0", "a number greater than 2")
  if (!is.null(ctrl$psi0) && !valid(ctrl$psi0, positive)) bad("psi0", "NULL or a positive number")
  if (!valid(ctrl$phi_min, function(x) x > 0 & x < 1)) bad("phi_min", "a number in (0, 1)")
  if (!valid(ctrl$phi_max, function(x) x > ctrl$phi_min & x <= 1)) {
    bad("phi_max", "a number in (phi_min, 1]")
  }
  p <- ctrl$prior_models
  if (!is.null(p) && !valid(p, function(x) x >= 0 & sum(x) > 0, seq_along(p))) {
    bad("prior_models", "NULL or non-negative numbers with a positive sum")
  }

  if (!(is.character(ctrl$integration) && length(ctrl$integration) == 1 &&
        ctrl$integration %in% c("auto", "quadrature", "ais"))) {
    bad("integration", "one of \"auto\", \"quadrature\", \"ais\"")
  }
  invisible(TRUE)
}

# Error message for unknown control entries: suggests the closest valid name
# for likely typos, and lists the valid names
unknown_control_message <- function(unknown, valid) {
  lines <- vapply(unknown, function(u) {
    d <- utils::adist(u, valid, ignore.case = TRUE)[1, ]
    msg <- sprintf("Unknown control entry '%s'.", u)
    if (min(d) <= max(2, nchar(u) %/% 3)) {
      msg <- sprintf("%s Did you mean '%s'?", msg, valid[which.min(d)])
    }
    msg
  }, character(1))
  paste(c(lines, paste("Valid entries:", paste(valid, collapse = ", "))), collapse = "\n")
}
