# Value of a control parameter for d smoothing parameters: a scalar applies to
# all d; from a vector, the element named d or else element d (the last one if
# the vector is shorter)
resolve_by_d <- function(param, d) {
  if (length(param) == 1L) return(as.integer(param))
  if (!is.null(names(param))) {
    key <- as.character(d)
    if (key %in% names(param)) return(as.integer(param[[key]]))
  }
  idx <- min(as.integer(d), length(param))
  as.integer(param[[idx]])
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

# Seasonal period of y, with the rules of forecast::ets(). 
# A non-integer frequency gives period 1
seasonal_period <- function(y) {
  if (stats::frequency(y) < 1) {
    warning("Frequency below 1 is treated as 1. Only non-seasonal models will be considered.",
            call. = FALSE)
    return(1L)
  }
  m <- stats::frequency(y)
  if (abs(m - round(m)) > 1e-4) {
    warning("Non-integer seasonal period. Only non-seasonal models will be considered.",
            call. = FALSE)
    return(1L)
  }
  as.integer(round(m))
}

# Components of a model as logicals: list(trend, seas, damped)
model_flags <- function(model_components) {
  list(trend  = model_components[[2]] == "A",
       seas   = model_components[[3]] == "A",
       damped = model_components[[4]] == "TRUE")
}

# Names of the smoothing parameters of a model
theta_names_of <- function(model_components) {
  flags <- model_flags(model_components)
  c("alpha", if (flags$trend) c("beta", if (flags$damped) "phi"), if (flags$seas) "gamma")
}

# Number of initial states of a model: l, [b,] [m seasonal states]
n_states <- function(model_components, m) {
  flags <- model_flags(model_components)
  1L + flags$trend + (if (flags$seas) m else 0L)
}

# Converts the model specification into a list of (error, trend, season, damped)
# vectors. model: an ETS code (e.g. "AAdN"; "Z" for any option, e.g. "ZZZ" or
# "AZN"), or a vector or list of codes. m: seasonal period, n: length of the series.
# As in forecast::ets(), seasonal models need 1 < m <= 24 and n > m: otherwise
# the seasonal options of "Z" are dropped, and explicit seasonal models are an error.
coerce_model_components <- function(model, m, additive.only = TRUE, n = Inf) {
  if (is.null(model)) model <- "ZZZ"
  if (is.list(model)) model <- unlist(model)
  if (!is.character(model) || length(model) == 0) {
    stop("model must be an ETS code (e.g. 'AAdN') or a vector of codes")
  }
  specs <- lapply(model, parse_model_spec)

  seasonal_ok <- m > 1 && m <= 24 && n > m
  if (m > 1 && !seasonal_ok && any(vapply(specs, function(s) s$season == "Z", logical(1)))) {
    warning(if (m > 24) "Seasonal models are not supported for frequency(y) > 24. "
            else "The number of observations is not sufficient for seasonal models. ",
            "Only non-seasonal models will be considered.", call. = FALSE)
  }

  models <- unlist(lapply(specs, expand_model_spec, m = m, n = n,
                          additive.only = additive.only, seasonal_ok = seasonal_ok),
                   recursive = FALSE)
  dup <- duplicated(models)
  if (any(dup)) {
    warning("Duplicate models removed: ",
            paste(unique(vapply(models[dup], ets_label, character(1))), collapse = ", "),
            call. = FALSE)
    models <- models[!dup]
  }
  models
}

# Components of one ETS code (e.g. "AAdN", "AZN"): list(error, trend, season,
# damped), with damped NA (both) for trend "Z"
parse_model_spec <- function(code) {
  code <- gsub("\\s+", "", toupper(code))
  parts <- if (nchar(code) == 3) {
    strsplit(code, "")[[1]]
  } else if (nchar(code) == 4 && substr(code, 3, 3) == "D") {
    c(substr(code, 1, 1), substr(code, 2, 3), substr(code, 4, 4))
  } else {
    stop("A model code must have 3 characters (e.g. 'ANN'), or 4 for a damped trend (e.g. 'AAdN')")
  }
  if (!parts[1] %in% c("A", "M", "Z")) {
    stop("Invalid error component: allowed values are 'A', 'M', 'Z'")
  }
  if (!parts[2] %in% c("N", "A", "M", "Z", "AD", "MD")) {
    stop("Invalid trend component: allowed values are 'N', 'A', 'M', 'Ad', 'Md', 'Z'")
  }
  if (!parts[3] %in% c("N", "A", "M", "Z")) {
    stop("Invalid seasonal component: allowed values are 'N', 'A', 'M', 'Z'")
  }
  trend <- substr(parts[2], 1, 1)
  damped <- if (trend == "Z") NA else if (nchar(parts[2]) == 2) "TRUE" else "FALSE"
  list(error = parts[1], trend = trend, season = parts[3], damped = damped)
}

# Models of one parsed code: each "Z" gives all its allowed options; explicit
# components that are not allowed are an error
expand_model_spec <- function(spec, m, n, additive.only, seasonal_ok) {
  if (additive.only && "M" %in% c(spec$error, spec$trend, spec$season)) {
    stop("Multiplicative components are not allowed when additive.only = TRUE")
  }
  if (spec$season %in% c("A", "M")) {
    if (m <= 1) stop("Seasonal models require frequency(y) to be an integer > 1")
    if (m > 24) stop("Seasonal models are not supported for frequency(y) > 24")
    if (n <= m) stop("Seasonal models require more than frequency(y) observations")
  }

  types <- if (additive.only) "A" else c("A", "M")
  options_of <- function(x, any) if (x == "Z") any else x
  # The first column varies fastest: "ZZZ" gives ANN, AAN, AAdN, ANA, AAA, AAdA
  grid <- expand.grid(error  = options_of(spec$error, types),
                      damped = if (is.na(spec$damped)) c("FALSE", "TRUE") else spec$damped,
                      trend  = options_of(spec$trend, c("N", types)),
                      season = options_of(spec$season, c("N", if (seasonal_ok) types)),
                      stringsAsFactors = FALSE)
  grid <- grid[!(grid$trend == "N" & grid$damped == "TRUE"), c("error", "trend", "season", "damped")]
  lapply(seq_len(nrow(grid)), function(i) unname(unlist(grid[i, ])))
}
