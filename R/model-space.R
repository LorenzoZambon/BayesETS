# Resolve a dimension-indexed control parameter.
# param may be:
#   - a scalar  → used for all theta dimensions
#   - an unnamed vector → param[min(d, length)] (last value repeated for d > length)
#   - a named vector with names "1","2","3","4" → look up by d
# Always returns a single integer.
resolve_by_d <- function(param, d) {
  if (length(param) == 1L) return(as.integer(param))
  if (!is.null(names(param))) {
    key <- as.character(d)
    if (key %in% names(param)) return(as.integer(param[[key]]))
  }
  idx <- min(as.integer(d), length(param))
  as.integer(param[[idx]])
}

bets_control_defaults <- function(freq = 1) {
  # N_draw and N_final are dimension-indexed vectors (index = theta dim d = 1..4).
  # Rule: N_draw = 256 * 2^(d-1); N_final = min(1000, 200 * 2^(d-1))
  # Pass a scalar to override uniformly; pass a length-4 vector for per-d control.
  # Both are resolved to a scalar inside each sampler via resolve_by_d(param, d).
  list(
    N_iter_max = 30,
    N_draw     = c(256L, 512L, 1024L, 2048L),
    N_draw_max = 1e5,
    N_final    = c(200L, 400L, 800L, 1000L),
    nu0 = 3,
    psi0 = NULL,
    phi_min = 0.8,
    phi_max = 0.98,
    min_ess = NULL,
    eta_df = 7,
    eta_df_incr_per_iter = 0,
    lr = 0.9,
    c_inflate_eta = 1,
    N_draw_mult = 1.5,
    first_iter_mult_N = 10,
    factor_inflate_Sigma = 1,
    verbose = 0,
    n_traj_forecast = 1000,
    prior_models = NULL,
    K_mix = 3,
    ridge_eps = 1e-4,
    em_iter = 5,
    jitter_scale = 0.5,
    n_sobol = NULL
  )
}

resolve_bets_control <- function(control = list(), freq = 1) {
  defaults <- bets_control_defaults(freq)
  if (!is.list(control)) {
    stop("control must be a named list")
  }
  # method, sampler live in control but are not tuning params
  bets_keys <- c("method", "sampler")
  unknown <- setdiff(names(control), c(names(defaults), bets_keys))
  if (length(unknown) > 0) {
    stop(sprintf("Unknown control entries: %s", paste(unknown, collapse = ", ")))
  }
  control <- control[setdiff(names(control), bets_keys)]

  out <- utils::modifyList(defaults, control)
  # min_ess and n_sobol are resolved per-model inside each sampler (after d is known):
  #   min_ess defaults to N_final[d] / 2
  #   n_sobol defaults to N_draw[d]
  out
}

bets_model_space <- function(m) {
  if (m > 1) {
    list(
      c("A", "N", "N", "FALSE"),
      c("A", "A", "N", "FALSE"),
      c("A", "A", "N", "TRUE"),
      c("A", "N", "A", "FALSE"),
      c("A", "A", "A", "FALSE"),
      c("A", "A", "A", "TRUE")
    )
  } else {
    list(
      c("A", "N", "N", "FALSE"),
      c("A", "A", "N", "FALSE"),
      c("A", "A", "N", "TRUE")
    )
  }
}

coerce_model_components <- function(model, m, additive.only = TRUE) {
  if (!is.logical(additive.only) || length(additive.only) != 1 || is.na(additive.only)) {
    stop("additive.only must be TRUE or FALSE")
  }

  normalize_damped <- function(damped) {
    if (is.logical(damped) && length(damped) == 1 && !is.na(damped)) {
      return(if (damped) "TRUE" else "FALSE")
    }
    if (!is.character(damped) || length(damped) != 1) {
      stop("Invalid damped component: use TRUE/FALSE (or T/F)")
    }
    damped_upper <- toupper(damped)
    switch(
      damped_upper,
      "TRUE" = "TRUE",
      "FALSE" = "FALSE",
      "T" = "TRUE",
      "F" = "FALSE",
      stop("Invalid damped component: use TRUE/FALSE (or T/F)")
    )
  }

  normalize_components <- function(error_comp, trend_comp, season_comp, damped = NULL) {
    error_comp <- toupper(error_comp)
    trend_comp <- toupper(trend_comp)
    season_comp <- toupper(season_comp)

    if (trend_comp %in% c("AD", "MD")) {
      trend_base <- substr(trend_comp, 1, 1)
      damped_from_trend <- "TRUE"
    } else if (trend_comp %in% c("A", "M", "N")) {
      trend_base <- trend_comp
      damped_from_trend <- "FALSE"
    } else {
      stop("Invalid trend component: allowed values are 'N', 'A', 'M', 'Ad', 'Md'")
    }

    if (!error_comp %in% c("A", "M")) {
      stop("Invalid error component: allowed values are 'A' or 'M'")
    }
    if (!season_comp %in% c("N", "A", "M")) {
      stop("Invalid seasonal component: allowed values are 'N', 'A', 'M'")
    }

    if (m <= 1 && season_comp != "N") {
      stop("Seasonal models require frequency(y) > 1")
    }

    damped_norm <- if (is.null(damped)) damped_from_trend else normalize_damped(damped)

    if (trend_base == "N" && damped_norm == "TRUE") {
      stop("Damped trend is only valid when trend component is 'A' or 'M'")
    }

    if (additive.only && any(c(error_comp, trend_base, season_comp) == "M")) {
      stop("Multiplicative components are not allowed when additive.only = TRUE")
    }

    c(error_comp, trend_base, season_comp, damped_norm)
  }

  normalize_one_model <- function(x) {
    if (is.character(x) && length(x) == 1) {
      code <- gsub("\\s+", "", toupper(x))
      if (identical(code, "ZZZ")) {
        stop("'ZZZ' is only allowed as top-level model input")
      }
      if (nchar(code) == 3) {
        return(normalize_components(
          substr(code, 1, 1),
          substr(code, 2, 2),
          substr(code, 3, 3)
        ))
      }
      if (nchar(code) == 4 && substr(code, 3, 3) == "D") {
        return(normalize_components(
          substr(code, 1, 1),
          paste0(substr(code, 2, 2), "D"),
          substr(code, 4, 4)
        ))
      }
      stop("model string must be 'ZZZ', a 3-character ETS code (e.g. 'ANN'), or a 4-character damped code (e.g. 'AAdN')")
    }

    if (is.character(x) && length(x) == 3) {
      return(normalize_components(x[[1]], x[[2]], x[[3]]))
    }

    if (is.character(x) && length(x) == 4) {
      return(normalize_components(x[[1]], x[[2]], x[[3]], x[[4]]))
    }

    stop("Each model specification must be a 3/4-character ETS string or a 3/4-component character vector")
  }

  if (is.null(model) || identical(model, "ZZZ")) {
    return(bets_model_space(m))
  }

  if (is.character(model) && length(model) %in% c(1, 3, 4)) {
    return(list(normalize_one_model(model)))
  }

  if (is.list(model)) {
    return(lapply(model, normalize_one_model))
  }

  stop("model must be 'ZZZ', a 3/4-character ETS code, a 3/4-component character vector, or a list of such specifications")
}
