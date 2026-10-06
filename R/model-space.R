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
  # Rule: N_draw = 128 * 2^(d-1). N_final (posterior draws kept per model) keeps
  # their share of the forecast-interval error at about 1% of the 95% interval
  # width for every d (well below the noise of the default 1000 trajectories).
  # Pass a scalar to override uniformly; pass a length-4 vector for per-d control.
  # Both are resolved to a scalar inside each sampler via resolve_by_d(param, d).
  # integration = "auto" uses quadrature for d <= 2 and AIS for d >= 3.
  # AIS (adaptive_is_rb) draws N_draw points per iteration from a Student-t
  # proposal with is_df degrees of freedom and scale is_scale * H^{-1}, and stops
  # once ESS >= min_ess (NULL -> N_draw / 4). n_quad = Gauss-Hermite nodes per
  # dimension for quadrature_rb, also indexed by d, so the grid has n_quad[d]^d
  # nodes (21, 441, 729, 2401). is_scale and n_quad[1:2] were tuned on M3 and
  # tourism series (2026-10). The combination (BMA/stacking) and verbose are
  # explicit arguments of bets().
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
    c_inflate_eta = 1,
    prior_models = NULL,
    integration = "auto",
    n_quad = c(21L, 21L, 9L, 7L)
  )
}

resolve_bets_control <- function(control = list(), freq = 1) {
  if (!is.list(control)) {
    stop("control must be a named list")
  }
  defaults <- bets_control_defaults(freq)

  # Check if control contains any unknown entries (not in defaults) and if so issue a warning and ignore them
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0) {
    warning(sprintf("Unknown control entries: %s. They will be ignored.", paste(unknown, collapse = ", ")))
    control <- control[setdiff(names(control), unknown)]
  }

  return(utils::modifyList(defaults, control))
}

bets_model_space <- function(freq) {
  if (freq > 1) {
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

coerce_model_components <- function(model, freq, additive.only = TRUE) {
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

    if (freq <= 1 && season_comp != "N") {
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
    return(bets_model_space(freq))
  }

  if (is.character(model) && length(model) %in% c(1, 3, 4)) {
    return(list(normalize_one_model(model)))
  }

  if (is.list(model)) {
    return(lapply(model, normalize_one_model))
  }

  stop("model must be 'ZZZ', a 3/4-character ETS code, a 3/4-component character vector, or a list of such specifications")
}
