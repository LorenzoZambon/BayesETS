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
bets_control_defaults <- function(freq = 1) {
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
resolve_bets_control <- function(control = list(), freq = 1) {
  if (!is.list(control) || (length(control) > 0 &&
                            (is.null(names(control)) || any(names(control) == "")))) {
    stop("control must be a named list")
  }
  defaults <- bets_control_defaults(freq)

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0) {
    stop(unknown_control_message(unknown, names(defaults)), call. = FALSE)
  }

  ctrl <- utils::modifyList(defaults, control)
  # nu0 > 2, so that the prior mean of \sigma^2 exists
  if (!is.numeric(ctrl$nu0) || length(ctrl$nu0) != 1 || !(ctrl$nu0 > 2)) {
    stop("control$nu0 must be a single number greater than 2")
  }
  ctrl
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

# Seasonal period of y, with the rules of forecast::ets(). A non-integer
# frequency gives period 1: rounding it would make the seasonal pattern drift.
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

# Models of "ZZZ" for seasonal period m
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

# Converts the model specification into a list of (error, trend, season,
# damped) vectors. m: seasonal period, n: length of the series. As in
# forecast::ets(), seasonal models need 1 < m <= 24 and n > m.
coerce_model_components <- function(model, m, additive.only = TRUE, n = Inf) {
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

    if (season_comp != "N") {
      if (m <= 1) stop("Seasonal models require frequency(y) to be an integer > 1")
      if (m > 24) stop("Seasonal models are not supported for frequency(y) > 24")
      if (n <= m) stop("Seasonal models require more than frequency(y) observations")
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
    if (m > 24) {
      warning("Seasonal models are not supported for frequency(y) > 24. ",
              "Only non-seasonal models will be considered.", call. = FALSE)
    }
    return(bets_model_space(if (m <= 24 && n > m) m else 1L))
  }

  if (is.character(model) && length(model) %in% c(1, 3, 4)) {
    return(list(normalize_one_model(model)))
  }

  if (is.list(model)) {
    return(lapply(model, normalize_one_model))
  }

  stop("model must be 'ZZZ', a 3/4-character ETS code, a 3/4-component character vector, or a list of such specifications")
}
