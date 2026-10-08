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

#' Label of an ETS model (e.g. "AAdN")
#'
#' @param model_components Vector (error, trend, season, damped).
#' @return Character string.
#' @keywords internal
ets_label <- function(model_components) {
  d <- if (model_components[[4]] == "TRUE") "d" else ""
  paste0(model_components[[1]], model_components[[2]], d, model_components[[3]])
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
