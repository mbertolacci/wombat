#' @export
to_month_start <- function(x) {
  x <- date(x)
  if_else(
    day(x) == 1,
    x + days(1) - months(1),
    x + days(2 - day(x))
  )
}

#' @export
gamma_quantile_prior <- function(
  q_lower, q_upper,
  p_lower = 0.05, p_upper = 0.95,
  interval = c(0.01, 10)
) {
  stopifnot(q_lower < q_upper)
  stopifnot(p_lower < p_upper)

  # Find the shape parameter that gives the appropriate ratio between the
  # quantiles
  ratio <- q_lower / q_upper
  shape <- uniroot(function(shape) {
    theoretical <- qgamma(c(p_lower, p_upper), shape = shape, rate = 1)
    theoretical[1] / theoretical[2] - ratio
  }, interval, tol = sqrt(.Machine$double.eps))$root

  # Find the rate parameter that gives the correct quantiles
  rate <- qgamma(p_upper, shape = shape, rate = 1) / q_upper

  list(shape = shape, rate = rate)
}

#' @export
generate <- function(x, ...) UseMethod('generate', x)

#' @export
calculate <- function(x, ...) {
  UseMethod('calculate', x)
}

#' @export
log_prior <- function(x, ...) UseMethod('log_prior', x)

#' @export
n_terms <- function(f) {
  length(attr(terms(f), 'term.labels'))
}

.rnorm_truncated <- function(
  n, mean = 0, sd = 1, lower = -Inf, upper = Inf, max_iterations = 1000
) {
  sapply(1 : n, function(i) {
    mean_i <- mean[1 + (i - 1) %% length(mean)]
    sd_i <- sd[1 + (i - 1) %% length(sd)]
    for (j in 1 : max_iterations) {
      output <- rnorm(1, mean_i, sd_i)
      if (output > lower && output < upper) {
        return(output)
      }
    }
    stop('max iterations exceeded')
  })
}

.log_debug <- function(...) {
  futile.logger::flog.debug(..., name = 'fluxcapacitor')
}

.log_trace <- function(...) {
  futile.logger::flog.trace(..., name = 'fluxcapacitor')
}

.chol_solve <- function(R, b) {
  if (is(R, 'triangularMatrix')) {
    solve(R, solve(t(R), b))
  } else if (is(R, 'CHMfactor')) {
    solve(R, b, system = 'A')
  } else {
    backsolve(R, backsolve(R, b, transpose = TRUE))
  }
}

# Given a series of lists as arguments, shallowly merge the lists from left to
# right (so that the right-most values override the left-most)
.extend_list <- function(...) {
  lists <- list(...)
  output <- lists[[1]]
  for (value in lists[2 : length(lists)]) {
    for (name in names(value)) {
      output[[name]] <- value[[name]]
    }
  }
  return(output)
}

.remove_nulls <- function(x) {
  Filter(Negate(is.null), x)
}

.remove_nulls_and_missing <- function(x) {
  Filter(Negate(is.symbol), Filter(Negate(is.null), x))
}

.strip_attributes <- function(x) {
  attributes(x) <- NULL
  x
}

.recycle_vector_to <- function(x, to_length) {
  x[
    ((seq_len(to_length) - 1) %% length(x)) + 1
  ]
}

.cache_memory_fifo <- function(algo = 'sha512', size = 10) {
  cache <- NULL
  keys <- NULL

  cache_reset <- function() {
    cache <<- new.env(TRUE, emptyenv())
    keys <<- NULL
  }

  cache_set <- function(key, value) {
    assign(key, value, envir = cache)
    keys <<- c(keys, key)
    if (length(keys) > size) {
      cache_drop_key(keys[[1]])
    }
  }

  cache_get <- function(key) {
    get(key, envir = cache, inherits = FALSE)
  }

  cache_has_key <- function(key) {
    exists(key, envir = cache, inherits = FALSE)
  }

  cache_drop_key <- function(key) {
    rm(list = key, envir = cache, inherits = FALSE)
    keys <<- keys[keys != key]
  }

  cache_reset()
  list(
    digest = function(...) digest::digest(..., algo = algo),
    reset = cache_reset,
    set = cache_set,
    get = cache_get,
    has_key = cache_has_key,
    drop_key = cache_drop_key,
    keys = function() ls(cache)
  )
}
