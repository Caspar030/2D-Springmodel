# R/fields.R
# Reference configuration; vector-valued displacement.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


# ------------------------------------------------------------
# Geometry
# ------------------------------------------------------------

check_vector <- function(x, name) {
  if (!is.numeric(x) ||
      !is.null(dim(x)) ||
      length(x) == 0L ||
      any(!is.finite(x))) {
    stop(name, " must be a finite numeric vector.")
  }
  
  invisible(TRUE)
}

A0 <- function(X, Y, parameters = list()) {
  check_vector(X, "X")
  check_vector(Y, "Y")
  
  if (length(X) != length(Y)) {
    stop("X and Y must have the same dimension.")
  }
  
  sqrt(sum((Y - X)^2))
}

current_position <- function(X, u) {
  check_vector(X, "X")
  check_vector(u, "u")
  
  if (length(X) != length(u)) {
    stop("X and u must have the same dimension.")
  }
  
  X + u
}

a_current <- function(X, Y, uX, uY,
                      parameters = list()) {
  x <- current_position(X, uX)
  y <- current_position(Y, uY)
  
  if (length(x) != length(y)) {
    stop("Both nodes must have the same spatial dimension.")
  }
  
  sqrt(sum((y - x)^2))
}


# ------------------------------------------------------------
# Material tensors
# ------------------------------------------------------------

validate_material <- function(C, dimension, name) {
  expected <- rep(as.integer(dimension), 4L)
  
  if (!is.numeric(C) ||
      !identical(dim(C), expected) ||
      any(!is.finite(C))) {
    stop(
      name, " must be a finite numeric array with dimensions ",
      paste(expected, collapse = " x "), "."
    )
  }
  
  C
}

# Isotropic example:
# C[i,j,k,l] =
#   lambda * delta[i,j] * delta[k,l]
#   + mu * (delta[i,k]*delta[j,l] + delta[i,l]*delta[j,k])
isotropic_tensor <- function(dimension, lambda, mu) {
  if (!is.numeric(dimension) ||
      length(dimension) != 1L ||
      !is.finite(dimension) ||
      !(dimension %in% 1:3)) {
    stop("dimension must be 1, 2 or 3.")
  }
  
  for (value in list(lambda, mu)) {
    if (!is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value)) {
      stop("lambda and mu must be finite numeric scalars.")
    }
  }
  
  C <- array(0, dim = rep(as.integer(dimension), 4L))
  
  for (i in seq_len(dimension)) {
    for (j in seq_len(dimension)) {
      for (k in seq_len(dimension)) {
        for (l in seq_len(dimension)) {
          C[i, j, k, l] <-
            lambda * (i == j) * (k == l) +
            mu * (
              (i == k) * (j == l) +
                (i == l) * (j == k)
            )
        }
      }
    }
  }
  
  C
}


# ------------------------------------------------------------
# Default fields
#
# parameters$E and parameters$D may each be:
#   - a fourth-order numeric array
#   - a function(X, parameters) returning such an array
#
# Scalar coefficients may be numbers or functions(X, parameters).
# P may be a vector or a function(X, t, parameters).
# ------------------------------------------------------------

default_fields <- function(dim) {
  if (!is.numeric(dim) ||
      length(dim) != 1L ||
      !is.finite(dim) ||
      !(dim %in% 1:3)) {
    stop("dim must be 1, 2 or 3.")
  }
  
  dimension <- as.integer(dim)
  
  check_position <- function(X) {
    check_vector(X, "X")
    
    if (length(X) != dimension) {
      stop("Position dimension does not match material dimension.")
    }
  }
  
  material_field <- function(name, lambda_default, mu_default) {
    force(name)
    force(lambda_default)
    force(mu_default)
    
    function(X, parameters = list()) {
      check_position(X)
      
      C <- parameters[[name]]
      
      if (is.null(C)) {
        C <- isotropic_tensor(
          dimension,
          lambda = parameters[[paste0("lambda_", name)]] %||%
            lambda_default,
          mu = parameters[[paste0("mu_", name)]] %||%
            mu_default
        )
      } else if (is.function(C)) {
        C <- C(X, parameters)
      }
      
      validate_material(C, dimension, name)
    }
  }
  
  scalar_field <- function(name, default) {
    force(name)
    force(default)
    
    function(X, parameters = list()) {
      check_position(X)
      
      value <- parameters[[name]] %||% default
      
      if (is.function(value)) {
        value <- value(X, parameters)
      }
      
      if (!is.numeric(value) ||
          length(value) != 1L ||
          !is.finite(value)) {
        stop(name, " must evaluate to a finite scalar.")
      }
      
      value
    }
  }
  
  list(
    A0 = function(X, Y, parameters = list()) {
      A0(X, Y, parameters)
    },
    
    # Illustrative isotropic defaults; replace for your material.
    E = material_field("E", lambda_default = 1, mu_default = 1),
    D = material_field("D", lambda_default = 1, mu_default = 1),
    
    rho = scalar_field("rho", 1),
    ceta = scalar_field("ceta", 0.8),
    kappa = scalar_field("kappa", 0.05),
    
    P = function(X, t, parameters = list()) {
      check_position(X)
      
      value <- parameters$P %||% rep(0, dimension)
      
      if (is.function(value)) {
        value <- value(X, t, parameters)
      }
      
      check_vector(value, "P")
      
      if (length(value) != dimension) {
        stop("P must have one component per spatial dimension.")
      }
      
      value
    }
  )
}