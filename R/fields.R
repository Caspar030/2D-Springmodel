# ============================================================
# Spatial fields and geometry
# Lagrangian description: reference configuration
# ============================================================

# ------------------------------------------------------------
# 1. Reference distance between two nodes
# ------------------------------------------------------------

A0 <- function(X, Y, parameters = list()) {
  
  stopifnot(
    is.numeric(X),
    is.numeric(Y),
    length(X) == length(Y)
  )
  
  sqrt(sum((Y - X)^2))
}


# ------------------------------------------------------------
# 2. Current node position
# ------------------------------------------------------------

current_position <- function(X, u) {
  
  stopifnot(
    is.numeric(X),
    is.numeric(u),
    length(X) == length(u)
  )
  
  X + u
}


# ------------------------------------------------------------
# 3. Current distance between two nodes
# ------------------------------------------------------------

a_current <- function(X, Y, uX, uY,
                      parameters = list()) {
  
  x <- current_position(X, uX)
  y <- current_position(Y, uY)
  
  sqrt(sum((y - x)^2))
}


# ------------------------------------------------------------
# 4. Default material and force fields
# ------------------------------------------------------------

default_fields <- function(dim) {
  
  if (length(dim) != 1 ||
      !is.numeric(dim) ||
      !dim %in% 1:3) {
    stop("dim must be 1, 2 or 3.")
  }
  
  list(
    
    # Reference distance
    A0 = function(X, Y, parameters) {
      A0(X, Y, parameters)
    },
    
    # Elasticity tensor E(X)
    E = function(X, parameters) {
      diag(
        vapply(
          seq_len(dim),
          function(k) {
            parameters[[paste0("Exyz", k)]] %||%
              c(1, 10, 10)[k]
          },
          numeric(1)
        ),
        nrow = dim
      )
    },
    
    # Viscosity tensor D(X)
    D = function(X, parameters) {
      diag(
        vapply(
          seq_len(dim),
          function(k) {
            parameters[[paste0("Dxyz", k)]] %||%
              c(1, 10, 10)[k]
          },
          numeric(1)
        ),
        nrow = dim
      )
    },
    
    # Local damping coefficient
    ceta = function(X, parameters) {
      parameters$ceta %||% 0.8
    },
    
    # Local retraction coefficient
    kappa = function(X, parameters) {
      parameters$kappa %||% 0.05
    },
    
    # Active force vector P(X, t)
    P = function(X, t, parameters) {
      parameters$P %||% rep(0, dim)
    }
  )
}


# ------------------------------------------------------------
# 5. Helper: use default if parameter is NULL
# ------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}