# ============================================================
# Operators for the expression tree
# Requires expression_tree.R
# ============================================================


# ------------------------------------------------------------
# 1. Gradient
# ------------------------------------------------------------

grad <- function(x) {
  
  stopifnot(inherits(x, "math_expr"))
  
  new_expr(
    op    = "grad",
    args  = list(x),
    shape = c(x$shape, 2)
  )
}


# ------------------------------------------------------------
# 2. Divergence
# Contracts the last tensor dimension with the spatial index
# ------------------------------------------------------------

div <- function(x) {
  
  stopifnot(inherits(x, "math_expr"))
  
  if (length(x$shape) == 0) {
    stop("Type error: divergence requires a vector or tensor.")
  }
  
  if (tail(x$shape, 1) != 2) {
    stop("Type error: divergence requires a spatial dimension of 2.")
  }
  
  new_expr(
    op    = "div",
    args  = list(x),
    shape = head(x$shape, -1)
  )
}


# ------------------------------------------------------------
# 3. Matrix multiplication
# ------------------------------------------------------------

matmul <- function(A, B) {
  
  stopifnot(
    inherits(A, "math_expr"),
    inherits(B, "math_expr")
  )
  
  if (length(A$shape) != 2) {
    stop("Type error: left operand must be a matrix.")
  }
  
  if (length(B$shape) == 0) {
    stop("Type error: right operand must be a vector or matrix.")
  }
  
  if (A$shape[2] != B$shape[1]) {
    stop("Type error: incompatible dimensions for multiplication.")
  }
  
  new_expr(
    op    = "matmul",
    args  = list(A, B),
    shape = c(A$shape[1], tail(B$shape, -1))
  )
}


# ------------------------------------------------------------
# 4. Addition
# ------------------------------------------------------------

add <- function(x, y) {
  
  stopifnot(
    inherits(x, "math_expr"),
    inherits(y, "math_expr")
  )
  
  if (!identical(x$shape, y$shape)) {
    stop("Type error: addition requires identical tensor shapes.")
  }
  
  new_expr(
    op    = "add",
    args  = list(x, y),
    shape = x$shape
  )
}


# ------------------------------------------------------------
# 5. Multiplication
# Scalar * scalar, or scalar * vector/matrix/tensor
# ------------------------------------------------------------

multiply <- function(x, y) {
  
  stopifnot(
    inherits(x, "math_expr"),
    inherits(y, "math_expr")
  )
  
  x_scalar <- length(x$shape) == 0
  y_scalar <- length(y$shape) == 0
  
  if (!x_scalar && !y_scalar) {
    stop("Type error: multiply requires at least one scalar.")
  }
  
  result_shape <- if (x_scalar) y$shape else x$shape
  
  new_expr(
    op    = "multiply",
    args  = list(x, y),
    shape = result_shape
  )
}