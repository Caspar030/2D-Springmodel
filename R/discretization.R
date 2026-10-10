# R/discretization.R
#
# Recursive symbolic discretization of differential operators.
# Coordinates refer to the reference configuration.
#
# Supported operators:
#   field, grad, div, matmul, add, multiply
#
# The output consists of symbolic expression strings.
#
# Finite differences:
#   Interior: centered difference, second-order accurate
#   Lower boundary: forward difference, second-order accurate
#   Upper boundary: backward difference, second-order accurate
#
# Boundary information is optional:
#   grid_bounds = list(min = c(1, 1), max = c(Nx, Ny))
#
# Without grid_bounds, centered differences are used everywhere.


# ------------------------------------------------------------
# 1. Grid utilities
# ------------------------------------------------------------

state_name <- function(variable, index) {
  paste(c(variable, as.integer(index)), collapse = "_")
}


neighbor_index <- function(index, dimension, direction) {
  result <- index
  result[dimension] <- result[dimension] + direction
  result
}


# ------------------------------------------------------------
# 2. Finite-difference stencils
# ------------------------------------------------------------

derivative_stencil <- function(
    index,
    dimension,
    grid_bounds = NULL
) {
  
  # Without boundary information, use centered differences.
  if (is.null(grid_bounds)) {
    return(list(
      points = list(
        neighbor_index(index, dimension, 1),
        neighbor_index(index, dimension, -1)
      ),
      coefficients = c(1, -1),
      denominator = 2
    ))
  }
  
  lower <- grid_bounds$min[dimension]
  upper <- grid_bounds$max[dimension]
  position <- index[dimension]
  
  if (position < lower || position > upper) {
    stop("Index lies outside grid_bounds.")
  }
  
  # Lower boundary: second-order forward difference.
  if (position == lower) {
    return(list(
      points = list(
        index,
        neighbor_index(index, dimension, 1),
        neighbor_index(index, dimension, 2)
      ),
      coefficients = c(-3, 4, -1),
      denominator = 2
    ))
  }
  
  # Upper boundary: second-order backward difference.
  if (position == upper) {
    return(list(
      points = list(
        index,
        neighbor_index(index, dimension, -1),
        neighbor_index(index, dimension, -2)
      ),
      coefficients = c(3, -4, 1),
      denominator = 2
    ))
  }
  
  # Interior: centered difference.
  list(
    points = list(
      neighbor_index(index, dimension, 1),
      neighbor_index(index, dimension, -1)
    ),
    coefficients = c(1, -1),
    denominator = 2
  )
}


# ------------------------------------------------------------
# 3. Symbolic field representation
# ------------------------------------------------------------

# Default symbolic name for a field at a reference-grid point.
#
# Examples:
#   u at (3,4)       -> u_3_4
#   E[1,2] at (3,4) -> E_1_2_3_4
#
# state_variables identifies scalar fields representing states.

default_field_symbol <- function(
    name,
    index,
    shape,
    component = integer(0),
    state_variables = c("u", "v"),
    fields = NULL,
    parameters = list()
) {
  
  if (name %in% state_variables && length(shape) == 0) {
    return(state_name(name, index))
  }
  
  # If the field is a known coefficient function, retain a
  # symbolic representation of its spatial evaluation.
  if (!is.null(fields) && name %in% names(fields)) {
    return(
      paste0(
        name,
        if (length(component)) {
          paste0("_", paste(component, collapse = "_"))
        },
        "_",
        paste(index, collapse = "_")
      )
    )
  }
  
  # Generic symbolic field.
  paste(
    c(name, component, as.integer(index)),
    collapse = "_"
  )
}


# ------------------------------------------------------------
# 4. Main discretization function
# ------------------------------------------------------------

discretize <- function(
    expr,
    index,
    spacing,
    fields = NULL,
    parameters = list(),
    field_symbol = default_field_symbol,
    state_variables = c("u", "v"),
    grid_bounds = NULL
) {
  
  if (!inherits(expr, "math_expr")) {
    stop("expr must be a math_expr.")
  }
  
  if (length(index) != length(spacing)) {
    stop("index and spacing must have the same dimension.")
  }
  
  if (any(!is.finite(spacing)) || any(spacing <= 0)) {
    stop("spacing must contain positive finite values.")
  }
  
  if (!is.null(grid_bounds)) {
    
    if (is.null(grid_bounds$min) ||
        is.null(grid_bounds$max)) {
      stop("grid_bounds must contain min and max.")
    }
    
    if (length(grid_bounds$min) != length(index) ||
        length(grid_bounds$max) != length(index)) {
      stop("grid_bounds dimensions must match index.")
    }
    
    if (any(grid_bounds$min >= grid_bounds$max)) {
      stop("Each grid dimension must contain at least two points.")
    }
    
    if (any(index < grid_bounds$min) ||
        any(index > grid_bounds$max)) {
      stop("index lies outside grid_bounds.")
    }
  }
  
  # Pass common arguments to recursive calls.
  recurse <- function(x, idx) {
    discretize(
      expr = x,
      index = idx,
      spacing = spacing,
      fields = fields,
      parameters = parameters,
      field_symbol = field_symbol,
      state_variables = state_variables,
      grid_bounds = grid_bounds
    )
  }
  
  
  # ----------------------------------------------------------
  # Field
  # ----------------------------------------------------------
  
  if (expr$op == "field") {
    
    shape <- expr$shape
    
    if (length(shape) == 0) {
      return(
        field_symbol(
          name = expr$name,
          index = index,
          shape = shape,
          component = integer(0),
          state_variables = state_variables,
          fields = fields,
          parameters = parameters
        )
      )
    }
    
    if (length(shape) == 1) {
      return(
        vapply(
          seq_len(shape[1]),
          function(k) {
            field_symbol(
              name = expr$name,
              index = index,
              shape = shape,
              component = k,
              state_variables = state_variables,
              fields = fields,
              parameters = parameters
            )
          },
          character(1)
        )
      )
    }
    
    if (length(shape) == 2) {
      
      result <- matrix(
        "",
        nrow = shape[1],
        ncol = shape[2]
      )
      
      for (i in seq_len(shape[1])) {
        for (j in seq_len(shape[2])) {
          result[i, j] <- field_symbol(
            name = expr$name,
            index = index,
            shape = shape,
            component = c(i, j),
            state_variables = state_variables,
            fields = fields,
            parameters = parameters
          )
        }
      }
      
      return(result)
    }
    
    stop("Fields with more than two tensor dimensions are unsupported.")
  }
  
  
  # ----------------------------------------------------------
  # Gradient
  # ----------------------------------------------------------
  
  if (expr$op == "grad") {
    
    child <- expr$args[[1]]
    
    if (length(child$shape) != 0) {
      stop("grad currently requires a scalar expression.")
    }
    
    return(
      vapply(
        seq_along(index),
        function(k) {
          
          stencil <- derivative_stencil(
            index = index,
            dimension = k,
            grid_bounds = grid_bounds
          )
          
          values <- vapply(
            seq_along(stencil$points),
            function(m) {
              recurse(child, stencil$points[[m]])
            },
            character(1)
          )
          
          terms <- paste0(
            "(",
            stencil$coefficients,
            ") * (",
            values,
            ")"
          )
          
          paste0(
            "(",
            paste(terms, collapse = " + "),
            ") / (",
            stencil$denominator,
            " * ",
            spacing[k],
            ")"
          )
        },
        character(1)
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Matrix-vector multiplication
  # ----------------------------------------------------------
  
  if (expr$op == "matmul") {
    
    left <- recurse(expr$args[[1]], index)
    right <- recurse(expr$args[[2]], index)
    
    if (is.null(dim(left)) || length(dim(left)) != 2) {
      stop("Left operand of matmul must be a matrix.")
    }
    
    if (length(right) != ncol(left)) {
      stop("Incompatible dimensions in matmul.")
    }
    
    return(
      vapply(
        seq_len(nrow(left)),
        function(i) {
          
          terms <- paste0(
            "(",
            left[i, ],
            ") * (",
            right,
            ")"
          )
          
          paste0("(", paste(terms, collapse = " + "), ")")
        },
        character(1)
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Divergence
  # ----------------------------------------------------------
  
  if (expr$op == "div") {
    
    child <- expr$args[[1]]
    dimension <- length(index)
    
    if (tail(child$shape, 1) != dimension) {
      stop(
        "Divergence requires a vector or tensor whose last ",
        "dimension matches the spatial dimension."
      )
    }
    
    return(
      paste0(
        "(",
        paste(
          vapply(
            seq_len(dimension),
            function(k) {
              
              stencil <- derivative_stencil(
                index = index,
                dimension = k,
                grid_bounds = grid_bounds
              )
              
              values <- lapply(
                stencil$points,
                function(idx) recurse(child, idx)
              )
              
              if (any(vapply(
                values,
                length,
                integer(1)
              ) != dimension)) {
                stop(
                  "Divergence currently requires a vector-valued flux."
                )
              }
              
              terms <- vapply(
                seq_along(values),
                function(m) {
                  paste0(
                    "(",
                    stencil$coefficients[m],
                    ") * (",
                    values[[m]][k],
                    ")"
                  )
                },
                character(1)
              )
              
              paste0(
                "(",
                paste(terms, collapse = " + "),
                ") / (",
                stencil$denominator,
                " * ",
                spacing[k],
                ")"
              )
            },
            character(1)
          ),
          collapse = " + "
        ),
        ")"
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Addition
  # ----------------------------------------------------------
  
  if (expr$op == "add") {
    
    left <- recurse(expr$args[[1]], index)
    right <- recurse(expr$args[[2]], index)
    
    if (!identical(dim(left), dim(right)) ||
        length(left) != length(right)) {
      stop("Incompatible dimensions in addition.")
    }
    
    return(paste0("(", left, " + ", right, ")"))
  }
  
  
  # ----------------------------------------------------------
  # Scalar multiplication
  # ----------------------------------------------------------
  
  if (expr$op == "multiply") {
    
    left <- recurse(expr$args[[1]], index)
    right <- recurse(expr$args[[2]], index)
    
    if (length(left) != 1 && length(right) != 1) {
      stop("multiply requires at least one scalar operand.")
    }
    
    return(paste0("(", left, " * ", right, ")"))
  }
  
  
  stop("Unknown expression operator: ", expr$op)
}