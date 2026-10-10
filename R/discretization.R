# R/discretization.R
#
# Supported operators:
#   field, grad, div, contract, matmul, add, multiply
#
# Convention:
#   grad(u)[i,j] = d u[i] / d X[j]
#   contract(C,G)[i,j] = sum_{k,l} C[i,j,k,l] * G[k,l]
#
# Discretization:
#   - uniform Cartesian grid
#   - divergence from shared face fluxes
#   - normal derivatives from neighboring node differences
#   - tangential face derivatives from averaged node derivatives
#   - arithmetic interpolation of coefficients to faces
#   - zero exterior flux when grid_bounds is supplied
#   - half-width control volumes at boundary nodes
#
# Output:
#   scalar: character string
#   vector: character vector
#   higher-order tensor: character array
#
# Material coefficients remain symbolic unless field_symbol
# explicitly replaces them by numerical values.


# ------------------------------------------------------------
# 1. Naming and grid utilities
# ------------------------------------------------------------

state_name <- function(variable, index) {
  paste(c(variable, as.integer(index)), collapse = "_")
}

neighbor_index <- function(index, dimension, direction) {
  result <- index
  result[dimension] <- result[dimension] + direction
  result
}

default_field_symbol <- function(
    name,
    index,
    shape,
    component = integer(0),
    state_variables = c("u", "v"),
    fields = NULL,
    parameters = list()
) {
  # Examples:
  #   scalar u:       u_3_4
  #   vector u[1]:    u_1_3_4
  #   E[1,2,1,2]:     E_1_2_1_2_3_4
  paste(
    c(name, as.integer(component), as.integer(index)),
    collapse = "_"
  )
}


# ------------------------------------------------------------
# 2. First-derivative stencils at nodes
# ------------------------------------------------------------

derivative_stencil <- function(
    index,
    dimension,
    grid_bounds = NULL
) {
  if (!is.null(grid_bounds)) {
    lower <- grid_bounds$min[dimension]
    upper <- grid_bounds$max[dimension]
    position <- index[dimension]
    
    if (position < lower || position > upper) {
      stop("Index lies outside grid_bounds.")
    }
    
    if (upper <= lower) {
      stop("Each axis needs at least two nodes.")
    }
    
    if (position == lower) {
      return(list(
        points = list(
          index,
          neighbor_index(index, dimension, 1L)
        ),
        coefficients = c(-1, 1),
        denominator = 1
      ))
    }
    
    if (position == upper) {
      return(list(
        points = list(
          neighbor_index(index, dimension, -1L),
          index
        ),
        coefficients = c(-1, 1),
        denominator = 1
      ))
    }
  }
  
  list(
    points = list(
      neighbor_index(index, dimension, -1L),
      neighbor_index(index, dimension, 1L)
    ),
    coefficients = c(-1, 1),
    denominator = 2
  )
}


# ------------------------------------------------------------
# 3. Recursive symbolic discretization
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
  
  if (!is.numeric(index) ||
      length(index) < 1L ||
      any(!is.finite(index)) ||
      any(index != floor(index))) {
    stop("index must contain finite integers.")
  }
  
  dimension <- length(index)
  
  if (!is.numeric(spacing) ||
      length(spacing) != dimension ||
      any(!is.finite(spacing)) ||
      any(spacing <= 0)) {
    stop("spacing must contain one positive finite value per axis.")
  }
  
  if (!is.function(field_symbol)) {
    stop("field_symbol must be a function.")
  }
  
  if (!is.null(grid_bounds)) {
    lower <- grid_bounds$min
    upper <- grid_bounds$max
    
    if (!is.numeric(lower) ||
        !is.numeric(upper) ||
        length(lower) != dimension ||
        length(upper) != dimension ||
        any(!is.finite(lower)) ||
        any(!is.finite(upper)) ||
        any(lower != floor(lower)) ||
        any(upper != floor(upper)) ||
        any(lower >= upper)) {
      stop("Invalid grid bounds: need at least two nodes per axis.")
    }
    
    if (any(index < lower) || any(index > upper)) {
      stop("index lies outside grid_bounds.")
    }
  }
  
  # ----------------------------------------------------------
  # Symbolic arithmetic
  # ----------------------------------------------------------
  
  number <- function(x) {
    sprintf("%.17g", x)
  }
  
  sum_strings <- function(x) {
    x <- x[x != "0"]
    
    if (length(x) == 0L) return("0")
    if (length(x) == 1L) return(x)
    
    paste0("(", paste(x, collapse = " + "), ")")
  }
  
  product <- function(a, b) {
    if (a == "0" || b == "0") return("0")
    if (a == "1") return(b)
    if (b == "1") return(a)
    
    paste0("(", a, ") * (", b, ")")
  }
  
  difference <- function(a, b) {
    if (identical(a, b)) return("0")
    
    paste0("((", a, ") - (", b, "))")
  }
  
  quotient <- function(a, denominator) {
    if (a == "0") return("0")
    
    paste0("(", a, ") / (", number(denominator), ")")
  }
  
  average <- function(a, b) {
    if (identical(a, b)) return(a)
    
    quotient(sum_strings(c(a, b)), 2)
  }
  
  # ----------------------------------------------------------
  # Evaluate one field component at a node
  # ----------------------------------------------------------
  
  field_at_node <- function(e, idx, component) {
    value <- field_symbol(
      name = e$name,
      index = idx,
      shape = e$shape,
      component = component,
      state_variables = state_variables,
      fields = fields,
      parameters = parameters
    )
    
    if (!is.character(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !nzchar(value)) {
      stop("field_symbol must return one nonempty expression string.")
    }
    
    value
  }
  
  # ----------------------------------------------------------
  # Differentiate one component at a node
  # ----------------------------------------------------------
  
  node_derivative <- function(e, idx, component, axis) {
    stencil <- derivative_stencil(
      index = idx,
      dimension = axis,
      grid_bounds = grid_bounds
    )
    
    terms <- vapply(seq_along(stencil$points), function(m) {
      product(
        number(stencil$coefficients[m]),
        value_at(e, stencil$points[[m]], component)
      )
    }, character(1))
    
    quotient(
      sum_strings(terms),
      stencil$denominator * spacing[axis]
    )
  }
  
  # ----------------------------------------------------------
  # Evaluate one expression component
  #
  # face = NULL: evaluate at node idx
  # face = axis: evaluate halfway between idx and idx + e_axis
  # ----------------------------------------------------------
  
  value_at <- function(e, idx, component = integer(0),
                       face = NULL) {
    
    # Field
    if (e$op == "field") {
      left <- field_at_node(e, idx, component)
      
      if (is.null(face)) {
        return(left)
      }
      
      right <- field_at_node(
        e,
        neighbor_index(idx, face, 1L),
        component
      )
      
      return(average(left, right))
    }
    
    # Addition
    if (e$op == "add") {
      return(sum_strings(c(
        value_at(e$args[[1]], idx, component, face),
        value_at(e$args[[2]], idx, component, face)
      )))
    }
    
    # Scalar multiplication
    if (e$op == "multiply") {
      a <- e$args[[1]]
      b <- e$args[[2]]
      
      a_scalar <- length(a$shape) == 0L
      b_scalar <- length(b$shape) == 0L
      
      if (!a_scalar && !b_scalar) {
        stop("multiply requires at least one scalar.")
      }
      
      return(product(
        value_at(
          a, idx,
          if (a_scalar) integer(0) else component,
          face
        ),
        value_at(
          b, idx,
          if (b_scalar) integer(0) else component,
          face
        )
      ))
    }
    
    # Gradient
    if (e$op == "grad") {
      if (tail(e$shape, 1L) != dimension) {
        stop("Gradient dimension does not match grid dimension.")
      }
      
      child <- e$args[[1]]
      axis <- tail(component, 1L)
      child_component <- head(component, -1L)
      
      if (is.null(face)) {
        return(node_derivative(
          child, idx, child_component, axis
        ))
      }
      
      right_idx <- neighbor_index(idx, face, 1L)
      
      # Normal derivative at a face.
      if (axis == face) {
        return(quotient(
          difference(
            value_at(child, right_idx, child_component),
            value_at(child, idx, child_component)
          ),
          spacing[axis]
        ))
      }
      
      # Tangential derivative at a face.
      # Uses immediate and diagonal neighboring nodes.
      return(average(
        node_derivative(child, idx, child_component, axis),
        node_derivative(child, right_idx, child_component, axis)
      ))
    }
    
    # Fourth-order tensor contracted with a matrix.
    if (e$op == "contract") {
      A <- e$args[[1]]
      B <- e$args[[2]]
      
      if (length(A$shape) != 4L ||
          length(B$shape) != 2L ||
          any(A$shape[3:4] != B$shape)) {
        stop("Invalid shapes for contract.")
      }
      
      terms <- character(0)
      
      for (k in seq_len(B$shape[1])) {
        for (l in seq_len(B$shape[2])) {
          terms <- c(terms, product(
            value_at(A, idx, c(component, k, l), face),
            value_at(B, idx, c(k, l), face)
          ))
        }
      }
      
      return(sum_strings(terms))
    }
    
    # Matrix-vector or matrix-matrix multiplication.
    if (e$op == "matmul") {
      A <- e$args[[1]]
      B <- e$args[[2]]
      
      if (length(A$shape) != 2L ||
          !(length(B$shape) %in% c(1L, 2L)) ||
          A$shape[2] != B$shape[1]) {
        stop("Invalid shapes for matmul.")
      }
      
      terms <- vapply(seq_len(A$shape[2]), function(k) {
        b_component <- if (length(B$shape) == 1L) {
          k
        } else {
          c(k, component[2])
        }
        
        product(
          value_at(A, idx, c(component[1], k), face),
          value_at(B, idx, b_component, face)
        )
      }, character(1))
      
      return(sum_strings(terms))
    }
    
    # Conservative divergence.
    if (e$op == "div") {
      if (!is.null(face)) {
        stop("A divergence inside a face expression is unsupported.")
      }
      
      child <- e$args[[1]]
      
      if (length(child$shape) == 0L ||
          tail(child$shape, 1L) != dimension) {
        stop("Last flux index must match spatial dimension.")
      }
      
      terms <- vapply(seq_len(dimension), function(axis) {
        at_lower <- !is.null(grid_bounds) &&
          idx[axis] == grid_bounds$min[axis]
        
        at_upper <- !is.null(grid_bounds) &&
          idx[axis] == grid_bounds$max[axis]
        
        flux_component <- c(component, axis)
        
        flux_plus <- if (at_upper) {
          "0"
        } else {
          value_at(child, idx, flux_component, face = axis)
        }
        
        flux_minus <- if (at_lower) {
          "0"
        } else {
          value_at(
            child,
            neighbor_index(idx, axis, -1L),
            flux_component,
            face = axis
          )
        }
        
        width <- spacing[axis]
        
        if (at_lower || at_upper) {
          width <- width / 2
        }
        
        quotient(
          difference(flux_plus, flux_minus),
          width
        )
      }, character(1))
      
      return(sum_strings(terms))
    }
    
    stop("Unknown expression operator: ", e$op)
  }
  
  # ----------------------------------------------------------
  # Assemble all output components in R array order
  # ----------------------------------------------------------
  
  shape <- expr$shape
  
  if (length(shape) == 0L) {
    return(value_at(expr, index))
  }
  
  components <- as.matrix(expand.grid(
    lapply(shape, seq_len),
    KEEP.OUT.ATTRS = FALSE
  ))
  
  result <- vapply(seq_len(nrow(components)), function(m) {
    value_at(
      expr,
      index,
      as.integer(components[m, ])
    )
  }, character(1))
  
  if (length(shape) == 1L) {
    return(unname(result))
  }
  
  array(result, dim = shape)
}



# 
# ########################################   Testing   #####################################
# 
# 
# source("R/expression_tree.R")
# source("R/operators.R")
# source("R/discretization.R")
# 
# u <- field("u", c(2))
# v <- field("v", c(2))
# E <- field("E", c(2, 2, 2, 2))
# D <- field("D", c(2, 2, 2, 2))
# 
# force <- div(add(
#   contract(E, grad(u)),
#   contract(D, grad(v))
# ))
# 
# bounds <- list(min = c(1, 1), max = c(7, 7))
# h <- c(0.5, 0.8)
# 
# # Anisotropic material:
# # stress = C %*% c(eps_xx, eps_yy, 2*eps_xy)
# C <- matrix(c(
#   4,   1,   0.5,
#   1,   3,   0.2,
#   0.5, 0.2, 2
# ), 3, 3, byrow = TRUE)
# 
# # Convert to fourth-order tensor.
# voigt <- matrix(c(1, 3, 3, 2), 2, 2)
# Ct <- array(0, c(2, 2, 2, 2))
# 
# for (i in 1:2) for (j in 1:2) {
#   for (k in 1:2) for (l in 1:2) {
#     Ct[i, j, k, l] <- C[voigt[i, j], voigt[k, l]]
#   }
# }
# 
# # Insert numerical material coefficients, retain state symbols.
# symbols <- function(name, index, shape, component,
#                     state_variables, fields, parameters) {
#   if (name %in% c("E", "D")) {
#     coefficient <- Ct[matrix(component, nrow = 1)]
#     if (name == "D") coefficient <- 0.25 * coefficient
#     return(sprintf("%.17g", coefficient))
#   }
#   
#   default_field_symbol(
#     name, index, shape, component,
#     state_variables, fields, parameters
#   )
# }
# 
# equations_at <- function(index) {
#   discretize(
#     force,
#     index = index,
#     spacing = h,
#     grid_bounds = bounds,
#     field_symbol = symbols
#   )
# }
# 
# zero <- function(x, y) c(0, 0)
# 
# evaluate <- function(rhs, displacement, velocity = zero) {
#   values <- new.env(parent = baseenv())
#   
#   for (i in 1:7) for (j in 1:7) {
#     x <- (i - 1) * h[1]
#     y <- (j - 1) * h[2]
#     
#     U <- displacement(x, y)
#     V <- velocity(x, y)
#     
#     for (a in 1:2) {
#       values[[paste("u", a, i, j, sep = "_")]] <- U[a]
#       values[[paste("v", a, i, j, sep = "_")]] <- V[a]
#     }
#   }
#   
#   vapply(rhs, function(s) {
#     eval(parse(text = s), envir = values)
#   }, numeric(1))
# }
# 
# check <- function(label, actual, expected) {
#   stopifnot(isTRUE(all.equal(
#     unname(actual), expected, tolerance = 1e-10
#   )))
#   cat("OK:", label, "\n")
# }
# 
# rhs <- equations_at(c(4, 4))
# 
# check(
#   "Quadratic displacement",
#   evaluate(rhs, function(x, y) c(x^2, y^2)),
#   c(8.4, 7)
# )
# 
# check(
#   "Mixed derivatives",
#   evaluate(rhs, function(x, y) c(x * y, 0)),
#   c(1, 3)
# )
# 
# check(
#   "Viscous contribution",
#   evaluate(rhs, zero, function(x, y) c(x^2, y^2)),
#   0.25 * c(8.4, 7)
# )
# 
# check(
#   "Constant displacement at free corner",
#   evaluate(equations_at(c(7, 7)), function(x, y) c(1, 2)),
#   c(0, 0)
# )
# 
# # Inspect the generated first force component:
# cat("Force component 1:\n", rhs[1], "\n")
# 
# 
# 
# 
# # Test 1
# # Spatially varying material factors
# a <- function(x, y) 1 + x + 0.5 * y
# b <- function(x, y) 2 + 0.3 * x - 0.2 * y
# 
# spatial_symbols <- function(
#     name, index, shape, component,
#     state_variables, fields, parameters
# ) {
#   if (name %in% c("E", "D")) {
#     X <- (index - bounds$min) * h
#     
#     factor <- if (name == "E") {
#       a(X[1], X[2])
#     } else {
#       b(X[1], X[2])
#     }
#     
#     coefficient <- factor * Ct[matrix(component, nrow = 1)]
#     return(sprintf("%.17g", coefficient))
#   }
#   
#   default_field_symbol(
#     name, index, shape, component,
#     state_variables, fields, parameters
#   )
# }
# 
# idx <- c(4, 4)
# X <- (idx - bounds$min) * h
# x <- X[1]
# y <- X[2]
# 
# rhs_spatial <- discretize(
#   force,
#   index = idx,
#   spacing = h,
#   grid_bounds = bounds,
#   field_symbol = spatial_symbols
# )
# 
# # Analytical stress for w = (x^2, y^2), using the test tensor C:
# sigma <- matrix(c(
#   8*x + 2*y,  x + 0.4*y,
#   x + 0.4*y,  2*x + 6*y
# ), nrow = 2, byrow = TRUE)
# 
# # div(f * sigma) = f * div(sigma) + sigma %*% grad(f)
# div_sigma <- c(8.4, 7)
# 
# expected_E <- a(x, y) * div_sigma +
#   drop(sigma %*% c(1, 0.5))
# 
# expected_D <- b(x, y) * div_sigma +
#   drop(sigma %*% c(0.3, -0.2))
# 
# quadratic <- function(x, y) c(x^2, y^2)
# 
# check(
#   "Spatially varying E",
#   evaluate(rhs_spatial, quadratic),
#   expected_E
# )
# 
# check(
#   "Spatially varying D",
#   evaluate(rhs_spatial, zero, quadratic),
#   expected_D
# )
# 
# 
# 
# 
# 
# #2
# # Uses equations_at(), evaluate(), check() and h
# # from the original test script with constant material.
# 
# affine <- function(x, y) c(x, 0)
# 
# sigma <- matrix(c(
#   4,   0.5,
#   0.5, 1
# ), nrow = 2, byrow = TRUE)
# 
# check(
#   "Constant stress: zero divergence inside",
#   evaluate(equations_at(c(4, 4)), affine),
#   c(0, 0)
# )
# 
# check(
#   "Right free boundary",
#   evaluate(equations_at(c(7, 4)), affine),
#   -sigma[, 1] / (h[1] / 2)
# )
# 
# check(
#   "Bottom free boundary",
#   evaluate(equations_at(c(4, 1)), affine),
#   sigma[, 2] / (h[2] / 2)
# )
# 
# check(
#   "Top-right free corner",
#   evaluate(equations_at(c(7, 7)), affine),
#   -sigma[, 1] / (h[1] / 2) -
#     sigma[, 2] / (h[2] / 2)
# )
# 
# 
# 
# 
# 
# 
# 
# 
# # 3
# # Requires Ct, force and symbols from the original test script.
# 
# convergence_test <- function() {
#   X <- c(0.5, 0.5)
#   x <- X[1]
#   y <- X[2]
#   
#   # H[k,j,l] = second derivative of u[k] wrt X[j], X[l]
#   H <- array(0, c(2, 2, 2))
#   
#   H[1, , ] <- matrix(c(
#     -sin(x)*cos(y), -cos(x)*sin(y),
#     -cos(x)*sin(y), -sin(x)*cos(y)
#   ), 2, 2)
#   
#   H[2, , ] <- matrix(c(
#     -cos(x)*sin(y), -sin(x)*cos(y),
#     -sin(x)*cos(y), -cos(x)*sin(y)
#   ), 2, 2)
#   
#   # Exact div(C : grad(u)) for constant C.
#   exact <- numeric(2)
#   
#   for (i in 1:2) for (j in 1:2) {
#     for (k in 1:2) for (l in 1:2) {
#       exact[i] <- exact[i] + Ct[i, j, k, l] * H[k, j, l]
#     }
#   }
#   
#   sizes <- c(9L, 17L, 33L, 65L)
#   
#   errors <- vapply(sizes, function(n) {
#     step <- 1 / (n - 1)
#     idx <- rep(as.integer((n + 1) / 2), 2)
#     
#     rhs <- discretize(
#       force,
#       index = idx,
#       spacing = rep(step, 2),
#       grid_bounds = list(min = c(1, 1), max = c(n, n)),
#       field_symbol = symbols
#     )
#     
#     values <- new.env(parent = baseenv())
#     
#     for (i in seq_len(n)) for (j in seq_len(n)) {
#       xx <- (i - 1) * step
#       yy <- (j - 1) * step
#       
#       U <- c(sin(xx)*cos(yy), cos(xx)*sin(yy))
#       
#       for (a in 1:2) {
#         values[[paste("u", a, i, j, sep = "_")]] <- U[a]
#         values[[paste("v", a, i, j, sep = "_")]] <- 0
#       }
#     }
#     
#     numerical <- vapply(rhs, function(s) {
#       eval(parse(text = s), envir = values)
#     }, numeric(1))
#     
#     sqrt(sum((numerical - exact)^2))
#   }, numeric(1))
#   
#   orders <- log2(head(errors, -1) / tail(errors, -1))
#   
#   print(data.frame(
#     nodes_per_axis = sizes,
#     spacing = 1 / (sizes - 1),
#     error = errors,
#     order = c(NA, orders)
#   ))
#   
#   stopifnot(
#     all(diff(errors) < 0),
#     all(abs(orders - 2) < 0.15)
#   )
#   
#   cat("OK: second-order convergence in the interior\n")
# }
# 
# convergence_test()
# 
# 
# 
# 
# 
# 
# # Requires u, v, E, D, symbols from the original test script.
# # Uses its constant material and D = 0.25 * E.
# 
# mechanics_test <- function() {
#   n <- 5L
#   h <- c(0.5, 0.8)
#   bounds <- list(min = c(1, 1), max = c(n, n))
#   points <- as.matrix(expand.grid(i = 1:n, j = 1:n))
#   ndof <- 2L * nrow(points)
#   
#   # Component order: (u1, u2) at each node.
#   state_names <- function(variable) {
#     unlist(lapply(seq_len(nrow(points)), function(p) {
#       vapply(1:2, function(a) {
#         paste(c(variable, a, points[p, ]), collapse = "_")
#       }, character(1))
#     }), use.names = FALSE)
#   }
#   
#   # Extract the linear operator by applying it to basis vectors.
#   operator_matrix <- function(expression, variable) {
#     expressions <- unlist(
#       lapply(seq_len(nrow(points)), function(p) {
#         discretize(
#           expression,
#           index = points[p, ],
#           spacing = h,
#           grid_bounds = bounds,
#           field_symbol = symbols
#         )
#       }),
#       use.names = FALSE
#     )
#     
#     parsed <- lapply(expressions, function(s) parse(text = s)[[1]])
#     names <- state_names(variable)
#     env <- list2env(
#       setNames(as.list(rep(0, ndof)), names),
#       parent = baseenv()
#     )
#     
#     L <- matrix(0, ndof, ndof)
#     
#     for (j in seq_len(ndof)) {
#       env[[names[j]]] <- 1
#       
#       L[, j] <- vapply(parsed, function(e) {
#         eval(e, envir = env)
#       }, numeric(1))
#       
#       env[[names[j]]] <- 0
#     }
#     
#     L
#   }
#   
#   K <- operator_matrix(div(contract(E, grad(u))), "u")
#   B <- operator_matrix(div(contract(D, grad(v))), "v")
#   
#   # Control-volume areas: half at edges, quarter at corners.
#   volumes <- apply(points, 1, function(idx) {
#     widths <- h
#     widths[idx == 1 | idx == n] <-
#       widths[idx == 1 | idx == n] / 2
#     prod(widths)
#   })
#   
#   weights <- rep(volumes, each = 2)
#   
#   # Volume-weighted force operators.
#   WK <- sweep(K, 1, weights, "*")
#   WB <- sweep(B, 1, weights, "*")
#   
#   # Infinitesimal rigid rotation: u = (-y, x).
#   rotation <- unlist(
#     lapply(seq_len(nrow(points)), function(p) {
#       X <- (points[p, ] - 1) * h
#       c(-X[2], X[1])
#     }),
#     use.names = FALSE
#   )
#   
#   max_eigenvalue <- function(A) {
#     max(eigen(
#       (A + t(A)) / 2,
#       symmetric = TRUE,
#       only.values = TRUE
#     )$values)
#   }
#   
#   scale_K <- max(1, norm(WK, "F"))
#   scale_B <- max(1, norm(WB, "F"))
#   
#   rotation_error <- max(abs(K %*% rotation)) /
#     max(1, norm(K, "I") * max(abs(rotation)))
#   
#   symmetry_error <- norm(WK - t(WK), "F") / scale_K
#   elastic_eigenvalue <- max_eigenvalue(WK) / scale_K
#   damping_eigenvalue <- max_eigenvalue(WB) / scale_B
#   
#   tol <- 1e-10
#   
#   results <- data.frame(
#     test = c(
#       "Rigid rotation: zero force",
#       "Elasticity: weighted symmetry",
#       "Elasticity: nonpositive eigenvalues",
#       "Damping: nonpositive power"
#     ),
#     value = c(
#       rotation_error,
#       symmetry_error,
#       elastic_eigenvalue,
#       damping_eigenvalue
#     ),
#     passed = c(
#       rotation_error <= tol,
#       symmetry_error <= tol,
#       elastic_eigenvalue <= tol,
#       damping_eigenvalue <= tol
#     )
#   )
#   
#   print(results, row.names = FALSE)
#   invisible(results)
# }
# 
# mechanics_test()
# 
# 
# 
# #4



