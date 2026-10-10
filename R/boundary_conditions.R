# R/boundary_conditions.R
#
# Boundary utilities for a two-dimensional grid.
#
# index[1] = x-direction
# index[2] = y-direction
#
# Default:
#   left boundary: u = v = (0, 0)
#   other boundaries: traction-free
#
# At corners, the fixed boundary takes precedence.
#
# This file classifies nodes and generates symbolic constraints.
# Enforcement is handled by discretization and the ODE system.
# No source() calls or test code are executed here.


# ------------------------------------------------------------
# 1. Validate grid bounds
# ------------------------------------------------------------

validate_grid_bounds <- function(grid_bounds) {
  
  if (!is.list(grid_bounds)) {
    stop("grid_bounds must be a list with min and max.")
  }
  
  lower <- grid_bounds$min
  upper <- grid_bounds$max
  
  valid_bound <- function(x) {
    is.numeric(x) &&
      is.null(dim(x)) &&
      length(x) == 2L &&
      all(is.finite(x)) &&
      all(x == floor(x))
  }
  
  if (!valid_bound(lower) || !valid_bound(upper)) {
    stop("min and max must contain two finite integer indices.")
  }
  
  if (any(lower >= upper)) {
    stop("Each axis must contain at least two grid points.")
  }
  
  invisible(TRUE)
}


# ------------------------------------------------------------
# 2. Validate a grid index
# ------------------------------------------------------------

validate_grid_index <- function(index, grid_bounds) {
  
  validate_grid_bounds(grid_bounds)
  
  if (!is.numeric(index) ||
      !is.null(dim(index)) ||
      length(index) != 2L ||
      any(!is.finite(index)) ||
      any(index != floor(index))) {
    stop("index must contain two finite integers.")
  }
  
  if (any(index < grid_bounds$min) ||
      any(index > grid_bounds$max)) {
    stop("index lies outside grid_bounds.")
  }
  
  invisible(TRUE)
}


# ------------------------------------------------------------
# 3. Identify boundary sides
# ------------------------------------------------------------

boundary_sides <- function(index, grid_bounds) {
  
  validate_grid_index(index, grid_bounds)
  
  c(
    left   = index[1] == grid_bounds$min[1],
    right  = index[1] == grid_bounds$max[1],
    bottom = index[2] == grid_bounds$min[2],
    top    = index[2] == grid_bounds$max[2]
  )
}


# ------------------------------------------------------------
# 4. Classify a node
# ------------------------------------------------------------

classify_boundary <- function(
    index,
    grid_bounds,
    fixed_side = "left"
) {
  
  fixed_side <- match.arg(
    fixed_side,
    c("left", "right", "bottom", "top")
  )
  
  sides <- boundary_sides(index, grid_bounds)
  
  # Fixed conditions also apply at the two adjoining corners.
  if (sides[[fixed_side]]) {
    return("fixed")
  }
  
  if (any(sides)) {
    return("free")
  }
  
  "interior"
}


# ------------------------------------------------------------
# 5. Boundary predicates
# ------------------------------------------------------------

is_fixed_boundary <- function(
    index,
    grid_bounds,
    fixed_side = "left"
) {
  classify_boundary(index, grid_bounds, fixed_side) == "fixed"
}

is_free_boundary <- function(
    index,
    grid_bounds,
    fixed_side = "left"
) {
  classify_boundary(index, grid_bounds, fixed_side) == "free"
}


# ------------------------------------------------------------
# 6. Names of constrained vector components
#
# Example:
#   u, component 1, node (1,3) -> u_1_1_3
#
# Matches default_field_symbol() in discretization.R.
# ------------------------------------------------------------

fixed_boundary_states <- function(
    index,
    grid_bounds,
    state_variables = c("u", "v"),
    fixed_side = "left"
) {
  
  if (!is.character(state_variables) ||
      length(state_variables) == 0L ||
      anyNA(state_variables) ||
      any(!nzchar(state_variables)) ||
      anyDuplicated(state_variables) > 0L) {
    stop("state_variables must contain unique, nonempty names.")
  }
  
  if (!is_fixed_boundary(index, grid_bounds, fixed_side)) {
    return(character(0))
  }
  
  unlist(
    lapply(state_variables, function(variable) {
      vapply(seq_len(2L), function(component) {
        paste(
          c(variable, component, index),
          collapse = "_"
        )
      }, character(1))
    }),
    use.names = FALSE
  )
}


# ------------------------------------------------------------
# 7. Symbolic fixed-boundary constraints
# ------------------------------------------------------------

fixed_boundary_equations <- function(
    index,
    grid_bounds,
    state_variables = c("u", "v"),
    fixed_side = "left"
) {
  
  states <- fixed_boundary_states(
    index = index,
    grid_bounds = grid_bounds,
    state_variables = state_variables,
    fixed_side = fixed_side
  )
  
  if (length(states) == 0L) {
    return(character(0))
  }
  
  paste0(states, " = 0")
}


# ------------------------------------------------------------
# 8. Enumerate all grid points
# First coordinate varies fastest.
# ------------------------------------------------------------

grid_points <- function(grid_bounds) {
  
  validate_grid_bounds(grid_bounds)
  
  points <- expand.grid(
    x = seq(
      from = grid_bounds$min[1],
      to = grid_bounds$max[1],
      by = 1
    ),
    y = seq(
      from = grid_bounds$min[2],
      to = grid_bounds$max[2],
      by = 1
    ),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  as.matrix(points)
}


# ------------------------------------------------------------
# 9. Classify the entire grid
# x and y below are grid indices, not physical coordinates.
# ------------------------------------------------------------

classify_grid <- function(
    grid_bounds,
    fixed_side = "left"
) {
  
  points <- grid_points(grid_bounds)
  
  classifications <- vapply(
    seq_len(nrow(points)),
    function(p) {
      classify_boundary(
        index = points[p, ],
        grid_bounds = grid_bounds,
        fixed_side = fixed_side
      )
    },
    character(1)
  )
  
  data.frame(
    index = seq_len(nrow(points)),
    x = points[, 1],
    y = points[, 2],
    boundary = classifications,
    stringsAsFactors = FALSE
  )
}