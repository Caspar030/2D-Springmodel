# ============================================================
# Grid geometry
# ============================================================

create_grid <- function(coords) {
  
  if (!is.list(coords) || length(coords) < 1 ||
      length(coords) > 3) {
    stop("coords must be a list of 1 to 3 coordinate vectors.")
  }
  
  if (is.null(names(coords)) || any(names(coords) == "")) {
    stop("Name the coordinate vectors, e.g. x, y, z.")
  }
  
  for (axis in coords) {
    if (!is.numeric(axis) || length(axis) < 1 ||
        any(!is.finite(axis))) {
      stop("Coordinates must be finite numeric vectors.")
    }
    
    if (length(axis) > 1 && any(diff(axis) <= 0)) {
      stop("Coordinates must be strictly increasing.")
    }
  }
  
  structure(
    list(
      dim   = length(coords),
      axes  = names(coords),
      coords = coords,
      N     = lengths(coords),
      nodes = expand.grid(
        coords,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    ),
    class = "model_grid"
  )
}


# Get the coordinates of a node by its grid indices.
# Example: node_position(grid, c(2, 3))

node_position <- function(grid, index) {
  
  stopifnot(inherits(grid, "model_grid"))
  
  if (length(index) != grid$dim ||
      any(index != as.integer(index)) ||
      any(index < 1) ||
      any(index > grid$N)) {
    stop("Invalid grid index.")
  }
  
  position <- vapply(
    seq_len(grid$dim),
    function(k) grid$coords[[k]][index[k]],
    numeric(1)
  )
  
  names(position) <- grid$axes
  position
}


# Generate a state name from grid indices.
# Example: state_name("u", c(2, 3)) -> "u_2_3"

state_name <- function(variable, index) {
  
  paste(
    c(variable, index),
    collapse = "_"
  )
}
