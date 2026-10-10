# ============================================================
# Expression tree: nodes and tensor shapes
# ============================================================

new_expr <- function(op, args = list(),
                     shape = integer(0), name = NULL) {
  
  structure(
    list(
      op    = op,
      args  = args,
      shape = shape,
      name  = name
    ),
    class = "math_expr"
  )
}


# Create a field
#
# shape = integer(0): scalar
# shape = c(2):       vector in 2D
# shape = c(2, 2):    matrix in 2D

field <- function(name, shape = integer(0)) {
  
  stopifnot(
    is.character(name),
    length(name) == 1,
    length(shape) >= 0,
    all(shape > 0)
  )
  
  new_expr(
    op    = "field",
    name  = name,
    shape = shape
  )
}


# Display expression type

shape_name <- function(shape) {
  
  if (length(shape) == 0) {
    return("scalar")
  }
  
  paste(shape, collapse = " x ")
}


show_type <- function(x) {
  
  stopifnot(inherits(x, "math_expr"))
  
  cat(
    "Operation:", x$op,
    "\nName:", if (is.null(x$name)) "-" else x$name,
    "\nShape:", shape_name(x$shape),
    "\n"
  )
}



# Recursively print the expression tree
print_expr <- function(x, indent = 0) {
  
  stopifnot(inherits(x, "math_expr"))
  
  prefix <- paste0(strrep("  ", indent))
  
  label <- if (x$op == "field") {
    paste0("Field: ", x$name)
  } else {
    paste0("Operator: ", x$op)
  }
  
  cat(
    prefix, label,
    " [", shape_name(x$shape), "]\n",
    sep = ""
  )
  
  # Recursively print all child expressions
  for (arg in x$args) {
    print_expr(arg, indent = indent + 1)
  }
  
  invisible(x)
}