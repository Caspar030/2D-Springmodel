# R/ode_system.R
#
# Build a symbolic ODE system from expression trees.
#
# Dependencies:
#   expression_tree.R
#   discretization.R
#   boundary_conditions.R
#
# Each named expression specifies a time derivative:
#   expressions$u -> du/dt
#   expressions$v -> dv/dt
#
# Supported states:
#   scalar fields or two-component vector fields
#
# Homogeneous fixed boundary:
#   - all supplied states are zero on fixed_side
#   - fixed states are omitted from the ODE system
#   - references to fixed states are replaced by literal zero
#
# Output:
#   data.frame with one row per dynamic state
#   rhs contains an R-compatible expression string
#
# Parameters remain symbolic by default.
# No numerical integration or dMod2 conversion is performed here.


build_ode_system <- function(
    expressions,
    grid_bounds,
    spacing,
    fields = NULL,
    parameters = list(),
    field_symbol = default_field_symbol,
    state_variables = c("u", "v"),
    fixed_side = "left"
) {
  
  # ----------------------------------------------------------
  # 1. Validate inputs
  # ----------------------------------------------------------
  
  validate_grid_bounds(grid_bounds)
  
  fixed_side <- match.arg(
    fixed_side,
    c("left", "right", "bottom", "top")
  )
  
  if (!is.numeric(spacing) ||
      length(spacing) != 2L ||
      any(!is.finite(spacing)) ||
      any(spacing <= 0)) {
    stop("spacing must contain two positive finite values.")
  }
  
  if (!is.function(field_symbol)) {
    stop("field_symbol must be a function.")
  }
  
  # A single expression is interpreted as du/dt.
  if (inherits(expressions, "math_expr")) {
    expressions <- list(u = expressions)
  }
  
  if (!is.list(expressions) || length(expressions) == 0L) {
    stop("expressions must be a nonempty named list.")
  }
  
  variables <- names(expressions)
  
  if (is.null(variables) ||
      anyNA(variables) ||
      any(!nzchar(variables)) ||
      anyDuplicated(variables) > 0L) {
    stop("Expression names must be nonempty and unique.")
  }
  
  if (!all(vapply(
    expressions, inherits, logical(1), "math_expr"
  ))) {
    stop("Every expression must be a math_expr.")
  }
  
  if (!is.character(state_variables) ||
      anyNA(state_variables) ||
      anyDuplicated(state_variables) > 0L ||
      !setequal(state_variables, variables)) {
    stop(
      "state_variables must match names(expressions). ",
      "Every dynamic field needs its own equation."
    )
  }
  
  # Each RHS has the same shape as its associated state.
  shapes <- lapply(expressions, function(e) e$shape)
  
  for (variable in variables) {
    shape <- shapes[[variable]]
    
    if (!(length(shape) == 0L ||
          (length(shape) == 1L && shape[1] == 2L))) {
      stop(
        "Equation for '", variable,
        "' must be scalar or a two-component vector."
      )
    }
  }
  
  # Check references to known states within the expression tree.
  validate_state_references <- function(e) {
    if (e$op == "field" && e$name %in% state_variables) {
      if (!identical(
        as.integer(e$shape),
        as.integer(shapes[[e$name]])
      )) {
        stop("Inconsistent shape for state field '", e$name, "'.")
      }
    }
    
    for (child in e$args) {
      validate_state_references(child)
    }
    
    invisible(TRUE)
  }
  
  for (e in expressions) {
    validate_state_references(e)
  }
  
  
  # ----------------------------------------------------------
  # 2. State naming
  #
  # Scalar: u_i_j
  # Vector: u_component_i_j
  # ----------------------------------------------------------
  
  state_symbol <- function(variable, index,
                           component = integer(0)) {
    paste(
      c(variable, as.integer(component), as.integer(index)),
      collapse = "_"
    )
  }
  
  
  # ----------------------------------------------------------
  # 3. Apply fixed values during symbolic discretization
  # ----------------------------------------------------------
  
  boundary_field_symbol <- function(
    name,
    index,
    shape,
    component = integer(0),
    state_variables = c("u", "v"),
    fields = NULL,
    parameters = list()
  ) {
    
    if (name %in% state_variables) {
      if (is_fixed_boundary(index, grid_bounds, fixed_side)) {
        return("0")
      }
      
      # Use one consistent naming convention for ODE states.
      return(state_symbol(name, index, component))
    }
    
    # Coefficients retain the supplied symbolic representation.
    field_symbol(
      name = name,
      index = index,
      shape = shape,
      component = component,
      state_variables = state_variables,
      fields = fields,
      parameters = parameters
    )
  }
  
  
  # ----------------------------------------------------------
  # 4. Generate dynamic equations
  # ----------------------------------------------------------
  
  points <- grid_points(grid_bounds)
  rows <- list()
  fixed_states <- character(0)
  
  for (p in seq_len(nrow(points))) {
    index <- as.integer(points[p, ])
    
    boundary_type <- classify_boundary(
      index,
      grid_bounds,
      fixed_side = fixed_side
    )
    
    for (variable in variables) {
      shape <- shapes[[variable]]
      is_scalar <- length(shape) == 0L
      n_components <- if (is_scalar) 1L else shape[1]
      
      names_at_node <- vapply(
        seq_len(n_components),
        function(a) {
          state_symbol(
            variable,
            index,
            component = if (is_scalar) integer(0) else a
          )
        },
        character(1)
      )
      
      # Eliminate fixed states rather than integrating du/dt = 0.
      if (boundary_type == "fixed") {
        fixed_states <- c(fixed_states, names_at_node)
        next
      }
      
      rhs <- discretize(
        expr = expressions[[variable]],
        index = index,
        spacing = spacing,
        fields = fields,
        parameters = parameters,
        field_symbol = boundary_field_symbol,
        state_variables = state_variables,
        grid_bounds = grid_bounds
      )
      
      if (!is.character(rhs) ||
          !is.null(dim(rhs)) ||
          length(rhs) != n_components ||
          anyNA(rhs) ||
          any(!nzchar(rhs))) {
        stop(
          "Unexpected discretization result for '", variable,
          "' at node (", paste(index, collapse = ","), ")."
        )
      }
      
      for (a in seq_len(n_components)) {
        # Check syntax without evaluating symbolic parameters.
        tryCatch(
          {
            parsed <- parse(text = rhs[a])
            
            if (length(parsed) != 1L) {
              stop("Expected exactly one expression.")
            }
          },
          error = function(e) {
            stop(
              "Invalid RHS for ", names_at_node[a], ": ",
              conditionMessage(e),
              call. = FALSE
            )
          }
        )
        
        rows[[length(rows) + 1L]] <- data.frame(
          index = paste(index, collapse = ","),
          variable = variable,
          component = if (is_scalar) NA_integer_ else a,
          boundary = boundary_type,
          state = names_at_node[a],
          rhs = unname(rhs[a]),
          equation = paste0(
            "d(", names_at_node[a], ")/dt = ", rhs[a]
          ),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  system <- do.call(rbind, rows)
  rownames(system) <- NULL
  
  if (anyDuplicated(system$state) > 0L) {
    stop("Generated state names are not unique.")
  }
  
  # Metadata for reconstruction and later model conversion.
  attr(system, "fixed_states") <- setNames(
    rep(0, length(fixed_states)),
    fixed_states
  )
  attr(system, "grid_bounds") <- grid_bounds
  attr(system, "spacing") <- spacing
  attr(system, "fixed_side") <- fixed_side
  attr(system, "state_shapes") <- shapes
  
  system
}






















#Testing



source("R/expression_tree.R")
source("R/operators.R")
source("R/discretization.R")
source("R/boundary_conditions.R")
source("R/ode_system.R")

# Reduced Kelvin–Voigt system:
# rho = 1, no active force or local retraction.
u <- field("u", c(2))
v <- field("v", c(2))
E <- field("E", c(2, 2, 2, 2))
D <- field("D", c(2, 2, 2, 2))

expressions <- list(
  u = v,
  v = div(add(
    contract(E, grad(u)),
    contract(D, grad(v))
  ))
)

bounds <- list(min = c(1, 1), max = c(7, 5))

make_system <- function(side) {
  build_ode_system(
    expressions = expressions,
    grid_bounds = bounds,
    spacing = c(0.5, 0.8),
    state_variables = c("u", "v"),
    fixed_side = side
  )
}

check <- function(label, condition) {
  if (!isTRUE(condition)) stop("FAILED: ", label)
  cat("OK:", label, "\n")
}

system <- make_system("left")

# 1. 35 nodes - 5 fixed nodes; 4 states per remaining node.
check(
  "Number of dynamic equations",
  nrow(system) == 120L &&
    sum(system$variable == "u") == 60L &&
    sum(system$variable == "v") == 60L
)

check(
  "Unique component names",
  anyDuplicated(system$state) == 0L &&
    all(system$component %in% 1:2) &&
    "u_1_2_3" %in% system$state &&
    "v_2_2_3" %in% system$state
)

check(
  "Boundary classification",
  sum(system$boundary == "free") == 60L &&
    sum(system$boundary == "interior") == 60L &&
    !any(system$boundary == "fixed")
)

# 2. du/dt must equal the corresponding velocity component.
u_rows <- system$variable == "u"

check(
  "Displacement-velocity coupling",
  all(
    system$rhs[u_rows] ==
      sub("^u_", "v_", system$state[u_rows])
  )
)

# 3. Fixed states must be absent from states AND right-hand sides.
fixed <- attr(system, "fixed_states")

expected_fixed <- unlist(lapply(1:5, function(j) {
  c(
    paste("u", 1, 1, j, sep = "_"),
    paste("u", 2, 1, j, sep = "_"),
    paste("v", 1, 1, j, sep = "_"),
    paste("v", 2, 1, j, sep = "_")
  )
}), use.names = FALSE)

check(
  "Correct fixed-state metadata",
  setequal(names(fixed), expected_fixed) &&
    all(fixed == 0)
)

# Parsing also checks that every RHS has valid R syntax.
parsed <- lapply(system$rhs, function(s) parse(text = s))
symbols_used <- unique(unlist(lapply(parsed, all.vars)))

check(
  "Fixed states eliminated everywhere",
  !any(expected_fixed %in% system$state) &&
    !any(expected_fixed %in% symbols_used)
)

# Every referenced u/v component needs a dynamic equation.
state_references <- grep(
  "^(u|v)_", symbols_used, value = TRUE
)

check(
  "No unresolved state references",
  all(state_references %in% system$state)
)

# Material coefficients must remain symbolic for later fitting.
check(
  "Material coefficients remain symbolic",
  any(grepl("^E_", symbols_used)) &&
    any(grepl("^D_", symbols_used))
)

# 4. Check that changing the fixed side reaches the generator.
top_system <- make_system("top")
top_fixed <- attr(top_system, "fixed_states")

check(
  "Fixed side can be changed",
  nrow(top_system) == 112L &&  # (35 - 7) * 4
    length(top_fixed) == 28L &&
    !("u_1_4_5" %in% top_system$state) &&
    "u_1_1_3" %in% top_system$state
)

cat("\nAll ODE-generator tests passed.\n")

# Compact inspection without printing the long force expressions.
print(head(system[, c(
  "index", "variable", "component", "boundary", "state"
)]))

