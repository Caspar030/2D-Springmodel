# R/equation_generator.R
#
# Assemble the full Kelvin–Voigt equation.
#
# Spatial fields and active-force definitions belong to fields.R.
# Fixed boundary states are eliminated by build_ode_system().
#
# Returns:
#   system           symbolic dynamic equations
#   parameter_values numerical defaults for coefficient symbols
#   initial_values   zero initial state
#   fixed_values     eliminated boundary states
#   expressions      expression trees
#
# No simulation or dMod2 conversion is performed here.


source("R/expression_tree.R")
source("R/operators.R")
source("R/fields.R")
source("R/discretization.R")
source("R/boundary_conditions.R")
source("R/ode_system.R")


generate_equations <- function(
    grid_bounds,
    spacing,
    fields = default_fields(2),
    parameters = list(),
    origin = c(0, 0),
    fixed_side = "left"
) {
  
  # ----------------------------------------------------------
  # 1. Validate geometry
  # ----------------------------------------------------------
  
  validate_grid_bounds(grid_bounds)
  
  if (!is.numeric(spacing) ||
      !is.null(dim(spacing)) ||
      length(spacing) != 2L ||
      any(!is.finite(spacing)) ||
      any(spacing <= 0)) {
    stop("spacing must contain two positive finite values.")
  }
  
  if (!is.numeric(origin) ||
      !is.null(dim(origin)) ||
      length(origin) != 2L ||
      any(!is.finite(origin))) {
    stop("origin must contain two finite coordinates.")
  }
  
  if (!is.list(parameters)) {
    stop("parameters must be a list.")
  }
  
  fixed_side <- match.arg(
    fixed_side,
    c("left", "right", "bottom", "top")
  )
  
  
  # ----------------------------------------------------------
  # 2. Bind fields to coefficient symbols
  # ----------------------------------------------------------
  
  binding <- make_field_binding(
    fields = fields,
    parameters = parameters,
    grid_bounds = grid_bounds,
    spacing = spacing,
    origin = origin
  )
  
  # Our current expression tree has no reciprocal or constant
  # operator. These two scalar placeholders implement the
  # algebra without changing operators.R or discretization.R.
  equation_symbols <- function(
    name,
    index,
    shape,
    component = integer(0),
    state_variables = c("u", "v"),
    fields = NULL,
    parameters = list()
  ) {
    if (name == "minus_one") {
      return("-1")
    }
    
    if (name == "inverse_rho") {
      rho_symbol <- binding$symbol(
        name = "rho",
        index = index,
        shape = integer(0),
        component = integer(0)
      )
      
      return(paste0("(1 / ", rho_symbol, ")"))
    }
    
    binding$symbol(
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
  # 3. Define symbolic fields
  # ----------------------------------------------------------
  
  u <- field("u", c(2))
  v <- field("v", c(2))
  
  E <- field("E", c(2, 2, 2, 2))
  D <- field("D", c(2, 2, 2, 2))
  
  inverse_rho <- field("inverse_rho")
  kappa <- field("kappa")
  minus_one <- field("minus_one")
  
  P <- field("P", c(2))
  
  
  # ----------------------------------------------------------
  # 4. Assemble the equation
  # ----------------------------------------------------------
  
  stress <- add(
    contract(E, grad(u)),
    contract(D, grad(v))
  )
  
  restoring_force <- multiply(
    minus_one,
    multiply(kappa, u)
  )
  
  total_force_density <- add(
    add(div(stress), restoring_force),
    P
  )
  
  expressions <- list(
    u = v,
    v = multiply(inverse_rho, total_force_density)
  )
  
  
  # ----------------------------------------------------------
  # 5. Generate componentwise ODEs
  # ----------------------------------------------------------
  
  system <- build_ode_system(
    expressions = expressions,
    grid_bounds = grid_bounds,
    spacing = spacing,
    fields = fields,
    parameters = parameters,
    field_symbol = equation_symbols,
    state_variables = c("u", "v"),
    fixed_side = fixed_side
  )
  
  
  # ----------------------------------------------------------
  # 6. Return model description
  # ----------------------------------------------------------
  
  list(
    system = system,
    
    parameter_values = binding$parameter_values(),
    
    initial_values = setNames(
      rep(0, nrow(system)),
      system$state
    ),
    
    fixed_values = attr(system, "fixed_states"),
    
    expressions = expressions,
    
    geometry = list(
      grid_bounds = grid_bounds,
      spacing = spacing,
      origin = origin
    ),
    
    fixed_side = fixed_side
  )
}