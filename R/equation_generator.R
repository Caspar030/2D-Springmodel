source("R/expression_tree.R")
source("R/operators.R")
source("R/fields.R")
source("R/discretization.R")
source("R/boundary_conditions.R")
source("R/ode_system.R")

# Fields
u <- field("u")
v <- field("v")
E <- field("E", c(2, 2))

# Example spatial operators
rhs_u <- div(matmul(E, grad(u)))
rhs_v <- div(matmul(E, grad(v)))

# Grid
grid <- list(
  min = c(1, 1),
  max = c(7, 5)
)

# Generate system
system <- build_ode_system(
  expressions = list(
    u = rhs_u,
    v = rhs_v
  ),
  grid_bounds = grid,
  spacing = c(1, 1),
  fields = default_fields(2)
)

# Inspect results
head(system)
table(system$boundary)
