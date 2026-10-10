# ============================================================
# Test and build the model
# Run from the project root: ~/2D Springmodel
# ============================================================
source("R/expression_tree.R")
source("R/operators.R")
source("R/grid.R")
source("R/fields.R")
source("R/discretization.R")

# ------------------------------------------------------------
# Define fields
# ------------------------------------------------------------

u <- field("u")                  # displacement
v <- field("v")                  # velocity

E <- field("E", c(2, 2))         # elasticity tensor
D <- field("D", c(2, 2))         # viscosity tensor


# ------------------------------------------------------------
# Build the expression:
#
# div(E * grad(u) + D * grad(v))
# ------------------------------------------------------------

elastic_term <- matmul(E, grad(u))
viscous_term <- matmul(D, grad(v))

stress <- add(elastic_term, viscous_term)

rhs <- div(stress)


# ------------------------------------------------------------
# Inspect types
# ------------------------------------------------------------

show_type(grad(u))
show_type(elastic_term)
show_type(stress)
show_type(rhs)


# ------------------------------------------------------------
# Test type checking
# ------------------------------------------------------------

# This should produce an error:
# div(u)