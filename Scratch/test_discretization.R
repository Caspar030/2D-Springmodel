# Expression
u <- field("u")
E <- field("E", c(2, 2))
rhs <- div(matmul(E, grad(u)))

grid <- list(
  min = c(1, 1),
  max = c(7, 5)
)

# Discretize at upper-right corner
result <- discretize(
  expr = rhs,
  index = c(7, 5),
  spacing = c(1, 1),
  grid_bounds = grid
)

# Evaluation environment
env <- new.env(parent = baseenv())

# u(x,y) = x^2 + y^2
for (i in 1:7) {
  for (j in 1:5) {
    assign(
      paste0("u_", i, "_", j),
      i^2 + j^2,
      envir = env
    )
  }
}

# E = identity matrix
for (i in 1:7) {
  for (j in 1:5) {
    assign(paste0("E_1_1_", i, "_", j), 1, envir = env)
    assign(paste0("E_1_2_", i, "_", j), 0, envir = env)
    assign(paste0("E_2_1_", i, "_", j), 0, envir = env)
    assign(paste0("E_2_2_", i, "_", j), 1, envir = env)
  }
}

# Evaluate
value <- eval(parse(text = result), envir = env)

cat("Numerical result:", value, "\n")
cat("Exact result:    ", 4, "\n")
cat("Error:           ", abs(value - 4), "\n")
