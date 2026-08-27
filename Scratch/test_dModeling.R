library(dMod)

# ============================================================
# 2D spring / elastic grid
# ============================================================

Nx <- 10
Ny <- 10

# Helper function for state names
u_name <- function(i, j) paste0("u_", i, "_", j)
v_name <- function(i, j) paste0("v_", i, "_", j)

eqns <- c()

# ============================================================
# Build ODE system
# Only interior points are dynamic
# Boundary: u = 0
# ============================================================

for (i in 2:(Nx - 1)) {
  
  for (j in 2:(Ny - 1)) {
    
    uij <- u_name(i, j)
    vij <- v_name(i, j)
    
    # Neighbours
    uR <- u_name(i + 1, j)
    uL <- u_name(i - 1, j)
    uU <- u_name(i, j + 1)
    uD <- u_name(i, j - 1)
    
    # Replace boundary states by 0
    if (i + 1 == Nx) uR <- "0"
    if (i - 1 == 1)  uL <- "0"
    
    if (j + 1 == Ny) uU <- "0"
    if (j - 1 == 1)  uD <- "0"
    
    
    # --------------------------------------------------------
    # 1. du/dt = v
    # --------------------------------------------------------
    
    eqns <- c(
      eqns,
      paste0("dot(", uij, ") = ", vij)
    )
    
    
    # --------------------------------------------------------
    # 2. 2D discrete Laplacian
    #
    # u_R + u_L + u_U + u_D - 4*u_ij
    # --------------------------------------------------------
    
    laplace <- paste0(
      "(", uR,
      " + ", uL,
      " + ", uU,
      " + ", uD,
      " - 4*", uij, ")"
    )
    
    
    # --------------------------------------------------------
    # 3. Active/passive forces
    #
    # lambda_x = 1 + (u_R - u_ij)/a0
    # lambda_y = 1 + (u_U - u_ij)/a0
    #
    # P_x = eta * lambda_x
    # P_y = eta * lambda_y
    #
    # For simplicity:
    # dP/dx + dP/dy
    # --------------------------------------------------------
    
    force_x <- paste0(
      "eta/a0 * (",
      "(", uR, " - ", uij, ")",
      " - (", uij, " - ", uL, ")",
      ")"
    )
    
    force_y <- paste0(
      "eta/a0 * (",
      "(", uU, " - ", uij, ")",
      " - (", uij, " - ", uD, ")",
      ")"
    )
    
    
    # --------------------------------------------------------
    # 4. dv/dt
    # --------------------------------------------------------
    
    acceleration <- paste0(
      "(E/(delta0*a0^2))*", laplace,
      " + (1/delta0)*(", force_x,
      " + ", force_y, ")",
      " - (c/delta0)*", vij,
      " - (kappa/delta0)*", uij
    )
    
    eqns <- c(
      eqns,
      paste0(
        "dot(", vij, ") = ",
        acceleration
      )
    )
  }
}


# ============================================================
# Create dMod model
# ============================================================

model <- odemodel(
  eqns,
  modelname = "spring2D"
)


# ============================================================
# Parameters
# ============================================================

parms <- c(
  E      = 1,
  delta0 = 1,
  a0     = 1,
  eta    = 0.2,
  c      = 0.1,
  kappa  = 0.05
)


# ============================================================
# Initial conditions
# ============================================================

x0 <- c()

# Gaussian displacement in the middle of the grid
for (i in 2:(Nx - 1)) {
  
  for (j in 2:(Ny - 1)) {
    
    uij <- u_name(i, j)
    vij <- v_name(i, j)
    
    # Initial Gaussian bump
    u0 <- exp(
      -(
        (i - Nx/2)^2 +
          (j - Ny/2)^2
      ) / 4
    )
    
    x0[uij] <- u0
    x0[vij] <- 0
  }
}


# ============================================================
# Simulation
# ============================================================

times <- seq(
  0,
  20,
  by = 0.1
)

sim <- simulate(
  model,
  times = times,
  parms = parms,
  x0 = x0
)