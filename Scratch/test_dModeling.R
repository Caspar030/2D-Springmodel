library(dMod)

# ============================================================
# GRID SIZE
# ============================================================

Nx <- 7
Ny <- 3


# ============================================================
# STATE NAMES
# ============================================================

u_name <- function(i, j) {
  paste0("u_", i, "_", j)
}

v_name <- function(i, j) {
  paste0("v_", i, "_", j)
}


# ============================================================
# EQUATION CONTAINER
# ============================================================

eqns <- c()


# ============================================================
# BUILD THE ODE SYSTEM
#
# Boundary conditions:
#
# LEFT:
#     u = 0  --> fixed
#
# RIGHT:
#     free
#
# TOP:
#     free
#
# BOTTOM:
#     free
# ============================================================

for (i in 2:Nx) {
  
  for (j in 1:Ny) {
    
    # --------------------------------------------------------
    # STATE NAMES
    # --------------------------------------------------------
    
    uij <- u_name(i, j)
    vij <- v_name(i, j)
    
    
    # --------------------------------------------------------
    # LEFT NEIGHBOUR
    #
    # If i = 2, the neighbour i = 1 is the fixed boundary:
    #
    #     u_(1,j) = 0
    # --------------------------------------------------------
    
    if (i == 2) {
      uL <- "0"
    } else {
      uL <- u_name(i - 1, j)
    }
    
    
    # --------------------------------------------------------
    # RIGHT NEIGHBOUR
    #
    # Free boundary:
    #
    #     u_(Nx+1,j) = u_(Nx,j)
    #
    # corresponding to du/dx = 0.
    # --------------------------------------------------------
    
    if (i == Nx) {
      uR <- uij
    } else {
      uR <- u_name(i + 1, j)
    }
    
    
    # --------------------------------------------------------
    # LOWER NEIGHBOUR
    #
    # Free boundary:
    #
    #     u_(i,0) = u_(i,1)
    # --------------------------------------------------------
    
    if (j == 1) {
      uD <- uij
    } else {
      uD <- u_name(i, j - 1)
    }
    
    
    # --------------------------------------------------------
    # UPPER NEIGHBOUR
    #
    # Free boundary:
    #
    #     u_(i,Ny+1) = u_(i,Ny)
    # --------------------------------------------------------
    
    if (j == Ny) {
      uU <- uij
    } else {
      uU <- u_name(i, j + 1)
    }
    
    
    # ========================================================
    # 2D DISCRETE LAPLACIAN
    #
    #     d²u/dx² + d²u/dy²
    #
    #     = uR + uL + uU + uD - 4*uij
    # ========================================================
    
    laplace <- paste0(
      "(",
      uR, " + ",
      uL, " + ",
      uU, " + ",
      uD, " - 4*", uij,
      ")"
    )
    
    
    # ========================================================
    # X-DIRECTION FORCE
    # ========================================================
    
    # --------------------------------------------------------
    # Outgoing force on right face
    # --------------------------------------------------------
    
    if (i == Nx) {
      
      # Free right boundary:
      #
      #     P_right = 0
      
      Px_right <- "0"
      
    } else {
      
      # Spring between (i,j) and (i+1,j):
      #
      # lambda = 1 + (uR - uij)/a0
      #
      # P = eta * lambda
      
      Px_right <- paste0(
        "eta*(1 + (",
        uR,
        " - ",
        uij,
        ")/a0)"
      )
    }
    
    
    # --------------------------------------------------------
    # Incoming force from left face
    # --------------------------------------------------------
    
    if (i == 2) {
      
      # Spring between fixed boundary and first dynamic point:
      #
      # left displacement = 0
      
      Px_left <- paste0(
        "eta*(1 + (",
        uij,
        " - 0)/a0)"
      )
      
    } else {
      
      # Spring between (i-1,j) and (i,j)
      
      Px_left <- paste0(
        "eta*(1 + (",
        uij,
        " - ",
        uL,
        ")/a0)"
      )
    }
    
    
    # ========================================================
    # Y-DIRECTION FORCE
    # ========================================================
    
    # --------------------------------------------------------
    # Upper face
    # --------------------------------------------------------
    
    if (j == Ny) {
      
      # Free upper boundary:
      #
      #     P_upper = 0
      
      Py_up <- "0"
      
    } else {
      
      # Spring between (i,j) and (i,j+1)
      
      Py_up <- paste0(
        "eta*(1 + (",
        uU,
        " - ",
        uij,
        ")/a0)"
      )
    }
    
    
    # --------------------------------------------------------
    # Lower face
    # --------------------------------------------------------
    
    if (j == 1) {
      
      # Free lower boundary:
      #
      #     P_lower = 0
      
      Py_down <- "0"
      
    } else {
      
      # Spring between (i,j-1) and (i,j)
      
      Py_down <- paste0(
        "eta*(1 + (",
        uij,
        " - ",
        uD,
        ")/a0)"
      )
    }
    
    
    # ========================================================
    # FORCE DIVERGENCE
    # ========================================================
    
    force_div <- paste0(
      "((",
      Px_right,
      ") - (",
      Px_left,
      "))",
      " + ",
      "((",
      Py_up,
      ") - (",
      Py_down,
      "))"
    )
    
    
    # ========================================================
    # ACCELERATION
    #
    #     dv/dt = acceleration
    # ========================================================
    
    acceleration <- paste0(
      "(E/(delta0*a0^2))*",
      laplace,
      " + (1/(delta0*a0))*(",
      force_div,
      ")",
      " - (c/delta0)*",
      vij,
      " - (kappa/delta0)*",
      uij
    )
    
    
    # ========================================================
    # ADD EQUATIONS TO eqns
    #
    # dMod uses the names of the equations to identify the
    # corresponding state variables.
    #
    #     du/dt = v
    #     dv/dt = acceleration
    #
    # Therefore:
    #
    #     u_2_1 = "v_2_1"
    #     v_2_1 = "acceleration expression"
    # ========================================================
    
    eqns <- c(
      eqns,
      setNames(
        vij,
        uij
      ),
      setNames(
        acceleration,
        vij
      )
    )
  }
}


# ============================================================
# OPTIONAL: INSPECT GENERATED EQUATIONS
# ============================================================

print(eqns)


# ============================================================
# CREATE THE DMOD MODEL
# ============================================================

model <- odemodel(
  eqns,
  modelname = "spring2D_mixedBC"
)


# ============================================================
# PARAMETERS
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
# INITIAL CONDITIONS
# ============================================================

x0 <- c()

for (i in 2:Nx) {
  
  for (j in 1:Ny) {
    
    uij <- u_name(i, j)
    vij <- v_name(i, j)
    
    # --------------------------------------------------------
    # Initial Gaussian displacement
    # --------------------------------------------------------
    
    x0[uij] <- exp(
      -(
        (i - 4)^2 +
          (j - Ny/2)^2
      ) / 4
    )
    
    # --------------------------------------------------------
    # Initially at rest
    # --------------------------------------------------------
    
    x0[vij] <- 0
  }
}


# ============================================================
# SIMULATION
# ============================================================



times <- seq(
  0,
  50,
  by = 0.1
)


x <- Xs(model)

pars <- c(x0, parms)

sim <- x(times, pars)



# ============================================================
# Plot-Grid
# ============================================================
# 
# 
# 
# 
# plot_grid <- function(sim, time_index) {
#   
#   sim_data <- sim[[1]]
#   
#   U <- matrix(
#     0,
#     nrow = Nx,
#     ncol = Ny
#   )
#   
#   for (i in 2:Nx) {
#     for (j in 1:Ny) {
#       
#       state <- u_name(i, j)
#       
#       U[i, j] <- sim_data[time_index, state]
#     }
#   }
#   
#   # Color scale
#   cols <- heat.colors(100)
#   
#   image(
#     x = 1:Nx,
#     y = 1:Ny,
#     z = U,
#     col = cols,
#     xlab = "x",
#     ylab = "y",
#     main = paste0(
#       "Displacement, t = ",
#       round(sim_data[time_index, "time"], 2)
#     )
#   )
#   
#   # Fixed boundary
#   points(
#     rep(1, Ny),
#     1:Ny,
#     pch = 19
#   )
#   
#   # Color legend
#   legend(
#     "topright",
#     legend = round(
#       seq(min(U), max(U), length.out = 5),
#       3
#     ),
#     fill = cols[
#       round(
#         seq(1, length(cols), length.out = 5)
#       )
#     ],
#     title = "u",
#     bty = "n"
#   )
# }


###############       different plot    ############
#
#     Getting a 2D Plot of the Grid, but over time    ###
#



plot_grid <- function(sim, time_index) {
  
  sim_data <- sim[[1]]
  
  # Original grid
  X <- matrix(rep(1:Nx, Ny), nrow = Ny, byrow = TRUE)
  Y <- matrix(rep(1:Ny, each = Nx), nrow = Ny, byrow = TRUE)
  
  # Add displacement
  X_def <- X
  Y_def <- Y
  
  for (j in 1:Ny) {
    for (i in 2:Nx) {
      
      X_def[j, i] <- X[j, i] +
        sim_data[time_index, u_name(i, j)]
      
      # If you also have y-displacement:
      # Y_def[j, i] <- Y[j, i] +
      #   sim_data[time_index, v_name(i, j)]
    }
  }
  
  plot(
    X_def,
    Y_def,
    type = "n",
    asp = 1,
    xlab = "x",
    ylab = "y",
    main = paste(
      "t =",
      sim_data[time_index, "time"]
    )
  )
  
  # Draw horizontal connections
  for (j in 1:Ny) {
    lines(X_def[j, ], Y_def[j, ])
  }
  
  # Draw vertical connections
  for (i in 1:Nx) {
    lines(X_def[, i], Y_def[, i])
  }
  
  # Draw mass points
  points(
    X_def,
    Y_def,
    pch = 19
  )
}

# ============================================================
# EXAMPLE PLOTS
# ============================================================

plot_grid(sim, 1)
plot_grid(sim, 10)
plot_grid(sim, 20)
plot_grid(sim, 40)
plot_grid(sim, 50)
plot_grid(sim, 60)
plot_grid(sim, 70)
plot_grid(sim, 80)
plot_grid(sim, 90)
plot_grid(sim, 100)
plot_grid(sim, 400)




