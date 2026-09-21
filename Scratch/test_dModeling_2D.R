library(dMod)
library(plotly)


# ============================================================
# GRID SIZE:
# u_1,1  ... u_1,Nx
# ...        ...
# u_Nx_1 ... u_Nx_Ny
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
# Buildung the coupled ODE System
#
# Boundary conditions:
#
# Left:
#     u = 0  --> fixed
#
# Right:
#     free
#
# Top:
#     free
#
# Bottom:
#     free
# ============================================================

for (i in 2:Nx) {
  
  for (j in 1:Ny) {
    
    # --------------------------------------------------------
    # Naming Variables:
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
    # This corresponds to du/dx = 0.
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
    # 2D Discrete Laplacian - "Fünfpunktlaplace"
    #
    #     d²u/dx² + d²u/dy²
    #
    #     ~ uR + uL + uU + uD - 4*uij
    # ========================================================
    
    laplace <- paste0(
      "(",
      uR, " + ",
      uL, " + ",
      uU, " + ",
      uD, " - 4*", uij,
      ")"
    )
    
    
    
    
    
    # Here, the P-Term is implemented. In this Version, Q(t), that is, the interesting function is not implemented yet.
    # ========================================================
    # X-Direction Force
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
    # Y-Direction Force
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
    # Force Divergence
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
    # Acceleration
    #
    #     The whole equation is as follows
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
    # Also:
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
# Printing generated equations
# ============================================================

print(eqns)


# ============================================================
# CCreating the dMod Model
# ============================================================

model <- odemodel(
  eqns,
  modelname = "spring2D_mixedBC"
)


# ============================================================
# Choice of parameters
# ============================================================

parms <- c(
  E      = 1,
  delta0 = 1,
  a0     = 1,
  eta    = 0.2,
  c      = 0.8,
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
# Simulation
# ============================================================



times <- seq(
  0,
  200,
  by = 0.1
)


x <- Xs(model)

pars <- c(x0, parms)

sim <- x(times, pars)



# ============================================================
# Plot-Grid
#     Getting a 2D Plot of the Grid, but over time    ###
# ============================================================



plot_grid_interactive <- function(sim) {
  
  sim_data <- sim[[1]]
  
  # ------------------------------------------------------------
  # Build data for ALL time points
  # ------------------------------------------------------------
  
  all_points <- data.frame()
  all_lines  <- data.frame()
  
  for (t in seq_len(nrow(sim_data))) {
    
    X_def <- matrix(1:Nx, nrow = Ny, ncol = Nx, byrow = TRUE)
    Y_def <- matrix(1:Ny, nrow = Ny, ncol = Nx, byrow = FALSE)
    
    # Apply x-displacement
    for (j in 1:Ny) {
      for (i in 2:Nx) {
        
        X_def[j, i] <- X_def[j, i] +
          sim_data[t, u_name(i, j)]
      }
    }
    
    # ----------------------------------------------------------
    # Points
    # ----------------------------------------------------------
    
    points_df <- data.frame(
      x = as.vector(X_def),
      y = as.vector(Y_def),
      time = t
    )
    
    all_points <- rbind(all_points, points_df)
    
    
    # ----------------------------------------------------------
    # Horizontal lines
    # ----------------------------------------------------------
    
    for (j in 1:Ny) {
      
      lines_df <- data.frame(
        x = X_def[j, ],
        y = Y_def[j, ],
        group = paste0("h", j),
        time = t
      )
      
      all_lines <- rbind(all_lines, lines_df)
    }
    
    
    # ----------------------------------------------------------
    # Vertical lines
    # ----------------------------------------------------------
    
    for (i in 1:Nx) {
      
      lines_df <- data.frame(
        x = X_def[, i],
        y = Y_def[, i],
        group = paste0("v", i),
        time = t
      )
      
      all_lines <- rbind(all_lines, lines_df)
    }
  }
  
  
  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------
  
  p <- plot_ly()
  
  
  # ------------------------------------------------------------
  # Add grid lines
  # ------------------------------------------------------------
  
  for (g in unique(all_lines$group)) {
    
    tmp <- all_lines[
      all_lines$group == g,
    ]
    
    p <- add_trace(
      p,
      data = tmp,
      x = ~x,
      y = ~y,
      type = "scatter",
      mode = "lines",
      frame = ~time,
      line = list(width = 1),
      showlegend = FALSE
    )
  }
  
  
  # ------------------------------------------------------------
  # Add mass points
  # ------------------------------------------------------------
  
  p <- add_trace(
    p,
    data = all_points,
    x = ~x,
    y = ~y,
    type = "scatter",
    mode = "markers",
    frame = ~time,
    marker = list(size = 8),
    showlegend = FALSE
  )
  
  
  # ------------------------------------------------------------
  # Animation controls
  # ------------------------------------------------------------
  
  p <- animation_opts(
    p,
    frame = 2,
    transition = 0,
    redraw = TRUE
  )
  
  p <- animation_slider(
    p,
    currentvalue = list(
      prefix = "Time: "
    )
  )
  
  p <- animation_button(
    p,
    x = 1,
    xanchor = "right",
    y = 0,
    yanchor = "top"
  )
  
  
  # ------------------------------------------------------------
  # Axis settings
  # ------------------------------------------------------------
  
  p <- layout(
    p,
    xaxis = list(
      title = "x",
      range = c(0, Nx + 1)
    ),
    yaxis = list(
      title = "y",
      range = c(0, Ny + 1),
      scaleanchor = "x",
      scaleratio = 1
    )
  )
  
  p
}


# Plotting the grid (may take some time)
plot_grid_interactive(sim)


