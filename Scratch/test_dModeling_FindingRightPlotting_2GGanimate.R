library(dMod)
library(ggplot2)
library(gganimate)

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
#     Getting a 2D Plot of the Grid, but over time    ###
# ============================================================



grid_data <- data.frame()

sim_data <- sim[[1]]

for (t in seq_len(nrow(sim_data))) {
  
  for (j in 1:Ny) {
    for (i in 1:Nx) {
      
      x <- i
      y <- j
      
      if (i > 1) {
        x <- x + sim_data[t, u_name(i, j)]
      }
      
      grid_data <- rbind(
        grid_data,
        data.frame(
          time = sim_data[t, "time"],
          i = i,
          j = j,
          x = x,
          y = y
        )
      )
    }
  }
}
########      hihihi    ###
p <- ggplot(
  grid_data,
  aes(x = x, y = y)
) +
  geom_point(size = 3) +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = "t = {frame_time}",
    x = "x",
    y = "y"
  ) +
  transition_time(time) +
  ease_aes("linear")

p


####### or try this #########

grid_data$u <- ...

ggplot(grid_data, aes(x, y)) +
  geom_point(
    aes(color = u),
    size = 4
  ) +
  coord_fixed() +
  scale_color_gradient2() +
  transition_time(time) +
  labs(
    title = "t = {frame_time}",
    color = "Displacement"
  )


########## and this here should show the undeformed grid underneath #######
ggplot() +
  
  # Original grid
  geom_point(
    data = original_grid,
    aes(x, y),
    alpha = 0.2
  ) +
  
  # Deformed grid
  geom_point(
    data = grid_data,
    aes(x, y),
    color = "black"
  ) +
  
  coord_fixed() +
  
  transition_time(time)



