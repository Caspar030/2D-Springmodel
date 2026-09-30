library(dMod)
library(plotly)
# ToDo:
# The 1D-Solution comes for free as a side dish by implementing the 2D Version (correct) and setting i=1 

# Properly define E


# Layout: 1. Implement Nodes -> 2. implement all interactions of nodes (E, P, c, kappa) as distributions. 3. Define the interactions for each term 4. Model.




                                          # 1. - Defining States/Nodes ---


# ============================================================
# GRID SIZE: A \in \mathbb{R}^{Nx x Ny}

Nx <- 7
Ny <- 3


# ============================================================
# Choice of parameters # CAREFUL - THESE SHOULD be changed according to Nx and Ny
# ============================================================

parms <- c(
  #for E
  Ex      = 1,
  sigma_x = 2,
  Ey      = 10,
  sigma_y = 10,
  
  # for damping
  ceta      = 0.8,
  sigma_ceta_x = 20,
  sigma_ceta_y = 20,
  
  #for kappa
  kappa  = 0.05,  
  sigma_k_x = 20,
  sigma_k_y = 20,
  
  #for Q
  r = 1,
  sigma_Q = 10,
  
  #for 
  a0     = 1
  )




# I - Grid Spacing function a
      # The functions a_x and a_y return the spacing between mass points.
      # a_x(i,j) gives the spacing between mass points a_(i+1,j) and a_(i,j)




# a) Wave/Elasticity Term E: 
       # as a function of i,j. Building the operator with the divergence.



E <- function(x, y, Nx, Ny, sigma_x, sigma_y, Ex, Ey) {
  
  mux <- (Nx + 1) / 2
  muy <- (Ny + 1) / 2
  
  G_x <- dnorm(x, mux, sigma_x) / dnorm(mux, mux, sigma_x)
  G_y <- dnorm(y, muy, sigma_y) / dnorm(muy, muy, sigma_y)
  
  matrix(c(
    Ex * G_x, 0,
    0, Ey * G_y
  ), nrow = 2, byrow = TRUE)
}


# Test at the center
E(
  x = 4,
  y = 7,
  Nx = Nx,
  Ny = Ny,
  sigma_x = parms["sigma_x"],
  sigma_y = parms["sigma_y"],
  Ex = parms["Ex"],
  Ey = parms["Ey"]
)


                                          

#b) Damping Term c - can be modified inhomogenously.
      # c returns a normal distributed variable, centered, with mean z and standard deviations sigma

c <- function(x, y, z, Nx, Ny, sigma_k_x, sigma_k_y) {
  
  mux <- (Nx + 1) / 2
  muy <- (Ny + 1) / 2
  
  Gx <- dnorm(x, mean = mux, sd = sigma_x)
  Gy <- dnorm(y, mean = muy, sd = sigma_y)
  
  z * Gx * Gy /
    (dnorm(mux, mux, sigma_x) * dnorm(muy, muy, sigma_y))
}

#Testing c

c(4, 2, parms[["ceta"]], Nx, Ny, 2, 1)
# 2. - d) Retraction Force  - kappa






#c) Retraction Force kappa
      # k returns a random variable, centered, with mean z and standard deviations sigma_k_x and sigma_k_y
k <- function(x, y, z, Nx, Ny, sigma_k_x, sigma_k_y) {
  
  mux <- (Nx + 1) / 2
  muy <- (Ny + 1) / 2
  
  Gx <- dnorm(x, mean = mux, sd = sigma_x)
  Gy <- dnorm(y, mean = muy, sd = sigma_y)
  
  z * Gx * Gy /
    (dnorm(mux, mux, sigma_x) * dnorm(muy, muy, sigma_y))
}

#Testing k

k(4, 2, parms[["kappa"]], Nx, Ny, 2, 1)







#d) Active Contraction Force P
    #P is a set, oscillating function. May be replaced by more fitting data. Its Peaks are also distributed normally (This time 1D), centered,
    # the Gaussian has peak r as a scaling factor.
    #     Implementing Q(t).


Q_active <- "(0.5*(1 + tanh(8*sin(2*pi*time)))) * sqrt(2*pi) * parms[['sigma_Q']] * dnorm(parms[['r']], mean = (Nx + 1)/2, sd = parms[['sigma_Q']])"

times <- seq(0, 1, length.out = 1000)

Q <- sapply(times, function(time) {
  eval(parse(text = Q_active))
})

#Testing Q
range(Q)
min(Q)
max(Q)








# ============================================================
# Functions to construct states
# ============================================================

u_name <- function(i, j) {
  paste0("u_", i, "_", j)
}

v_name <- function(i, j) {
  paste0("v_", i, "_", j)
}


# ============================================================
# EQUATION CONTAINER

eqns <- c()

# ============================================================



#Beginning For Loop
for (i in 1:Nx) {

  for (j in 1:Ny) {
    
    # --------------------------------------------------------
    # Naming Variables with above functions.
    # --------------------------------------------------------
    
    uij <- u_name(i, j)
    vij <- v_name(i, j)
    
    
    # Left Boundary is fixed:
    
  if (j == 1) {
      uij <- "0"
      vij <- "0"
  }
   
    

    

    
    
    
    
# Constructing discrete Laplace. Ensure that edge cases with i==1, j==1 are not violated.



# ========================================================
# 2D Discrete Laplacian - "Fünfpunktlaplace"
#
#     d²u/dx² + d²u/dy²
#
#     ~ uR + uL + uU + uD - 4*uij
# ========================================================

    
    # ========================================================
    # Force Divergence: P
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
  by = 0.2
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
    frame = 50,
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


