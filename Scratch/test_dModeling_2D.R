library(dMod2)
library(plotly)
# ToDo:
# The 1D-Solution comes for free as a side dish by implementing the 2D Version (correct) and setting i=1 

# Properly define E


# Layout: 1. Implement Nodes -> 2. implement all interactions of nodes (E, P, c, kappa) as distributions. 3. Define the interactions for each term 4. Model.





# ============================================================


                                          # 1. - Defining States/Nodes ---

#        i = 1      i = 2      i = 3      ...      i = Nx
#
# j = 1  a_1_1      a_2_1      a_3_1      ...      a_Nx_1
# j = 2  a_1_2      a_2_2      a_3_2      ...      a_Nx_2
# j = 3  a_1_3      a_2_3      a_3_3      ...      a_Nx_3
#  ...     ...        ...        ...        ...             ...
# j = Ny  a_1_Ny     a_2_Ny     a_3_Ny     ...      a_Nx_Ny




Nx <- 7               # for j
Ny <- 3               # for i





# ============================================================
# Choice of parameters # CAREFUL - THESE SHOULD be changed according to Nx and Ny
# ============================================================

parms <- c(
  #spacing function a
  sigma_upper = 10,
  sigma_right = 20,
  
  
  #for E
  Ex      = 1,
  sigma_E_x = 2,
  Ey      = 10,
  sigma_E_y = 10,
  
  #for D
  Dx      = 1,
  sigma_D_x = 2,
  Dy      = 10,
  sigma_D_y = 10,
  
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


#                   Implementing all sorts of functions necessary for constructing the PDE     ###




#Mass density function      #here also Gaussian implement
      mass_density <- function(x, y) {
        dnorm(x, mean = (Nx+1)/2, sd = 20) *
          dnorm(y, mean = (Ny + 1)/2, sd = 10)/ (dnorm((Nx+1)/2 , mean = (Nx+1)/2, sd = 20) *
          dnorm((Ny + 1)/2, mean = (Ny + 1)/2, sd = 10))
      }





# I - Grid Spacing function a
      # The functions a_x and a_y return the spacing between mass points.
      # a_x(i,j,o) gives the spacing between mass points a_(i+1,j) 

a <- function(x, y, orientation) {
  
  if (orientation == 0) {
    # right neighbor
    exp(-((x)^2 + (y - 1)^2) / (2 * parms["sigma_right"]^2))
    
  } else if (orientation == 1) {
    # upper neighbor
    exp(-((x - 1)^2 + (y)^2) / (2 * parms["sigma_upper"]^2))
  }
}







# a) Wave/Elasticity Term E: 
       # as a function of i,j. Building the operator with the divergence.



E <- function(x, y, Nx, Ny, sigma_E_x, sigma_E_y, Ex, Ey) {
  
  mux <- (Nx + 1) / 2
  muy <- (Ny + 1) / 2
  
  G_x <- dnorm(x, mux, sigma_E_x) / dnorm(mux, mux, sigma_E_x)
  G_y <- dnorm(y, muy, sigma_E_y) / dnorm(muy, muy, sigma_E_y)
  
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
  sigma_E_x = parms["sigma_E_x"],
  sigma_E_y = parms["sigma_E_y"],
  Ex = parms["Ex"],
  Ey = parms["Ey"]
)


                                          

# b) Klevin Voigt - Damping Term E: 
# as a function of i,j. Building the operator with the divergence.



D <- function(x, y, Nx, Ny, sigma_D_x, sigma_D_y, Dx, Dy) {
  
  mux <- (Nx + 1) / 2
  muy <- (Ny + 1) / 2
  
  G_x <- dnorm(x, mux, sigma_D_x) / dnorm(mux, mux, sigma_D_x)
  G_y <- dnorm(y, muy, sigma_D_y) / dnorm(muy, muy, sigma_D_y)
  
  matrix(c(
    Dx * G_x, 0,
    0, Dy * G_y
  ), nrow = 2, byrow = TRUE)
}


# Test at the center
D(
  x = 4,
  y = 7,
  Nx = Nx,
  Ny = Ny,
  sigma_D_x = parms["sigma_D_x"],
  sigma_D_y = parms["sigma_D_y"],
  Dx = parms["Dx"],
  Dy = parms["Dy"]
)





#b) Damping Term ceta - can be modified inhomogenously.
      # c returns a normal distributed variable, centered, with mean z and standard deviations sigma

ceta <- function(x, y, z, Nx, Ny, sigma_k_x, sigma_k_y) {
  
  mux <- (Nx + 1) / 2
  muy <- (Ny + 1) / 2
  
  Gx <- dnorm(x, mean = mux, sd = sigma_k_x)
  Gy <- dnorm(y, mean = muy, sd = sigma_k_y)
  
  z * Gx * Gy /
    (dnorm(mux, mux, sigma_k_x) * dnorm(muy, muy, sigma_k_y))
}

#Testing c

ceta(4, 2, parms[["ceta"]], Nx, Ny, 2, 1)
# 2. - d) Retraction Force  - kappa






#c) Retraction Force kappa
      # k returns a random variable, centered, with mean z and standard deviations sigma_k_x and sigma_k_y
k <- function(x, y, z, Nx, Ny, sigma_k_x, sigma_k_y) {
  
  mux <- (Nx + 1) / 2
  muy <- (Ny + 1) / 2
  
  Gx <- dnorm(x, mean = mux, sd = sigma_k_x)
  Gy <- dnorm(y, mean = muy, sd = sigma_k_y)
  
  z * Gx * Gy /
    (dnorm(mux, mux, sigma_k_x) * dnorm(muy, muy, sigma_k_y))
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











####################    Functions to construct states    ####################

# These functions already inherit the boundary condition: Left == 0.

u_name <- function(i, j) {
  if (j == 1) {
    return("0")
  } else {
    return(paste0("u_", i, "_", j))
  }
}

v_name <- function(i, j) {
  if (j == 1) {
    return("0")
  } else {
    return(paste0("v_", i, "_", j))
  }
}








####################    Gradient functions    #################### (might it be prettier to introduce one gradient function that takes other functions as input?)

grad_u <- function(i, j) {
  c(
    paste0(
      "(", u_name(i + 1, j), " - ", u_name(i, j), ") / ",
      a(i, j, 0)   # upper neighbour
    ),
    
    paste0(
      "(", u_name(i, j + 1), " - ", u_name(i, j), ") / ",
      a(i, j, 1)   # right neighbour
    )
  )
}

grad_v <- function(i, j) {
  c(
    paste0(
      "(", v_name(i + 1, j), " - ", v_name(i, j), ") / ",
      a(i, j, 0)   # upper neighbour
    ),
    
    paste0(
      "(", v_name(i, j + 1), " - ", v_name(i, j), ") / ",
      a(i, j, 1)   # right neighbour
    )
  )
}



# v and w specify the entry of E
grad_E <- function(i, j, v, w) {
  c(
    paste0(
      "(", E(i + 1, j, Nx, Ny, sigma_E_x, sigma_E_y, Ex, Ey)[v, w],
      " - ",
      E(i, j, Nx, Ny, sigma_E_x, sigma_E_y, Ex, Ey)[v, w],
      ") / ",
      a(i, j, 0)
    ),  #Upper neighbour
      
    
    paste0(
      "(", E(i, j + 1, Nx, Ny, sigma_E_x, sigma_E_y, Ex, Ey)[v, w],
      " - ",
      E(i, j, Nx, Ny, sigma_E_x, sigma_E_y, Ex, Ey)[v, w],
      ") / ",
      a(i, j, 1)
    )  # Right neighbour
  )
}

# v and w specify the entry of D
grad_D <- function(i, j, v, w) {
  c(
    paste0(
      "(", D(i + 1, j, Nx, Ny, sigma_D_x, sigma_D_y, Dx, Dy)[v, w],
      " - ",
      D(i, j, Nx, Ny, sigma_D_x, sigma_D_y, Dx, Dy)[v, w],
      ") / ",
      a(i, j, 0)
    ),  #Upper neighbour
    
    
    paste0(
      "(", D(i, j + 1, Nx, Ny, sigma_D_x, sigma_D_y, Dx, Dy)[v, w],
      " - ",
      D(i, j, Nx, Ny, sigma_D_x, sigma_D_y, Dx, Dy)[v, w],
      ") / ",
      a(i, j, 1)
    )  # Right neighbour
  )
}

######### Second Derivative Function #### works only in the interior!!!!###

partial_x2_u <- function(i, j) {
  paste0(
    "(",
    u_name(i, j + 1 ), " + ",
    u_name(i, j - 1), " - 2 * ",
    u_name(i, j),
    ") / ",
    a(i, j, 0), "^2"
  )
}

partial_y2_u <- function(i, j) {
  paste0(
    "(",
    u_name(i + 1, j), " + ",
    u_name(i - 1, j), " - 2 * ",
    u_name(i, j),
    ") / ",
    a(i, j, 1), "^2"
  )
}

# Cross-Terms: partial_x partial_y


####################    E * grad(u)    ####################

E_grad_u <- function(i, j) {
  
  grad <- grad_u(i, j)
  
  Eij <- E(
    i, j,
    Nx, Ny,
    parms["sigma_E_x"],
    parms["sigma_E_y"],
    parms["Ex"],
    parms["Ey"]
  )
  
  c(
    paste0(
      Eij[1,1], " * (", grad[1], ") + ",
      Eij[1,2], " * (", grad[2], ")"
    ),
    
    paste0(
      Eij[2,1], " * (", grad[1], ") + ",
      Eij[2,2], " * (", grad[2], ")"
    )
  )
}

####################    D * grad(v)    ####################

D_grad_v <- function(i, j) {
  
  grad <- grad_v(i, j)
  
  Dij <- D(
    i, j,
    Nx, Ny,
    parms["sigma_D_x"],
    parms["sigma_D_y"],
    parms["Dx"],
    parms["Dy"]
  )
  
  c(
    paste0(
      Dij[1,1], " * (", grad[1], ") + ",
      Dij[1,2], " * (", grad[2], ")"
    ),
    
    paste0(
      Dij[2,1], " * (", grad[1], ") + ",
      Dij[2,2], " * (", grad[2], ")"
    )
  )
}

########## Divergence Operator          ###############




####################    Equation container    ####################

eqns <- c()














########                    FOR LOOP STARTS     #################

####################    Building equations    ####################

for (i in 1:Nx) {
  
  for (j in 1:Ny) {
    
    # --------------------------------------------------------
    # Naming variables
    # --------------------------------------------------------
    
    uij <- u_name(i, j)
    vij <- v_name(i, j)
    
    
    # --------------------------------------------------------
    # Left boundary is fixed
    # --------------------------------------------------------
    
    if (j == 1) {
      uij <- "0"
      vij <- "0"
    }
    
    
    # --------------------------------------------------------
    # Differential operators
    # Only calculate if upper AND right neighbour exist
    # --------------------------------------------------------
    
    if (i !=1 && i != Nx && j != 1 && j != Ny) {
      
      grad_uij <- grad_u(i, j)
      grad_vij <- grad_v(i, j)
      
      E_grad_uij <- E_grad_u(i, j)
      D_grad_vij <- D_grad_v(i, j)
      
      # Testing
      print(E_grad_uij)
      print(D_grad_uij)
    }
    
    
  } #######       For loop ending brackets.
}
    
    
    

    
    
    
  # And now we only need the divergence of these two :)



    
    
    
    
    
    




    
    
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
# Startbedingungen
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


