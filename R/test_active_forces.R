source("R/model_script.R")

test_active_forces <- function() {
  
  check <- function(label, actual, expected) {
    if (!isTRUE(all.equal(
      as.numeric(actual),
      as.numeric(expected),
      tolerance = 1e-10
    ))) {
      stop(
        "FAILED: ", label,
        "\nActual: ", paste(actual, collapse = ", "),
        "\nExpected: ", paste(expected, collapse = ", ")
      )
    }
    
    cat("OK:", label, "\n")
  }
  
  fields <- default_fields(2)
  points <- grid_points(bounds)
  
  position <- function(index) {
    origin + (index - bounds$min) * spacing
  }
  
  volume <- function(index) {
    widths <- spacing
    on_boundary <- index == bounds$min | index == bounds$max
    widths[on_boundary] <- widths[on_boundary] / 2
    
    prod(widths) * thickness
  }
  
  # Convert peak force densities back to actual nodal forces.
  # Includes the fixed nodes.
  nodal_forces <- function(pars) {
    t(vapply(seq_len(nrow(points)), function(p) {
      index <- points[p, ]
      
      fields$P_amplitude(position(index), pars) * volume(index)
    }, numeric(2)))
  }
  
  
  # ----------------------------------------------------------
  # 1. Total active force must vanish
  # ----------------------------------------------------------
  
  forces <- nodal_forces(parameters)
  
  check(
    "Internal active forces balance over all nodes",
    colSums(forces),
    c(0, 0)
  )
  
  
  # ----------------------------------------------------------
  # 2. Weak sarcomere: equal and opposite force changes
  #
  # Matches model_script.R:
  # weak bond between (3,3) and (4,3).
  # ----------------------------------------------------------
  
  uniform_parameters <- parameters
  uniform_parameters$weak_fraction <- 1
  
  uniform_forces <- nodal_forces(uniform_parameters)
  difference <- forces - uniform_forces
  
  left_node <- which(points[, 1] == 3 & points[, 2] == 3)
  right_node <- which(points[, 1] == 4 & points[, 2] == 3)
  
  stopifnot(length(left_node) == 1L, length(right_node) == 1L)
  
  reduction <- parameters$Tmax * (1 - parameters$weak_fraction)
  
  expected_difference <- matrix(0, nrow(points), 2)
  
  # Weaker rightward pull on the left endpoint.
  expected_difference[left_node, ] <- c(-reduction, 0)
  
  # Weaker leftward pull on the right endpoint.
  expected_difference[right_node, ] <- c(reduction, 0)
  
  check(
    "Weak sarcomere changes only its two endpoint forces",
    difference,
    expected_difference
  )
  
  
  # ----------------------------------------------------------
  # 3. At u = v = 0:
  # du/dt = 0 and dv/dt = P(X,t) / rho(X)
  #
  # Evaluate the actual generated symbolic equations.
  # ----------------------------------------------------------
  
  system <- model$system
  
  expressions <- lapply(system$rhs, function(s) {
    parse(text = s)[[1]]
  })
  
  indices <- lapply(strsplit(system$index, ",", fixed = TRUE),
                    as.numeric)
  
  stopifnot(
    !anyDuplicated(names(model$parameter_values)),
    !any(names(model$parameter_values) %in% system$state)
  )
  
  env <- list2env(
    c(
      as.list(model$parameter_values),
      setNames(as.list(rep(0, nrow(system))), system$state)
    ),
    parent = baseenv()
  )
  
  pace <- pacing_parameters(parameters)
  period <- pace[["pace_period"]]
  duty <- pace[["pace_duty"]]
  start <- pace[["pace_start"]]
  
  test_times <- c(
    before_start = start - period / 10,
    pulse_peak = start + duty * period / 2,
    inactive = start + (duty + 1) * period / 2,
    next_peak = start + period + duty * period / 2
  )
  
  for (label in names(test_times)) {
    time <- unname(test_times[label])
    env$t <- time
    
    actual <- vapply(expressions, function(e) {
      eval(e, envir = env)
    }, numeric(1))
    
    expected <- vapply(seq_len(nrow(system)), function(r) {
      if (system$variable[r] == "u") {
        return(0)
      }
      
      X <- position(indices[[r]])
      
      acceleration <- fields$P(X, time, parameters) /
        fields$rho(X, parameters)
      
      acceleration[system$component[r]]
    }, numeric(1))
    
    check(
      paste("Generated ODE at zero state:", label),
      actual,
      expected
    )
  }
  
  cat("\nAll active-force tests passed.\n")
}

test_active_forces()