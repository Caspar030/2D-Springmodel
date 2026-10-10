# R/dmod2_model.R
#
# Requires the object "model" from R/model_script.R.
#
# Steps:
#   1. Replace the symbolic pacing pulse by an external forcing.
#   2. Generate and compile the dMod2 model.
#   3. Assign parameters and initial values.
#   4. Simulate without sensitivities.
#
# The forcing specifies activation, NOT displacement.


library(dMod2)


# ------------------------------------------------------------
# 1. Check the generated model
# ------------------------------------------------------------

if (!exists("model") ||
    !is.list(model) ||
    !all(c(
      "system", "parameter_values", "initial_values"
    ) %in% names(model))) {
  stop("First run: source('R/model_script.R')")
}

# Preserve our model separately from the dMod2 ODE object.
kv_model <- model

rhs <- setNames(
  kv_model$system$rhs,
  kv_model$system$state
)

if (anyDuplicated(names(rhs))) {
  stop("ODE state names must be unique.")
}


# ------------------------------------------------------------
# 2. Replace the inline pulse by an external forcing
# ------------------------------------------------------------

pulse_text <- pacing_expression()

if (!any(grepl(pulse_text, rhs, fixed = TRUE))) {
  stop("The expected pacing expression was not found in the ODEs.")
}

rhs[] <- gsub(
  pattern = pulse_text,
  replacement = "pace_input",
  x = rhs,
  fixed = TRUE
)

used_symbols <- unique(unlist(lapply(rhs, function(s) {
  all.vars(parse(text = s))
})))

if (any(c(
  "pace_period", "pace_duty", "pace_start"
) %in% used_symbols)) {
  stop("Some inline pacing parameters remain after replacement.")
}


# ------------------------------------------------------------
# 3. Simulation times and forcing table
# ------------------------------------------------------------

times <- seq(0, 5, by = 0.01)

pace_names <- c(
  "pace_period", "pace_duty", "pace_start"
)

if (!all(pace_names %in% names(kv_model$parameter_values))) {
  stop("Pacing parameter values are missing.")
}

pace_parameters <- as.list(
  kv_model$parameter_values[pace_names]
)

# Validate the pacing configuration.
pace <- pacing_parameters(pace_parameters)

period <- unname(pace["pace_period"])
duty <- unname(pace["pace_duty"])
start <- unname(pace["pace_start"])

# At least approximately 100 intervals per active phase.
forcing_step <- min(
  min(diff(times)),
  period * duty / 100
)

cycles <- seq(
  floor((min(times) - start) / period) - 1,
  ceiling((max(times) - start) / period) + 1,
  by = 1
)

# Include pulse starts and ends explicitly.
edges <- c(
  start + cycles * period,
  start + (cycles + duty) * period
)

edges <- edges[
  edges >= min(times) & edges <= max(times)
]

forcing_times <- sort(unique(c(
  seq(min(times), max(times), by = forcing_step),
  times,
  edges
)))

forcing_data <- data.frame(
  name = "pace_input",
  time = forcing_times,
  value = pacing_pulse(
    forcing_times,
    parameters = pace_parameters
  )
)


# ------------------------------------------------------------
# 4. Generate and compile the ODE model
# ------------------------------------------------------------

outdir <- file.path(tempdir(), "kelvin_voigt_dmod2")

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

ode_model <- dMod2::odemodel(
  f = rhs,
  deriv = FALSE,
  forcings = "pace_input",
  modelname = "kelvin_voigt_ode",
  backend = "cppDE",
  compile = FALSE,
  outdir = outdir
)

x <- dMod2::Xs(
  ode_model,
  forcing_data,
  condition = "pacing",
  optionsOde = list(
    atol = 1e-9,
    rtol = 1e-7
  )
)

dMod2::compile(x)


# ------------------------------------------------------------
# 5. Match parameters and initial values
# ------------------------------------------------------------

available <- c(
  kv_model$parameter_values,
  kv_model$initial_values
)

if (anyDuplicated(names(available))) {
  stop("Parameter and initial-value names must be unique.")
}

required <- dMod2::getParameters(x)

missing_parameters <- setdiff(
  required,
  names(available)
)

if (length(missing_parameters)) {
  stop(
    "No numerical value available for: ",
    paste(missing_parameters, collapse = ", ")
  )
}

simulation_parameters <- available[required]

if (any(!is.finite(simulation_parameters))) {
  stop("All simulation parameters must be finite.")
}


# ------------------------------------------------------------
# 6. Simulate
# ------------------------------------------------------------

prediction <- x(
  times,
  simulation_parameters,
  deriv = FALSE
)

prediction_df <- as.data.frame(prediction)

cat("\nSimulation returned successfully.\n")
print(head(prediction_df))