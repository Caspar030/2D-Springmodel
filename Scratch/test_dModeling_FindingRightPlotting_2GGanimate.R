library(dMod)
library(ggplot2)
library(gganimate)


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



