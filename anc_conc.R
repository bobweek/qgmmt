microbial_lineage_grid_plot <- function(n_hosts = 10, n_gens = 20,
                                        lineage_path,
                                        host_lineage = NULL,
                                        environment_col = n_hosts + 1,
                                        title = "Microbial Lineage Ancestry") {
  
  plot(NA, xlim = c(0.5, environment_col + 0.5), ylim = c(0.5, n_gens + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = title, asp = 1)
  
  axis(2, at = 1:n_gens, labels = paste0("", 1:n_gens), las = 1, cex.axis = 0.6)
  
  for (gen in 1:n_gens) {
    points(1:n_hosts, rep(gen, n_hosts), pch = 21, bg = "gray90", col = "black")
    rect(environment_col - 0.4, gen - 0.4, environment_col + 0.4, gen + 0.4,
         col = rgb(0.9, 0.95, 1, 0.3), border = "lightblue")
  }
  text(environment_col, n_gens + 0.75, "Environment", cex = 0.8)
  
  if (!is.null(host_lineage)) {
    host_x <- sapply(host_lineage, function(pos) pos[2])
    host_y <- sapply(host_lineage, function(pos) pos[1])
    lines(host_x, host_y, col = "darkgreen", lwd = 2)
    points(host_x, host_y, pch = 21, bg = "darkgreen", col = "black", cex = 1.2)
  }
  
  x_vals <- sapply(lineage_path, function(pos) pos[2])
  y_vals <- sapply(lineage_path, function(pos) pos[1])
  lines(x_vals, y_vals, col = "black", lwd = 2, lty = 2)
  points(x_vals, y_vals, pch = 21, bg = "tomato", col = "black", cex = 1.2)
  
  text(0.5, 1, "", pos = 4, cex = 0.8)
}

compare_microbial_lineages_stochastic <- function() {
  par(mfrow = c(1, 3), mar = c(2, 2, 2, 1))
  
  # Host lineage as random walk
  host_positions <- numeric(20)
  host_positions[1] <- sample(3:8, 1)
  for (g in 2:20) {
    move <- sample(c(-1, 0, 1), 1)
    host_positions[g] <- min(max(host_positions[g - 1] + move, 1), 10)
  }
  host_lineage <- lapply(1:20, function(g) c(g, host_positions[g]))
  offspring_col <- host_positions[1]
  
  # Lineal: tracks host with rare deviations
  lineage_lineal <- lapply(1:20, function(g) {
    if (runif(1) < 0.15) c(g, sample(setdiff(1:10, host_positions[g]), 1))
    else c(g, host_positions[g])
  })
  
  microbial_lineage_grid_plot(lineage_path = lineage_lineal,
                              host_lineage = host_lineage,
                              title = "Lineal")
  
  # Non-lineal: random host trajectory, ends in offspring
  lineage_non_lineal <- lapply(20:2, function(g) c(g, sample(1:10, 1)))
  lineage_non_lineal <- append(lineage_non_lineal, list(c(1, offspring_col)))
  microbial_lineage_grid_plot(lineage_path = lineage_non_lineal,
                              host_lineage = host_lineage,
                              title = "Non-lineal")
  
  # Discordant: mostly environment, ends in offspring
  lineage_discordant <- list(c(1, offspring_col))
  for (g in 2:20) {
    if (runif(1) < 0.85) {
      lineage_discordant[[g]] <- c(g, 11)  # environment
    } else {
      lineage_discordant[[g]] <- c(g, sample(1:10, 1))  # host dip
    }
  }
  lineage_discordant <- rev(lineage_discordant)
  microbial_lineage_grid_plot(lineage_path = lineage_discordant,
                              host_lineage = host_lineage,
                              title = "Novel")
}

compare_microbial_lineages_stochastic()

