microbial_lineage_grid_plot <- function(n_hosts = 10, n_gens = 4,
                                        lineage_path,
                                        host_lineage = NULL,
                                        environment_col = n_hosts + 1,
                                        title = "Microbial Lineage Ancestry",
                                        y_axis_label = FALSE) {
  
  plot(NA, xlim = c(0.5, environment_col + 0.5), ylim = c(0.5, n_gens + 0.5),
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = if (y_axis_label) "Host Generation" else "",
       main = title, asp = 1)
  
  axis(2, at = 1:n_gens, labels = paste0("", 1:n_gens), las = 1, cex.axis = 0.6)
  
  for (gen in 1:n_gens) {
    points(1:n_hosts, rep(gen, n_hosts), pch = 21, bg = "gray90", col = "black")
    rect(environment_col - 0.4, gen - 0.4, environment_col + 0.4, gen + 0.4,
         col = rgb(0.9, 0.95, 1, 0.3), border = "lightblue")
  }
  text(environment_col, n_gens + 0.75, "Environment", cex = 0.8)
  
  if (!is.null(host_lineage)) {
    for (line in host_lineage) {
      host_x <- sapply(line, function(pos) pos[2])
      host_y <- sapply(line, function(pos) pos[1])
      lines(host_x, host_y, col = rgb(0, 0.5, 0, 0.4), lwd = 2, alpha=0.1)
      points(host_x, host_y, pch = 21, bg = "darkgreen", col = "black", cex = 1.2)
    }
  }
  
  x_vals <- sapply(lineage_path, function(pos) pos[2])
  y_vals <- sapply(lineage_path, function(pos) pos[1])
  lines(x_vals, y_vals, col = "black", lwd = 2, lty = 2)
  points(x_vals, y_vals, pch = 21, bg = "tomato", col = "black", cex = 1.2)
  
  text(0.5, 1, "", pos = 4, cex = 0.8)
}

compare_microbial_lineages_bifurcating <- function() {
  # set.seed(123)
  par(mfrow = c(1, 3), mar = c(2, 2, 2, 1))
  n_gens <- 4
  n_hosts <- 10
  env_col <- n_hosts + 1
  
  # Track hosts per generation
  host_nodes_by_gen <- list()
  host_edges <- list()
  
  # Start with offspring
  offspring_col <- sample(1:n_hosts, 1)
  host_nodes_by_gen[[1]] <- list(c(1, offspring_col))
  
  for (gen in 2:n_gens) {
    nodes_this_gen <- list()
    for (child in host_nodes_by_gen[[gen - 1]]) {
      parents <- sample(1:n_hosts, 2)
      for (p in parents) {
        parent_node <- c(gen, p)
        nodes_this_gen[[length(nodes_this_gen) + 1]] <- parent_node
        host_edges[[length(host_edges) + 1]] <- list(child, parent_node)
      }
    }
    host_nodes_by_gen[[gen]] <- nodes_this_gen
  }
  
  # Convert edges to path segments for plotting
  host_lineage <- lapply(host_edges, function(edge) list(edge[[1]], edge[[2]]))
  
  # --- Lineal Microbe: follows one ancestral path only ---
  # Build a single host lineage path back in time
  lineage_lineal <- list(c(1, offspring_col))
  current_node <- c(1, offspring_col)
  
  for (gen in 2:n_gens) {
    # Find parent(s) of current node in host_edges
    parent_edges <- Filter(function(edge) {
      all(edge[[1]] == current_node)
    }, host_edges)
    
    if (length(parent_edges) > 0) {
      chosen_parent <- sample(parent_edges, 1)[[1]][[2]]
      if (runif(1) < 0.1) {
        lineage_lineal[[gen]] <- c(gen, env_col)
      } else {
        lineage_lineal[[gen]] <- chosen_parent
        current_node <- chosen_parent
      }
    } else {
      # If no parent found (shouldn't happen), stay in environment
      lineage_lineal[[gen]] <- c(gen, env_col)
    }
  }
  
  microbial_lineage_grid_plot(lineage_path = lineage_lineal,
                              host_lineage = host_lineage,
                              title = "Lineal",
                              y_axis_label = TRUE)
  
  # --- Non-lineal Microbe ---
  lineage_non_lineal <- list(c(1, offspring_col))
  for (gen in 2:n_gens) {
    if (gen == 3) {
      lineage_non_lineal[[gen]] <- c(gen, env_col)
    } else {
      col <- if (runif(1) < 0.1) env_col else sample(1:n_hosts, 1)
      lineage_non_lineal[[gen]] <- c(gen, col)
    }
  }
  
  microbial_lineage_grid_plot(lineage_path = lineage_non_lineal,
                              host_lineage = host_lineage,
                              title = "Non-lineal")
  
  # --- Discordant Microbe ---
  lineage_discordant <- list(c(1, offspring_col))
  for (gen in 2:n_gens) {
    lineage_discordant[[gen]] <- c(gen, env_col)
  }
  lineage_discordant <- rev(lineage_discordant)
  
  microbial_lineage_grid_plot(lineage_path = lineage_discordant,
                              host_lineage = host_lineage,
                              title = "Novel")
}

compare_microbial_lineages_bifurcating()



compare_microbial_lineages_bifurcating <- function() {
  set.seed(123)
  par(mgp = c(1.75, 1, 0), mfrow = c(1, 3), mar = c(4, 4, 2, 1))
  par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))  # More space for y-axis label
  n_gens <- 4
  n_hosts <- 10
  env_col <- n_hosts + 1
  
  # Track hosts per generation
  host_nodes_by_gen <- list()
  host_edges <- list()
  
  # Start with offspring
  offspring_col <- sample(1:n_hosts, 1)
  host_nodes_by_gen[[1]] <- list(c(1, offspring_col))
  
  for (gen in 2:n_gens) {
    nodes_this_gen <- list()
    for (child in host_nodes_by_gen[[gen - 1]]) {
      parents <- sample(1:n_hosts, 2)
      for (p in parents) {
        parent_node <- c(gen, p)
        nodes_this_gen[[length(nodes_this_gen) + 1]] <- parent_node
        host_edges[[length(host_edges) + 1]] <- list(child, parent_node)
      }
    }
    host_nodes_by_gen[[gen]] <- nodes_this_gen
  }
  
  # Convert edges to path segments for plotting
  host_lineage <- lapply(host_edges, function(edge) list(edge[[1]], edge[[2]]))
  
  # --- Lineal Microbe: follows one ancestral path only ---
  lineage_lineal <- list(c(1, offspring_col))
  current_node <- c(1, offspring_col)
  
  for (gen in 2:n_gens) {
    parent_edges <- Filter(function(edge) all(edge[[1]] == current_node), host_edges)
    if (length(parent_edges) > 0) {
      chosen_parent <- sample(parent_edges, 1)[[1]][[2]]
      if (runif(1) < 0.1) {
        lineage_lineal[[gen]] <- c(gen, env_col)
      } else {
        lineage_lineal[[gen]] <- chosen_parent
        current_node <- chosen_parent
      }
    } else {
      lineage_lineal[[gen]] <- c(gen, env_col)
    }
  }
  
  microbial_lineage_grid_plot(lineage_path = lineage_lineal,
                              host_lineage = host_lineage,
                              title = "Lineal",
                              y_axis_label = TRUE)
  
  # --- Non-lineal Microbe: in environment in generation 3 ---
  lineage_non_lineal <- list(c(1, offspring_col))
  for (gen in 2:n_gens) {
    if (gen == 3) {
      lineage_non_lineal[[gen]] <- c(gen, env_col)
    } else {
      col <- if (runif(1) < 0.1) env_col else sample(1:n_hosts, 1)
      lineage_non_lineal[[gen]] <- c(gen, col)
    }
  }
  
  microbial_lineage_grid_plot(lineage_path = lineage_non_lineal,
                              host_lineage = host_lineage,
                              title = "Non-lineal",
                              y_axis_label = FALSE)
  
  # --- Discordant Microbe ---
  lineage_discordant <- list(c(1, offspring_col))
  for (gen in 2:n_gens) {
    lineage_discordant[[gen]] <- c(gen, env_col)
  }
  lineage_discordant <- rev(lineage_discordant)
  
  microbial_lineage_grid_plot(lineage_path = lineage_discordant,
                              host_lineage = host_lineage,
                              title = "Novel",
                              y_axis_label = FALSE)
}

compare_microbial_lineages_bifurcating()
