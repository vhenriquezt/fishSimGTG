#' Extract units from simulation result
#' @param simulation_result Output from runProjection
#' @return List of units for different metrics
#' @keywords internal

get_metric_units <- function(simulation_result) {
  lh <- simulation_result$LifeHistoryObj

  # Extract units from life history object
  length_units <- if(length(lh@L_units) > 0) lh@L_units else "cm"
  weight_units <- if(length(lh@Walpha_units) > 0) lh@Walpha_units else "g"

  # Define units for each metric
  units <- list(
    SB = weight_units,
    VB = weight_units,
    RB = weight_units,
    catchB = weight_units,
    catchN = "numbers",
    Ftotal = "year^-1",
    discB = weight_units,
    discN = "numbers",
    recN = "numbers",
    SPR = "proportion",
    TAC = weight_units
  )

  return(units)
}


#' Internal helper function to ensure observation data is in matrix format
#' @param simulation_result Output object from runProjection()
#' @return simulation_result with observation data in matrix format
#' @keywords internal

ensure_matrix_format <- function(simulation_result) {
  if(!is.null(simulation_result$HCR$decisionData) &&
     is.null(simulation_result$observation_matrices)) {

    obs_data <- simulation_result$HCR$decisionData
    has_indices <- any(grepl("^IDX_(Survey_|CPUE_)", names(obs_data)))
    has_catch <- any(grepl("^(true_catch|observed_catch|fleet_\\d+_(true_catch|observed_catch))$", names(obs_data)))
    has_lc <- any(grepl("^LC_.*_count_bin_\\d+$", names(obs_data)))

    if(has_indices || has_catch || has_lc) {
      return(adapt_observation_data(simulation_result, verbose = FALSE))
    }
  }
  return(simulation_result)
}

#Convert observation data from long vector format to matrix format
adapt_observation_data <- function(simulation_result,
                                   reshape_observations = TRUE,
                                   keep_original = FALSE,
                                   verbose = FALSE) {

  if(is.null(simulation_result$HCR$decisionData)) {
    if(verbose) cat("no observation data found in simulation result\n")
    return(simulation_result)
  }

  obs_data <- simulation_result$HCR$decisionData

  total_years <- unique(obs_data$years_total)[1]
  total_iterations <- unique(obs_data$iterations_total)[1]
  historical_end <- unique(obs_data$historical_end)[1]

  if(is.na(total_years) || is.na(total_iterations)) {
    if(verbose) cat("cannot determine simulation dimensions from observation data\n")
    return(simulation_result)
  }

  if(verbose) {
    cat("converting observation data to matrix format:\n")
    cat("  Years:", total_years, "\n")
    cat("  Iterations:", total_iterations, "\n")
    cat("  Historical end:", historical_end, "\n")
  }

  obs_columns <- find_observation_columns(obs_data, verbose)

  if(length(obs_columns) == 0) {
    if(verbose) cat("no observation columns found to reshape\n")
    return(simulation_result)
  }

  matrices <- list()

  for(col_name in obs_columns) {
    if(verbose) cat("  converting", col_name, "...")

    matrix_result <- convert_vector_to_matrix(
      obs_data, col_name, total_years, total_iterations, verbose
    )

    if(!is.null(matrix_result)) {
      matrices[[col_name]] <- matrix_result
      if(verbose) cat(" yes\n")
    } else {
      if(verbose) cat(" no (failed)\n")
    }
  }

  if(length(matrices) > 0) {
    if(reshape_observations) {
      simulation_result$observation_matrices <- matrices

      if(!keep_original) {
        for(col_name in names(matrices)) {
          if(!grepl("_(areas|indexYears|fleet_id|indextype|selectivity_|survey_timing|sample_size|total_catch)", col_name)) {
            simulation_result$HCR$decisionData[[col_name]] <- NULL
          }
        }
      }

      if(verbose) {
        cat("created", length(matrices), "observation matrices\n")
        cat("access via: simulation_result$observation_matrices$[column_name]\n")
      }
    }
  }

  return(simulation_result)
}

find_observation_columns <- function(obs_data, verbose = FALSE) {
  index_pattern <- "^IDX_(Survey_|CPUE_)"
  catch_pattern <- "^(true_catch|observed_catch|fleet_\\d+_(true_catch|observed_catch))$"
  lc_pattern <- "^LC_.*_count_bin_\\d+$"

  index_cols <- grep(index_pattern, names(obs_data), value = TRUE)
  catch_cols <- grep(catch_pattern, names(obs_data), value = TRUE)
  lc_cols <- grep(lc_pattern, names(obs_data), value = TRUE)

  exclude_pattern <- "_(areas|indexYears|fleet_id|indextype|selectivity_|survey_timing|sample_size|total_catch|cv_min|cv_max|obs_CV|R_t|n_areas|areas_included)$"

  all_obs_cols <- c(index_cols, catch_cols, lc_cols)
  all_obs_cols <- all_obs_cols[!grepl(exclude_pattern, all_obs_cols)]

  if(verbose && length(all_obs_cols) > 0) {
    cat("found observation columns to convert:\n")
    cat("  Indices:", length(index_cols), "\n")
    cat("  Catch:", length(catch_cols), "\n")
    cat("  Length comp:", length(lc_cols), "\n")
  }

  return(all_obs_cols)
}

convert_vector_to_matrix <- function(obs_data, col_name, total_years, total_iterations, verbose = FALSE) {
  obs_vector <- obs_data[[col_name]]

  if(is.null(obs_vector)) {
    if(verbose) cat("column", col_name, "not found\n")
    return(NULL)
  }

  valid_data <- obs_data[!is.na(obs_vector), ]

  if(nrow(valid_data) == 0) {
    if(verbose) cat("no valid data for", col_name, "\n")
    return(NULL)
  }

  obs_matrix <- matrix(NA, nrow = total_years, ncol = total_iterations)

  for(i in 1:nrow(valid_data)) {
    year_idx <- valid_data$j[i]
    iter_idx <- valid_data$k[i]
    value <- valid_data[[col_name]][i]

    if(year_idx >= 1 && year_idx <= total_years &&
       iter_idx >= 1 && iter_idx <= total_iterations &&
       !is.na(value)) {
      obs_matrix[year_idx, iter_idx] <- value
    }
  }

  rownames(obs_matrix) <- paste0("Year_", 1:total_years)
  colnames(obs_matrix) <- paste0("Iter_", 1:total_iterations)

  return(obs_matrix)
}

#' Plot population dynamics metrics
#'
#' Creates time series plots for various population dynamics metrics from fishSimGTG simulation results.
#' Supports both single-fleet and multifleet simulations with flexible visualization options.
#'
#' @param simulation_result Output object from \code{runProjection()}
#' @param metric Character. The metric to plot. One of: "SB", "VB", "RB", "catchB", "catchN", "Ftotal", "discB", "discN", "recN", "SPR"
#' @param areas Character "all" or numeric vector specifying which areas to plot. Default is "all"
#' @param iterations Character "all" or numeric vector specifying which simulation iterations to include. Default is "all"
#' @param show_median Logical. Whether to show median line across iterations. Default is TRUE
#' @param show_quantiles Logical. Whether to show 25th-75th percentile ribbon. Default is TRUE
#' @param show_individual Logical. Whether to show individual iteration lines. Default is FALSE
#' @param show_fleets Logical. Whether to show fleet-specific data in multifleet mode. Default is FALSE
#' @param color_palette Character vector of colors. If NULL, uses default colors
#' @param title Character. Custom plot title. If NULL, generates automatic title
#' @param save_plot Logical. Whether to save plot to file. Default is FALSE
#' @param filename Character. Filename for saved plot. If NULL, generates automatic filename
#' @param width Numeric. Plot width in inches for saved plot. Default is 12
#' @param height Numeric. Plot height in inches for saved plot. Default is 8
#'
#' @return A ggplot object
#' @export

plot_population_metric <- function(simulation_result,
                                   metric,
                                   areas = "all",
                                   iterations = "all",
                                   show_median = TRUE,
                                   show_quantiles = TRUE,
                                   show_individual = FALSE,
                                   show_fleets = FALSE,
                                   color_palette = NULL,
                                   title = NULL,
                                   save_plot = FALSE,
                                   filename = NULL,
                                   width = 12,
                                   height = 8) {


  #check if valid metric names exist
  valid_metrics <- c("SB", "VB", "RB", "catchB", "catchN", "Ftotal", "discB", "discN", "recN", "SPR")
  if(!metric %in% valid_metrics) {
    stop("metric must be one of: ", paste(valid_metrics, collapse = ", "))
  }
  dynamics <- simulation_result$dynamics
  is_multifleet <- !is.null(dynamics$multifleet)
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1

  #3 paths:
  #Population-level metrics (recN, SPR) - no areas dimension
  #Area-specific metrics with fleet
  #Area-specific metrics without fleet

  #dimensions
  if(metric %in% c("recN", "SPR")) {
    #population-level metrics (no areas)
    plot_data <- prepare_population_data(dynamics, metric, iterations, historical_end)
    p <- create_population_plot(plot_data, metric, show_median, show_quantiles,
                                show_individual, color_palette, title, historical_end, simulation_result)
  } else {
    #area-specific
    array_data <- dynamics[[metric]]
    total_areas <- dim(array_data)[3]

    #handle areas
    if(length(areas) == 1 && areas == "all") {
      areas <- 1:total_areas
    }

    if(any(areas > total_areas)) {
      stop("areas exceed available areas (", total_areas, ")")
    }

    #data preaparartion
    if(show_fleets && is_multifleet && metric %in% c("Ftotal", "catchB", "catchN", "discB", "discN")) {
      plot_data <- prepare_multifleet_data(dynamics, metric, areas, iterations, historical_end)
      p <- create_multifleet_plot(plot_data, metric, areas, show_median, show_quantiles,
                                  show_individual, color_palette, title, historical_end, simulation_result)
    } else {
      plot_data <- prepare_area_data(dynamics, metric, areas, iterations, historical_end)
      p <- create_area_plot(plot_data, metric, areas, show_median, show_quantiles,
                            show_individual, color_palette, title, historical_end,simulation_result)
    }
  }
  #save plot
  if(save_plot) {
    if(is.null(filename)) {
      area_suffix <- if(metric %in% c("recN", "SPR")) "population" else paste0("areas_", paste(areas, collapse = "_"))
      filename <- paste0(metric, "_", area_suffix, ".jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


prepare_population_data <- function(dynamics, metric, iterations, historical_end) {

  array_data <- dynamics[[metric]]
  years <- dim(array_data)[1]
  total_iterations <- dim(array_data)[2]

  #iterations
  if(length(iterations) == 1 && iterations == "all") {
    iterations <- 1:total_iterations
  }

  plot_data <- data.frame()

  for(iter in iterations) {
    iter_data <- data.frame(
      year = 1:years,
      user_year = 0:(years-1),
      value = array_data[, iter],
      iteration = iter,
      period = ifelse(1:years <= historical_end, "Historical", "Projection")
    )
    plot_data <- rbind(plot_data, iter_data)
  }

  return(plot_data)
}



prepare_area_data <- function(dynamics, metric, areas, iterations, historical_end) {

  array_data <- dynamics[[metric]]
  years <- dim(array_data)[1]
  total_iterations <- dim(array_data)[2]

  #iterations
  if(length(iterations) == 1 && iterations == "all") {
    iterations <- 1:total_iterations
  }

  plot_data <- data.frame()

  for(area in areas) {
    for(iter in iterations) {
      iter_data <- data.frame(
        year = 1:years,
        user_year = 0:(years-1),
        value = array_data[, iter, area],
        iteration = iter,
        area = paste("Area", area),
        period = ifelse(1:years <= historical_end, "Historical", "Projection")
      )
      plot_data <- rbind(plot_data, iter_data)
    }
  }

  return(plot_data)
}


prepare_multifleet_data <- function(dynamics, metric, areas, iterations, historical_end) {

  #array_data <- dynamics$multifleet[[paste0(metric, "_by_fleet")]]
  #to handle the correct array name for multifleet
  array_name <- paste0(metric, "_by_fleet")
  array_data <- dynamics$multifleet[[array_name]]

  #add a check
  # check if array exists
  if(is.null(array_data)) {
    stop("multifleet array '", array_name, "' it was not found in dynamics$multifleet")
  }

  years <- dim(array_data)[1]
  total_iterations <- dim(array_data)[2]
  nfleets <- dim(array_data)[4]

  #iterations
  if(length(iterations) == 1 && iterations == "all") {
    iterations <- 1:total_iterations
  }

  plot_data <- data.frame()

  for(area in areas) {
    for(fleet in 1:nfleets) {
      for(iter in iterations) {
        iter_data <- data.frame(
          year = 1:years,
          user_year = 0:(years-1),
          value = array_data[, iter, area, fleet],
          iteration = iter,
          area = paste("Area", area),
          fleet = paste("Fleet", fleet),
          area_fleet = paste("Area", area, "- Fleet", fleet),
          period = ifelse(1:years <= historical_end, "Historical", "Projection")
        )
        plot_data <- rbind(plot_data, iter_data)
      }
    }
  }

  return(plot_data)
}


create_population_plot <- function(plot_data, metric, show_median, show_quantiles,
                                   show_individual, color_palette, title, historical_end, simulation_result = NULL) {

  #calculate median and quantiles
  summary_data <- plot_data %>%
    group_by(user_year, period) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"                 #ungroup the data after summarizing
    )

  #generate title if not provided
  if(is.null(title)) {
    title <- paste(get_metric_label(metric), "(Total Population)")
  }

  if(is.null(color_palette)) {
    main_color <- get_metric_color(metric)
  } else {
    main_color <- color_palette[1]
  }

  p <- ggplot()

  #add iter individual lines
  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = user_year, y = value, group = iteration),
                       color = main_color, alpha = 0.3, size = 0.5)
  }

  #quantile ribbon
  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75),
                         fill = main_color, alpha = 0.3)
  }

  #add median
  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value),
                       color = main_color, size = 1.5)
  }

  #formatting
  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(title = title, x = "Year", y = get_metric_ylabel(metric, simulation_result)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  #SPR reference line
  if(metric == "SPR") {
    p <- p + geom_hline(yintercept = 0.3, linetype = "dotted", color = "red", alpha = 0.7)
  }

  return(p)
}


create_area_plot <- function(plot_data, metric, areas, show_median, show_quantiles,
                             show_individual, color_palette, title, historical_end, simulation_result = NULL) {

  #summary statistics
  summary_data <- plot_data %>%
    group_by(user_year, area, period) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  #generate title
  if(is.null(title)) {
    if(length(areas) == 1) {
      title <- paste(get_metric_label(metric), "- Area", areas)
    } else {
      title <- paste(get_metric_label(metric), "by Area")
    }
  }
# format colors
  if(is.null(color_palette)) {
    if(length(areas) == 1) {
      colors <- get_metric_color(metric)
    } else {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:length(areas)]
    }
  } else {
    colors <- color_palette[1:length(areas)]
  }

  p <- ggplot()

  # iteration lines
  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = user_year, y = value, group = interaction(area, iteration), color = area),
                       alpha = 0.3, size = 0.5)
  }

  #quantiles
  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75, fill = area),
                         alpha = 0.3)
  }

  #median
  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value, color = area),
                       size = 1.5)
  }

  # fixing for zero values - extend y axis under zero to show the line of zero
  # this fix the discB plot when discards are zero
  y_range <- range(plot_data$value, na.rm = TRUE)
  y_limits <- NULL

  if(all(y_range == 0)) {
    # All values are zero - create small range around zero
    y_limits <- c(-0.01, 0.01)
  } else if(y_range[1] == 0 && y_range[2] > 0) {
    # Min is zero - extend below to show zero line
    y_limits <- c(-y_range[2] * 0.05, y_range[2] * 1.05)
  }



  #formatting
  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(title = title, x = "Year", y = get_metric_ylabel(metric, simulation_result),
         color = "Area", fill = "Area")  +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = if(length(areas) > 1) "bottom" else "none"
    )


  # y-limits and add zero reference line for discard metrics
  if(!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }

  if(metric %in% c("discB", "discN")) {
    p <- p + geom_hline(yintercept = 0, linetype = "solid", color = "gray50", alpha = 0.7)
  }




  #using facets for multiple areas (e.g., more than 3 areas)
  if(length(areas) > 3) {
    p <- p + facet_wrap(~ area, scales = "free_y") +
      theme(legend.position = "none")
  }

  return(p)
}


create_multifleet_plot <- function(plot_data, metric, areas, show_median, show_quantiles,
                                   show_individual, color_palette, title, historical_end,simulation_result = NULL) {


#statistical summary (speratate stats for each area-fleet-period combination)
summary_data <- plot_data %>%
  group_by(user_year, area, fleet, area_fleet, period) %>%
  summarise(
    median_value = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
#title
if(is.null(title)) {
  title <- paste(get_metric_label(metric), "by Area and Fleet")
}
nfleets <- length(unique(plot_data$fleet))
if(is.null(color_palette)) {
  colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
} else {
  colors <- color_palette[1:nfleets]
}

p <- ggplot()
#iterations
if(show_individual) {
  p <- p + geom_line(data = plot_data,
                     aes(x = user_year, y = value, group = interaction(area_fleet, iteration), color = fleet),
                     alpha = 0.3, size = 0.5)
}
#quantiles
if(show_quantiles) {
  p <- p + geom_ribbon(data = summary_data,
                       aes(x = user_year, ymin = q25, ymax = q75, fill = fleet),
                       alpha = 0.3)
#median
  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value, color = fleet),
                       size = 1.5)
  }
}

#formatting
p <- p +
  geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = function(x) pretty(x, n = 6)) +
  #facet_wrap:creates separate panels for each area
  facet_wrap(~ area, scales = "free_y") +
  labs(title = title, x = "Year", y = get_metric_ylabel(metric, simulation_result),
       color = "Fleet", fill = "Fleet") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

return(p)
}



get_metric_label <- function(metric) {
  labels <- list(
    "SB" = "Spawning Biomass",
    "VB" = "Vulnerable Biomass",
    "RB" = "Retained Biomass",
    "catchB" = "Catch Biomass",
    "catchN" = "Catch Numbers",
    "Ftotal" = "Fishing Mortality",
    "discB" = "Discard Biomass",
    "discN" = "Discard Numbers",
    "recN" = "Recruitment",
    "SPR" = "Spawning Potential Ratio"
  )
  return(labels[[metric]])
}



get_metric_ylabel <- function(metric, simulation_result = NULL) {
  base_labels <- list(
    "SB" = "Spawning Biomass",
    "VB" = "Vulnerable Biomass",
    "RB" = "Retained Biomass",
    "catchB" = "Catch Biomass",
    "catchN" = "Catch Numbers",
    "Ftotal" = "Fishing Mortality (F)",
    "discB" = "Discard Biomass",
    "discN" = "Discard Numbers",
    "recN" = "Recruitment (Numbers)",
    "SPR" = "SPR",
    "TAC" = "TAC"
  )
  base_label <- base_labels[[metric]]

  #add units
  if(!is.null(simulation_result)) {
    units <- get_metric_units(simulation_result)
    unit <- units[[metric]]

    if(!is.null(unit) && unit != "proportion") {
      base_label <- paste0(base_label, " (", unit, ")")
    }
  }
  return(base_label)
}



get_metric_color <- function(metric) {
  colors <- list(
    "SB" = "steelblue",
    "VB" = "darkgreen",
    "RB" = "darkblue",
    "catchB" = "brown",
    "catchN" = "chocolate",
    "Ftotal" = "orange",
    "discB" = "red",
    "discN" = "darkred",
    "recN" = "cyan4",
    "SPR" = "purple"
  )
  return(colors[[metric]])
}


#' Plot spawning biomass
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_SB     <- function(result,areas="all", ...) plot_population_metric(result, "SB", areas = areas, ...)

#' Plot vulnerable biomass
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_VB     <- function(result,areas="all", ...) plot_population_metric(result, "VB", areas = areas, ...)

#' Plot total fishing mortality
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_Ftotal <- function(result,areas="all", ...) plot_population_metric(result, "Ftotal", areas = areas, ...)

#' Plot catch biomass
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_catchB <- function(result,areas="all", ...) plot_population_metric(result, "catchB", areas = areas, ...)

#' Plot catch numbers
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_catchN <- function(result,areas="all", ...) plot_population_metric(result, "catchN", areas = areas, ...)

#' Plot dicard catch biomass
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_discB  <- function(result,areas="all", ...) plot_population_metric(result, "discB", areas = areas, ...)

#' Plot dicard catch numbers
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_discN  <- function(result,areas="all", ...) plot_population_metric(result, "discN", areas = areas, ...)

#' Plot spawning potential ratio
#' @param result Simulation result from runProjection
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_SPR    <- function(result, ...) plot_population_metric(result, "SPR", ...)

#' Plot recruitment
#' @param result Simulation result from runProjection
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_recN   <- function(result, ...) plot_population_metric(result, "recN", ...)


#' Plot multifleet fishing mortality
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_Ftotal_multi <- function(result, areas = "all", ...) plot_population_metric(result, "Ftotal", areas = areas, show_fleets = TRUE, ...)

#' Plot multifleet catch biomass
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_catchB_multi <- function(result, areas = "all", ...) plot_population_metric(result, "catchB", areas = areas, show_fleets = TRUE, ...)

#' Plot multifleet catch numbers
#' @param result Simulation result from runProjection
#' @param areas Areas to plot ("all" or numeric vector)
#' @param ... Additional arguments passed to plot_population_metric
#' @export
plot_catchN_multi <- function(result, areas = "all", ...) plot_population_metric(result, "catchN", areas = areas, show_fleets = TRUE, ...)


#' Plot observation model indices
#'
#' Creates time series plots of survey indices and CPUE data from fishSimGTG observation models.
#'
#' @param simulation_result Output object from \code{runProjection()} containing observation data
#' @param index_pattern Character. Pattern to match index column names. Default is "IDX_"
#' @param show_median Logical. Whether to show median points. Default is TRUE
#' @param show_quantiles Logical. Whether to show quantile ranges. Default is TRUE
#' @param show_individual Logical. Whether to show individual iteration points. Default is FALSE
#' @param point_size Numeric. Size of points in plot. Default is 1.5
#' @param line_alpha Numeric. Transparency of individual iteration points (0-1). Default is 0.5
#' @param color_palette Character vector of colors. If NULL, uses default colors
#' @param title Character. Custom plot title. If NULL, generates automatic title
#' @param save_plot Logical. Whether to save plot to file. Default is FALSE
#' @param filename Character. Filename for saved plot. If NULL, generates automatic filename
#'
#' @return A ggplot object with faceted panels for each index
#' @export

plot_indices <- function(simulation_result,
                         index_pattern = "IDX_",
                         show_median = TRUE,
                         show_quantiles = TRUE,
                         show_individual = FALSE,
                         point_size = 1.5,
                         line_alpha = 0.5,
                         color_palette = NULL,
                         title = NULL,
                         save_plot = FALSE,
                         filename = NULL) {

  obs_data <- simulation_result$HCR$decisionData
  if(is.null(obs_data)) {
    stop("no observation data found in simulation result")
  }

  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1

  #find index columns
  #grep(): finds column names matching a pattern
  #!grepl(): excludes columns containing metadata suffix
  #value = TRUE: returns column names, not positions
  index_cols <- grep(paste0("^", index_pattern), names(obs_data), value = TRUE)
  index_cols <- index_cols[!grepl("_(areas|indexYears|fleet_id|indextype|selectivity_|survey_timing)", index_cols)]

  if(length(index_cols) == 0) {
    stop("No index columns found matching pattern: ", index_pattern)
  }

  #data prep
  plot_data <- prepare_index_data_with_gaps(obs_data, index_cols, historical_end)

  if(nrow(plot_data) == 0) {
    stop("no valid data for plotting")
  }

  #create plot
  p <- create_index_plot_with_gaps(plot_data, show_median, show_quantiles, show_individual,
                                   point_size, line_alpha, color_palette, title, historical_end)


  #save
  if(save_plot) {
    if(is.null(filename)) {
      pattern_clean <- gsub("_$", "", gsub("^IDX_", "", index_pattern))
      filename <- paste0("indices_", pattern_clean, ".jpeg")
    }
    ggsave(filename, p, width = 12, height = 8, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


prepare_index_data_with_gaps <- function(obs_data, index_cols, historical_end) {
  plot_data_list <- list()

  for(i in seq_along(index_cols)) {
    index_col <- index_cols[i]

    #get data
    areas_col <- paste0(index_col, "_areas")
    fleet_col <- paste0(index_col, "_fleet_id")

    areas_str <- if(areas_col %in% names(obs_data)) {
      unique_areas <- unique(obs_data[[areas_col]][!is.na(obs_data[[areas_col]])])
      if(length(unique_areas) > 0) unique_areas[1] else "Unknown"
    } else {
      "Unknown"
    }

    fleet_str <- if(fleet_col %in% names(obs_data)) {
      unique_fleets <- unique(obs_data[[fleet_col]][!is.na(obs_data[[fleet_col]])])
      if(length(unique_fleets) > 0 && !is.na(unique_fleets[1])) {
        paste("Fleet", unique_fleets[1])
      } else {
        NA
      }
    } else {
      NA
    }

    #panel name
    if(!is.na(fleet_str)) {
      panel_name <- paste(index_col, "-", fleet_str, "- Area(s)", areas_str)
    } else {
      panel_name <- paste(index_col, "- Area(s)", areas_str)
    }

    #extract data - !!sym(index_col for dynamic col reference
    #sym(index_col)  # converts string to symbol
    # "IDX_Survey_1" becomes IDX_Survey_1 (without quotes) - !! remove quote operator
    # tells dplyr: "treat this as an actual column name, not a string"

    index_data <- obs_data %>%
      filter(!is.na(!!sym(index_col))) %>%
      mutate(
        user_year = j - 1,
        iteration = k,
        value = !!sym(index_col),
        panel = panel_name,
        index_name = index_col,
        period = ifelse(j <= historical_end, "Historical", "Projection")
      ) %>%
      select(user_year, iteration, value, panel, index_name, period)

    if(nrow(index_data) > 0) {
      plot_data_list[[i]] <- index_data
    }
  }

  if(length(plot_data_list) == 0) {
    return(data.frame())
  }

  return(dplyr::bind_rows(plot_data_list))
}



create_index_plot_with_gaps <- function(plot_data, show_median, show_quantiles, show_individual,
                                        point_size, line_alpha, color_palette, title, historical_end) {

  #calculate stats
  summary_data <- plot_data %>%
    group_by(user_year, panel, index_name, period) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      n_obs = n(),
      .groups = "drop"
    )

  n_panels <- length(unique(plot_data$panel))
  if(is.null(color_palette)) {
    colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink", "cyan4", "red", "darkgray")[1:n_panels]
  } else {
    colors <- rep(color_palette, length.out = n_panels)
  }

  p <- ggplot()

  #iteration lines (connected only where data exists)
  if(show_individual) {
    p <- p + geom_point(data = plot_data,
                        aes(x = user_year, y = value),
                        alpha = line_alpha, size = point_size * 0.7, color = "lightblue")
  }

  #quantile
  if(show_quantiles) {
    p <- p + geom_pointrange(data = summary_data,
                             aes(x = user_year, y = median_value,
                                 ymin = q25, ymax = q75),
                             alpha = 0.6, size = 0.8, color = "gray60")
  }

  #median points and lines (only connect where data exists)
  if(show_median) {
    p <- p + geom_point(data = summary_data,
                        aes(x = user_year, y = median_value),
                        color = "black", size = point_size * 1.2)
  }

  #formatting
  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_wrap(~ panel, scales = "free_y") +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(
      title = if(is.null(title)) "Index/CPUE Time Series" else title,
      subtitle = "Points show observations (gaps = no data), black = median, gray = quantiles",
      x = "Year",
      y = "Index Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 9, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    )

  return(p)
}

#' Plot catch observations
#'
#' Creates time series plots comparing true catch vs observed catch from fishSimGTG observation models.
#'
#' @param simulation_result Output object from \code{runProjection()} containing catch observation data
#' @param catch_type Character. Type of catch to plot: "true", "observed", or "both". Default is "both"
#' @param fleet_specific Logical. Whether to create fleet-specific plots in multifleet mode. Default is TRUE
#' @param show_median Logical. Whether to show median lines. Default is TRUE
#' @param show_quantiles Logical. Whether to show quantile ribbons. Default is TRUE
#' @param show_individual Logical. Whether to show individual iteration lines. Default is FALSE
#' @param point_size Numeric. Size of points in plot. Default is 1.5
#' @param title Character. Custom plot title. If NULL, generates automatic title
#' @param save_plot Logical. Whether to save plot to file. Default is FALSE
#' @param filename Character. Filename for saved plot. If NULL, generates automatic filename
#'
#' @return A ggplot object showing catch time series
#' @export

plot_catch_observations <- function(simulation_result,
                                    catch_type = "both",
                                    areas = "all",
                                    fleet_specific = TRUE,
                                    show_median = TRUE,
                                    show_quantiles = TRUE,
                                    show_individual = FALSE,
                                    point_size = 1.5,
                                    title = NULL,
                                    save_plot = FALSE,
                                    filename = NULL) {

  obs_data <- simulation_result$HCR$decisionData
  if(is.null(obs_data)) {
    stop("no observation data found in simulation result")
  }

  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  total_areas <- unique(obs_data$total_areas)[1]

  #area filter
  if(length(areas) == 1 && areas == "all") {
    areas_to_plot <- 1:total_areas
  } else {
    areas_to_plot <- areas
  }


  #detect if multifleet
  has_multifleet <- any(grepl("fleet_\\d+_", names(obs_data)))

  if(has_multifleet && fleet_specific) {
    plot_data <- prepare_multifleet_catch_data(obs_data, catch_type, historical_end, areas_to_plot)
    p <- create_multifleet_catch_plot(plot_data, catch_type, show_median, show_quantiles,
                                      show_individual, point_size, title, historical_end)
  } else {
    plot_data <- prepare_single_catch_data(obs_data, catch_type, historical_end,areas_to_plot)
    p <- create_single_catch_plot(plot_data, catch_type, show_median, show_quantiles,
                                  show_individual, point_size, title, historical_end)
  }

  if(save_plot) {
    if(is.null(filename)) {
      mode <- if(has_multifleet && fleet_specific) "multifleet" else "single"
      area_suffix <- paste0("areas_", paste(areas_to_plot, collapse = "_"))
      filename <- paste0("catch_observations_", mode, "_", catch_type, ".jpeg")
    }
    ggsave(filename, p, width = 12, height = 8, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


prepare_multifleet_catch_data <- function(obs_data, catch_type, historical_end,areas_to_plot = NULL) {

  #find fleet catch columns
  if(catch_type %in% c("true", "both")) {
    true_cols <- grep("fleet_\\d+_true_catch$", names(obs_data), value = TRUE)
  }
  if(catch_type %in% c("observed", "both")) {
    obs_cols <- grep("fleet_\\d+_observed_catch$", names(obs_data), value = TRUE)
  }

  plot_data_list <- list()

  #extract fleet numbers
  if(catch_type == "true") {
    fleet_nums <- unique(gsub("fleet_(\\d+)_true_catch", "\\1", true_cols))
  } else if(catch_type == "observed") {
    fleet_nums <- unique(gsub("fleet_(\\d+)_observed_catch", "\\1", obs_cols))
  } else {
    fleet_nums <- unique(c(
      gsub("fleet_(\\d+)_true_catch", "\\1", true_cols),
      gsub("fleet_(\\d+)_observed_catch", "\\1", obs_cols)
    ))
  }

  for(fleet_num in fleet_nums) {
    true_col <- paste0("fleet_", fleet_num, "_true_catch")
    obs_col <- paste0("fleet_", fleet_num, "_observed_catch")
    areas_col <- paste0("fleet_", fleet_num, "_areas_included")

    #get areas info
    areas_str <- if(areas_col %in% names(obs_data)) {
      unique_areas <- unique(obs_data[[areas_col]][!is.na(obs_data[[areas_col]])])
      if(length(unique_areas) > 0) unique_areas[1] else "Unknown"
    } else {
      "Unknown"
    }

    fleet_data <- data.frame()

    if(catch_type %in% c("true", "both") && true_col %in% names(obs_data)) {
      true_data <- obs_data %>%
        #true_col contains a string ("fleet_1_true_catch") - convert to fleet_1_true_catch
        filter(!is.na(!!sym(true_col))) %>%
        mutate(
          user_year = j - 1,
          iteration = k,
          value = !!sym(true_col),
          catch_type = "True Catch",
          fleet = paste("Fleet", fleet_num),
          panel = paste("Fleet", fleet_num, "- Area(s)", areas_str),
          period = ifelse(j <= historical_end, "Historical", "Projection")
        ) %>%
        select(user_year, iteration, value, catch_type, fleet, panel, period)

      fleet_data <- rbind(fleet_data, true_data)
    }

    if(catch_type %in% c("observed", "both") && obs_col %in% names(obs_data)) {
      obs_data_sub <- obs_data %>%
        filter(!is.na(!!sym(obs_col))) %>%
        mutate(
          user_year = j - 1,
          iteration = k,
          value = !!sym(obs_col),
          catch_type = "Observed Catch",
          fleet = paste("Fleet", fleet_num),
          panel = paste("Fleet", fleet_num, "- Area(s)", areas_str),
          period = ifelse(j <= historical_end, "Historical", "Projection")
        ) %>%
        select(user_year, iteration, value, catch_type, fleet, panel, period)

      fleet_data <- rbind(fleet_data, obs_data_sub)
    }

    if(nrow(fleet_data) > 0) {
      plot_data_list[[length(plot_data_list) + 1]] <- fleet_data
    }
  }

  if(length(plot_data_list) == 0) {
    return(data.frame())
  }

  return(dplyr::bind_rows(plot_data_list))
}


prepare_single_catch_data <- function(obs_data, catch_type, historical_end, areas_to_plot) {

  plot_data <- data.frame()

  #if no area filtering, use total catch columns
  if(is.null(areas_to_plot)) {

  if(catch_type %in% c("true", "both") && "true_catch" %in% names(obs_data)) {
    true_data <- obs_data %>%
      filter(!is.na(true_catch)) %>%
      mutate(
        user_year = j - 1,
        iteration = k,
        value = true_catch,
        catch_type = "True Catch",
        panel = "Total Catch",
        period = ifelse(j <= historical_end, "Historical", "Projection")
      ) %>%
      select(user_year, iteration, value, catch_type, panel, period)

    plot_data <- rbind(plot_data, true_data)
  }

  if(catch_type %in% c("observed", "both") && "observed_catch" %in% names(obs_data)) {
    obs_data_sub <- obs_data %>%
      filter(!is.na(observed_catch)) %>%
      mutate(
        user_year = j - 1,
        iteration = k,
        value = observed_catch,
        catch_type = "Observed Catch",
        panel = "Total Catch",
        period = ifelse(j <= historical_end, "Historical", "Projection")
      ) %>%
      select(user_year, iteration, value, catch_type, panel, period)

    plot_data <- rbind(plot_data, obs_data_sub)
  }

  } else {
    # area-specific catch
    for(area in areas_to_plot) {
      if(catch_type %in% c("true", "both")) {
        true_col <- paste0("true_catch_area_", area)
        if(true_col %in% names(obs_data)) {
          true_data <- obs_data %>%
            filter(!is.na(!!sym(true_col))) %>%
            mutate(
              user_year = j - 1,
              iteration = k,
              value = !!sym(true_col),
              catch_type = "True Catch",
              panel = paste("Area", area),
              period = ifelse(j <= historical_end, "Historical", "Projection")
            ) %>%
            select(user_year, iteration, value, catch_type, panel, period)

          plot_data <- rbind(plot_data, true_data)
        }
      }

      if(catch_type %in% c("observed", "both")) {
        obs_col <- paste0("observed_catch_area_", area)
        if(obs_col %in% names(obs_data)) {
          obs_data_sub <- obs_data %>%
            filter(!is.na(!!sym(obs_col))) %>%
            mutate(
              user_year = j - 1,
              iteration = k,
              value = !!sym(obs_col),
              catch_type = "Observed Catch",
              panel = paste("Area", area),
              period = ifelse(j <= historical_end, "Historical", "Projection")
            ) %>%
            select(user_year, iteration, value, catch_type, panel, period)

          plot_data <- rbind(plot_data, obs_data_sub)
        }
      }
    }
  }

  return(plot_data)
}


create_multifleet_catch_plot <- function(plot_data, catch_type, show_median, show_quantiles,
                                         show_individual, point_size, title, historical_end) {

  if(nrow(plot_data) == 0) return(NULL)

  #calculate stats
  summary_data <- plot_data %>%
    group_by(user_year, panel, catch_type, period) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  t
  p <- ggplot()

  #iterations
  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = user_year, y = value, color = catch_type,
                           group = interaction(catch_type, iteration)),
                       alpha = 0.4, size = 0.5)
    p <- p + geom_point(data = plot_data,
                        aes(x = user_year, y = value, color = catch_type),
                        alpha = 0.6, size = point_size * 0.8)
  }

  #quantiles
  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75, fill = catch_type),
                         alpha = 0.3)
  }

  #median
  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value, color = catch_type),
                       size = 1.5)
    p <- p + geom_point(data = summary_data,
                        aes(x = user_year, y = median_value, color = catch_type),
                        size = point_size)
  }

  #formatting
  colors <- c("True Catch" = "darkblue", "Observed Catch" = "orange")

  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    facet_wrap(~ panel, scales = "free_y") +
    scale_x_continuous(breaks = function(x) pretty(x, n = 6)) +
    labs(
      title = if(is.null(title)) "Multifleet Catch Observations" else title,
      subtitle = "Fleet-specific catch data with observation model effects",
      x = "Year", y = "Catch", color = "Catch Type", fill = "Catch Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  return(p)
}



create_single_catch_plot <- function(plot_data, catch_type, show_median, show_quantiles,
                                     show_individual, point_size, title, historical_end) {

  if(nrow(plot_data) == 0) return(NULL)

  #calculate stats
  summary_data <- plot_data %>%
    group_by(user_year,panel, catch_type, period) %>%   #panel contains area info
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  colors <- c("True Catch" = "darkblue", "Observed Catch" = "orange")

  p <- ggplot()

  #iterations
  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = user_year, y = value, color = catch_type,
                           group = interaction(catch_type, iteration)),
                       alpha = 0.4, size = 0.5)
    p <- p + geom_point(data = plot_data,
                        aes(x = user_year, y = value, color = catch_type),
                        alpha = 0.6, size = point_size * 0.8)
  }

  #quantiles
  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75, fill = catch_type),
                         alpha = 0.3)
  }

  #median
  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value, color = catch_type),
                       size = 1.5)
    p <- p + geom_point(data = summary_data,
                        aes(x = user_year, y = median_value, color = catch_type),
                        size = point_size)
  }

  #check if multiple panels (areas)
  has_multiple_panels <- length(unique(plot_data$panel)) > 1

  #formatting
  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(
      title = if(is.null(title)) "Catch Observations" else title,
      subtitle = "True vs observed catch showing observation model effects",
      x = "Year", y = "Catch", color = "Catch Type", fill = "Catch Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  #add faceting if multiple areas
  if(has_multiple_panels) {
    p <- p + facet_wrap(~ panel, scales = "free_y")
  }

  return(p)
}



#' Plot length composition by area
#'
#' Creates length composition plots from fishSimGTG observation models with filtering by area, fleet, and time.
#'
#' @param simulation_result Output object from \code{runProjection()} containing length composition data
#' @param program_pattern Character. Pattern to match LC program names. Default is "LC_"
#' @param area_filter Character "all" or numeric vector specifying areas to include. Default is "all"
#' @param fleet_filter Character "all" or numeric vector specifying fleets to include. Default is "all"
#' @param years_to_plot Character "all", "auto", or numeric vector of years to plot. Default is "all"
#' @param max_programs Numeric. Maximum number of programs to plot. Default is 8
#' @param show_individual Logical. Whether to show individual iteration lines. Default is FALSE
#' @param show_median Logical. Whether to show median bars. Default is TRUE
#' @param separate_by_area Logical. Whether to separate plots by area. Default is TRUE
#' @param title Character. Custom plot title. If NULL, generates automatic title
#' @param save_plot Logical. Whether to save plot to file. Default is FALSE
#' @param filename Character. Filename for saved plot. If NULL, generates automatic filename
#' @param auto_display Logical. Whether to automatically display multi-panel plots. Default is TRUE
#'
#' @return A ggplot object or grid of plots showing length composition histograms
#' @export

plot_length_composition_by_area <- function(simulation_result,
                                            program_pattern = "LC_",
                                            area_filter = "all",
                                            fleet_filter = "all",
                                            years_to_plot = "all",
                                            max_programs = 8,
                                            show_individual = FALSE,
                                            show_median = TRUE,
                                            separate_by_area = TRUE,
                                            title = NULL,
                                            save_plot = FALSE,
                                            filename = NULL,
                                            auto_display = TRUE) {

  obs_data <- simulation_result$HCR$decisionData
  if(is.null(obs_data)) {
    stop("no observation data found in simulation result")
  }

  #length bin info
  n_bins <- unique(obs_data$n_length_bins)[1]
  if(is.na(n_bins) || n_bins == 0) {
    stop("no length composition data found")
  }

  length_bin_width <- unique(obs_data$length_bin_width)[1]

  #get total areas in simulation
  total_areas <- unique(obs_data$total_areas)[1]

  #area filtering
  if(length(area_filter) == 1 && area_filter == "all") {
    areas_to_plot <- 1:total_areas
  } else {
    areas_to_plot <- area_filter
  }

  #find LC programs matching pattern
  sample_size_cols <- grep(paste0("^", program_pattern, ".*_sample_size$"), names(obs_data), value = TRUE)
  program_names <- gsub("_sample_size$", "", sample_size_cols)

  if(length(program_names) == 0) {
    stop("no LC programs found matching pattern: ", program_pattern)
  }

  #filter programs by fleet if specified
  if(length(fleet_filter) == 1 && fleet_filter != "all") {
    fleet_programs <- grep(paste0("_Fleet_", fleet_filter), program_names, value = TRUE)
    if(length(fleet_programs) > 0) {
      program_names <- fleet_programs
    }
  } else if(length(fleet_filter) > 1) {
    fleet_pattern <- paste0("_Fleet_(", paste(fleet_filter, collapse = "|"), ")")
    fleet_programs <- grep(fleet_pattern, program_names, value = TRUE)
    if(length(fleet_programs) > 0) {
      program_names <- fleet_programs
    }
  }

  #limit number of programs
  if(length(program_names) > max_programs) {
    program_names <- program_names[1:max_programs]
    cat("limiting to first", max_programs, "programs:", paste(program_names, collapse = ", "), "\n")
  }

  plot_list <- list()

  for(program_name in program_names) {
    cat("processing program:", program_name, "\n")

    #get program areas info
    areas_col <- paste0(program_name, "_areas")
    program_areas_str <- if(areas_col %in% names(obs_data)) {
      unique_areas_str <- unique(obs_data[[areas_col]][!is.na(obs_data[[areas_col]])])
      if(length(unique_areas_str) > 0) unique_areas_str[1] else "Unknown"
    } else {
      "Unknown"
    }

    #parse program areas (convert "1_2" to c(1,2))
    if(program_areas_str != "Unknown") {
      program_areas <- as.numeric(unlist(strsplit(program_areas_str, "_")))
    } else {
      program_areas <- 1:total_areas  # default to all areas
    }

    #filter areas for this program
    relevant_areas <- intersect(program_areas, areas_to_plot)

    if(length(relevant_areas) == 0) {
      cat("no relevant areas for program:", program_name, "\n")
      next
    }

    #years with data for this program
    sample_size_col <- paste0(program_name, "_sample_size")
    available_years <- sort(unique(obs_data$j[!is.na(obs_data[[sample_size_col]])]))

    if(length(available_years) == 0) {
      cat("No data for program:", program_name, "\n")
      next
    }

    #select years to plot
    if(length(years_to_plot) == 1 && years_to_plot == "all") {
      selected_years <- available_years
    } else if(length(years_to_plot) == 1 && years_to_plot == "auto") {
      if(length(available_years) <= 8) {
        selected_years <- available_years
      } else {
        indices <- round(seq(1, length(available_years), length.out = 8))
        selected_years <- available_years[indices]
      }
    } else {
      specified_sim_years <- years_to_plot + 1
      selected_years <- intersect(specified_sim_years, available_years)
    }

    if(length(selected_years) == 0) {
      next
    }

    #prepare data with area information
    lc_data <- prepare_length_comp_data_with_areas(obs_data, program_name, selected_years,
                                                   n_bins, length_bin_width,
                                                   program_areas, relevant_areas)

    if(nrow(lc_data) > 0) {
      #create plot for this program
      p <- create_length_comp_plot_with_areas(lc_data, program_name, show_individual,
                                              show_median, length_bin_width,
                                              separate_by_area, relevant_areas)
      plot_list[[program_name]] <- p
    }
  }

  if(length(plot_list) == 0) {
    stop("no valid length composition data found for plotting")
  }

  #combine plots and display
  if(length(plot_list) == 1) {
    final_plot <- plot_list[[1]]
  } else {
    final_plot <- gridExtra::arrangeGrob(grobs = plot_list, ncol = 2)

    if(auto_display) {
      gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
    }
  }

  #save
  if(save_plot) {
    if(is.null(filename)) {
      pattern_clean <- gsub("_$", "", gsub("^LC_", "", program_pattern))
      area_suffix <- if(length(areas_to_plot) == 1) paste0("_area", areas_to_plot) else "_all_areas"
      filename <- paste0("length_composition_", pattern_clean, area_suffix, ".png")
    }

    if(length(plot_list) == 1) {
      ggsave(filename, final_plot, width = 12, height = 8, dpi = 300)
    } else {
      ggsave(filename, final_plot, width = 16, height = 12, dpi = 300)
    }
    cat("Plot saved as:", filename, "\n")
  }

  return(final_plot)
}


prepare_length_comp_data_with_areas <- function(obs_data, program_name, selected_years,
                                                n_bins, length_bin_width,
                                                program_areas, relevant_areas) {

  plot_data <- data.frame()

  for(sim_year in selected_years) {
    year_data <- obs_data[obs_data$j == sim_year, ]

    if(nrow(year_data) == 0) next

    for(iter in unique(year_data$k)) {
      iter_data <- year_data[year_data$k == iter, ]
      if(nrow(iter_data) == 0) next

      #extract length composition for this iteration
      numbers <- numeric(n_bins)
      has_data <- FALSE

      #extracts length composition data bin by bin
      for(bin in 1:n_bins) {
        col_name <- paste0(program_name, "_count_bin_", bin)
        if(col_name %in% names(iter_data) && !is.na(iter_data[[col_name]][1])) {
          numbers[bin] <- iter_data[[col_name]][1]
          if(numbers[bin] > 0) has_data <- TRUE
        }
      }

      if(has_data && sum(numbers) > 0) {
        #create length bins (midpoints)
        length_bins <- seq(length_bin_width/2,
                           n_bins * length_bin_width - length_bin_width/2,
                           by = length_bin_width)

        #determine area label based on program configuration
        if(length(program_areas) == 1) {
          area_label <- paste("Area", program_areas[1])
        } else {
          area_label <- paste("Areas", paste(program_areas, collapse = "+"))
        }

        iteration_df <- data.frame(
          length_bin = length_bins,
          numbers = numbers,
          user_year = sim_year - 1,
          iteration = iter,
          year_label = paste("Year", sim_year - 1),
          program = gsub("^LC_", "", program_name),
          area_label = area_label,
          program_areas = paste(program_areas, collapse = "_")
        )

        plot_data <- rbind(plot_data, iteration_df)
      }
    }
  }

  return(plot_data)
}

create_length_comp_plot_with_areas <- function(lc_data, program_name, show_individual,
                                               show_median, length_bin_width,
                                               separate_by_area, relevant_areas) {

  if(nrow(lc_data) == 0) {
    return(NULL)
  }

  #calculate median across iterations
  if(show_median) {
    median_data <- lc_data %>%
      group_by(user_year, length_bin, year_label, program, area_label) %>%
      summarise(median_numbers = median(numbers, na.rm = TRUE), .groups = "drop")
  }

  #order year labels properly
  year_order <- sort(unique(lc_data$user_year))
  lc_data$year_label <- factor(lc_data$year_label,
                               levels = paste("Year", year_order))

  if(show_median) {
    median_data$year_label <- factor(median_data$year_label,
                                     levels = paste("Year", year_order))
  }

  p <- ggplot()

  #individual iterations
  if(show_individual) {
    p <- p + geom_line(data = lc_data,
                       aes(x = length_bin, y = numbers,
                           group = interaction(user_year, iteration, area_label)),
                       color = "lightblue", alpha = 0.3, size = 0.3)
  }

  #median bars
  if(show_median) {
    p <- p + geom_col(data = median_data,
                      aes(x = length_bin, y = median_numbers),
                      fill = "darkblue", alpha = 0.7, width = length_bin_width * 0.8)
  }

  #create title with area information
  clean_program_name <- gsub("^LC_", "", program_name)
  area_info <- unique(lc_data$area_label)[1]
  plot_title <- paste("Length Composition -", clean_program_name, "-", area_info)

  #formatting
  p <- p +
    facet_wrap(~ year_label, scales = "free_y") +
    scale_x_continuous(breaks = function(x) pretty(x, n = 6)) +
    labs(
      title = plot_title,
      subtitle = if(show_individual) "Light blue = individual iterations, dark blue = median" else "Dark blue bars = median across iterations",
      x = "Length (cm)",
      y = "Frequency (Numbers)"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}


#' Plot survey indices
#'
#' Creates time series plots of fishery-independent survey indices from fishSimGTG observation models.
#' This is a wrapper function that filters for survey indices (IDX_Survey pattern).
#'
#' @param simulation_result Output object from \code{runProjection()} containing observation data
#' @param ... Additional arguments passed to \code{plot_indices()}
#' @return A ggplot object with faceted panels for each survey index
#' @export
#' @seealso \code{\link{plot_indices}}
plot_survey_indices <- function(simulation_result, ...) {
  plot_indices(simulation_result, index_pattern = "IDX_Survey",
               title = "Survey Indices (Fishery Independent)", ...)
}

#' Plot CPUE indices
#'
#' Creates time series plots of fishery-dependent CPUE indices from fishSimGTG observation models.
#' This is a wrapper function that filters for CPUE indices (IDX_CPUE pattern).
#'
#' @param simulation_result Output object from \code{runProjection()} containing observation data
#' @param ... Additional arguments passed to \code{plot_indices()}
#' @return A ggplot object with faceted panels for each CPUE index
#' @export
#' @seealso \code{\link{plot_indices}}
plot_cpue_indices <- function(simulation_result, ...) {
  plot_indices(simulation_result, index_pattern = "IDX_CPUE",
               title = "CPUE Indices (Fishery Dependent)", ...)
}

#' Plot all indices
#'
#' Creates time series plots of all indices (both survey and CPUE) from fishSimGTG observation models.
#' This is a wrapper function that includes all index types.
#'
#' @param simulation_result Output object from \code{runProjection()} containing observation data
#' @param ... Additional arguments passed to \code{plot_indices()}
#' @return A ggplot object with faceted panels for each index
#' @export
#' @seealso \code{\link{plot_indices}}
plot_all_indices <- function(simulation_result, ...) {
  plot_indices(simulation_result, index_pattern = "IDX_",
               title = "All Indices (Survey + CPUE)", ...)
}


#' Plot both true and observed catch
#'
#' Creates time series plots comparing true catch vs observed catch from fishSimGTG catch observation models.
#' Shows both true and observed values on the same plot.
#'
#' @param simulation_result Output object from \code{runProjection()} containing catch observation data
#' @param areas Character "all" or numeric vector of areas. Default is "all"
#' @param ... Additional arguments passed to \code{plot_catch_observations()}
#' @return A ggplot object showing catch time series
#' @export
#' @seealso \code{\link{plot_catch_observations}}
plot_catch_observations_both <- function(simulation_result, ...) {
  plot_catch_observations(simulation_result, catch_type = "both", ...)
}

#' Plot true catch only
#'
#' Creates time series plots of true catch from fishSimGTG catch observation models.
#' Shows only the true (simulated) catch values.
#'
#' @param simulation_result Output object from \code{runProjection()} containing catch observation data
#' @param ... Additional arguments passed to \code{plot_catch_observations()}
#' @return A ggplot object showing true catch time series
#' @export
#' @seealso \code{\link{plot_catch_observations}}
plot_catch_observations_true <- function(simulation_result, ...) {
  plot_catch_observations(simulation_result, catch_type = "true", ...)
}

#' Plot observed catch only
#'
#' Creates time series plots of observed catch from fishSimGTG catch observation models.
#' Shows only the observed catch values with observation model effects.
#'
#' @param simulation_result Output object from \code{runProjection()} containing catch observation data
#' @param ... Additional arguments passed to \code{plot_catch_observations()}
#' @return A ggplot object showing observed catch time series
#' @export
#' @seealso \code{\link{plot_catch_observations}}
plot_catch_observations_observed <- function(simulation_result, ...) {
  plot_catch_observations(simulation_result, catch_type = "observed", ...)
}

#' Plot catch observations (single fleet mode)
#'
#' Creates time series plots of catch observations using single fleet aggregation.
#' Forces single fleet display even in multifleet simulations.
#'
#' @param simulation_result Output object from \code{runProjection()} containing catch observation data
#' @param ... Additional arguments passed to \code{plot_catch_observations()}
#' @return A ggplot object showing aggregated catch time series
#' @export
#' @seealso \code{\link{plot_catch_observations}}
plot_catch_observations_single <- function(simulation_result, ...) {
  plot_catch_observations(simulation_result, fleet_specific = FALSE, ...)
}

#' Plot catch observations (multifleet mode)
#'
#' Creates time series plots of catch observations with fleet-specific breakdown.
#' Shows separate panels for each fleet in multifleet simulations.
#'
#' @param simulation_result Output object from \code{runProjection()} containing catch observation data
#' @param ... Additional arguments passed to \code{plot_catch_observations()}
#' @return A ggplot object showing fleet-specific catch time series
#' @export
#' @seealso \code{\link{plot_catch_observations}}
plot_catch_observations_multifleet <- function(simulation_result, ...) {
  plot_catch_observations(simulation_result, fleet_specific = TRUE, ...)
}



#' Plot fishery length composition by area and fleet
#'
#' Creates length composition plots from fishery-dependent length composition data.
#' Allows filtering by specific areas and fleets.
#'
#' @param result Output object from \code{runProjection()} containing length composition data
#' @param areas Character "all" or numeric vector specifying areas to include. Default is "all"
#' @param fleets Character "all" or numeric vector specifying fleets to include. Default is "all"
#' @param ... Additional arguments passed to \code{plot_length_composition_by_area()}
#' @return A ggplot object or grid of plots showing length composition histograms
#' @export
#' @seealso \code{\link{plot_length_composition_by_area}}
plot_fishery_length_comp <- function(result, areas = "all", fleets = "all", ...) {
  plot_length_composition_by_area(result, program_pattern = "LC_Fishery",
                                  area_filter = areas, fleet_filter = fleets, ...)
}

#' Plot survey length composition by area
#'
#' Creates length composition plots from fishery-independent survey length composition data.
#' Allows filtering by specific areas.
#'
#' @param result Output object from \code{runProjection()} containing length composition data
#' @param areas Character "all" or numeric vector specifying areas to include. Default is "all"
#' @param ... Additional arguments passed to \code{plot_length_composition_by_area()}
#' @return A ggplot object or grid of plots showing length composition histograms
#' @export
#' @seealso \code{\link{plot_length_composition_by_area}}
plot_survey_length_comp <- function(result, areas = "all", ...) {
  plot_length_composition_by_area(result, program_pattern = "LC_Survey",
                                  area_filter = areas, ...)
}

#' Plot TAC (Total Allowable Catch) by area
#'
#' @param simulation_result Output from runProjection
#' @param areas Character "all" or numeric vector
#' @param show_fleets Logical. Show fleet breakdown in multifleet
#' @param ... Additional plotting arguments
#' @export
plot_TAC <- function(simulation_result,
                     areas = "all",
                     show_fleets = FALSE,
                     show_median = TRUE,
                     show_quantiles = TRUE,
                     show_individual = FALSE,
                     title = NULL,
                     save_plot = FALSE,
                     filename = NULL) {

  tac_data <- simulation_result$HCR$decisionAnnual

  if(is.null(tac_data) || !"TAC" %in% names(tac_data)) {
    stop("no TAC data found in simulation result")
  }

  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  units <- get_metric_units(simulation_result)

  # Get total areas
  total_areas <- max(tac_data$area)

  # Handle areas
  if(length(areas) == 1 && areas == "all") {
    areas <- 1:total_areas
  }

  # Filter areas
  tac_data <- tac_data[tac_data$area %in% areas, ]

  # Check for fleet column
  has_fleets <- "fleet" %in% names(tac_data)

  if(has_fleets && show_fleets) {
    p <- create_TAC_plot_multifleet(tac_data, areas, show_median, show_quantiles,
                                    show_individual, title, historical_end, units)
  } else {
    p <- create_TAC_plot_single(tac_data, areas, show_median, show_quantiles,
                                show_individual, title, historical_end, units)
  }

  if(save_plot) {
    if(is.null(filename)) {
      area_suffix <- paste0("areas_", paste(areas, collapse = "_"))
      filename <- paste0("TAC_", area_suffix, ".jpeg")
    }
    ggsave(filename, p, width = 12, height = 8, dpi = 300)
  }

  return(p)
}

create_TAC_plot_single <- function(tac_data, areas, show_median, show_quantiles,
                                   show_individual, title, historical_end, units) {

  # Prepare data
  plot_data <- tac_data %>%
    mutate(
      user_year = year - 1,
      area_label = paste("Area", area),
      period = ifelse(year <= historical_end, "Historical", "Projection")
    )

  # Summary statistics
  summary_data <- plot_data %>%
    group_by(user_year, area, area_label, period) %>%
    summarise(
      median_TAC = median(TAC, na.rm = TRUE),
      q25 = quantile(TAC, 0.25, na.rm = TRUE),
      q75 = quantile(TAC, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  colors <- c("steelblue", "darkgreen", "orange", "purple")[1:length(areas)]

  p <- ggplot()

  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = user_year, y = TAC, color = area_label,
                           group = interaction(area, iteration)),
                       alpha = 0.4, size = 0.5)
  }

  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75, fill = area_label),
                         alpha = 0.3)
  }

  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_TAC, color = area_label),
                       size = 1.5)
    p <- p + geom_point(data = summary_data,
                        aes(x = user_year, y = median_TAC, color = area_label),
                        size = 2)
  }

  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed",
               color = "red", alpha = 0.7) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(
      title = if(is.null(title)) "TAC by Area" else title,
      x = "Year",
      y = paste0("TAC (", units$TAC, ")"),
      color = "Area",
      fill = "Area"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # Facet if multiple areas
  if(length(areas) > 1) {
    p <- p + facet_wrap(~ area_label, scales = "free_y")
  }

  return(p)
}

# multifleet TAC helper (area × fleet)
create_TAC_plot_multifleet <- function(tac_data, areas, show_median, show_quantiles,
                                       show_individual, title, historical_end, units) {
  # expected cols in tac_data: year, area, fleet, TAC, (optional) iteration

  # prep
  plot_data <- tac_data %>%
    dplyr::mutate(
      user_year   = .data$year - 1,
      area_label  = paste("Area", .data$area),
      fleet_label = paste("Fleet", .data$fleet),
      period      = ifelse(.data$year <= historical_end, "Historical", "Projection")
    )

  has_iteration <- "iteration" %in% names(plot_data)

  # summaries per area × fleet × year
  summary_data <- plot_data %>%
    dplyr::group_by(.data$user_year, .data$area_label, .data$fleet_label, .data$period) %>%
    dplyr::summarise(
      median_TAC = stats::median(.data$TAC, na.rm = TRUE),
      q25        = stats::quantile(.data$TAC, 0.25, na.rm = TRUE),
      q75        = stats::quantile(.data$TAC, 0.75, na.rm = TRUE),
      .groups    = "drop"
    )

  # color palette per fleet present in filtered data
  fleets <- unique(summary_data$fleet_label)
  cols <- grDevices::rainbow(length(fleets))
  names(cols) <- fleets

  p <- ggplot2::ggplot()

  # individual iterations
  if (show_individual) {
    p <- p + ggplot2::geom_line(
      data = plot_data,
      ggplot2::aes(
        x = .data$user_year, y = .data$TAC, color = .data$fleet_label,
        group = if (has_iteration) interaction(.data$fleet_label, .data$iteration) else .data$fleet_label
      ),
      alpha = 0.35, linewidth = 0.4
    )
  }

  # quantile ribbons by fleet
  if (show_quantiles) {
    p <- p + ggplot2::geom_ribbon(
      data = summary_data,
      ggplot2::aes(x = .data$user_year, ymin = .data$q25, ymax = .data$q75, fill = .data$fleet_label),
      alpha = 0.25
    )
  }

  # medians
  if (show_median) {
    p <- p +
      ggplot2::geom_line(
        data = summary_data,
        ggplot2::aes(x = .data$user_year, y = .data$median_TAC, color = .data$fleet_label),
        linewidth = 1.2
      ) +
      ggplot2::geom_point(
        data = summary_data,
        ggplot2::aes(x = .data$user_year, y = .data$median_TAC, color = .data$fleet_label),
        size = 1.8
      )
  }

  y_lab <- if (!is.null(units$TAC)) paste0("TAC (", units$TAC, ")") else "TAC"

  p <- p +
    ggplot2::geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                        color = "red", alpha = 0.7) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(
      title = if (is.null(title)) "TAC by Area and Fleet" else title,
      x = "Year", y = y_lab, color = "Fleet", fill = "Fleet"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  # facet by area
  if (length(areas) > 1) {
    p <- p + ggplot2::facet_wrap(~ area_label, scales = "free_y")
  }

  p
}





#' Plot TAC by area (wrapper)
#' @param result Simulation result
#' @param areas Areas to plot
#' @param ... Additional arguments
#' @export
plot_TAC_by_area <- function(result, areas = "all", ...) {
  plot_TAC(result, areas = areas, show_fleets = FALSE, ...)
}

#' Plot TAC by fleet (wrapper)
#' @param result Simulation result
#' @param ... Additional arguments
#' @export
plot_TAC_by_fleet <- function(result, ...) {
  plot_TAC(result, show_fleets = TRUE, ...)
}


#' Plot total spawning biomass (sum across all areas)
#'
#' Creates time series plot of spawning biomass summed across all areas
#'
#' @param simulation_result Output object from \code{runProjection()}
#' @param show_median Logical. Whether to show median line. Default is TRUE
#' @param show_quantiles Logical. Whether to show quantile ribbon. Default is TRUE
#' @param show_individual Logical. Whether to show individual iteration lines. Default is FALSE
#' @param color_palette Character vector of colors. If NULL, uses default
#' @param title Character. Custom plot title. If NULL, generates automatic title
#' @param save_plot Logical. Whether to save plot. Default is FALSE
#' @param filename Character. Filename for saved plot
#' @param width Numeric. Plot width in inches. Default is 12
#' @param height Numeric. Plot height in inches. Default is 8
#'
#' @return A ggplot object
#' @export

plot_SB_total <- function(simulation_result,
                          show_median = TRUE,
                          show_quantiles = TRUE,
                          show_individual = FALSE,
                          color_palette = NULL,
                          title = NULL,
                          save_plot = FALSE,
                          filename = NULL,
                          width = 12,
                          height = 8) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1

  #get SB array: [year, iter, area]
  sb_array <- dynamics$SB

  #sum across areas (dimension 3)
  sb_total <- apply(sb_array, c(1, 2), sum, na.rm = TRUE)

  years <- dim(sb_total)[1]
  total_iterations <- dim(sb_total)[2]

  #prepare plot data
  plot_data <- data.frame()
  for(iter in 1:total_iterations) {
    iter_data <- data.frame(
      year = 1:years,
      user_year = 0:(years-1),
      value = sb_total[, iter],
      iteration = iter,
      period = ifelse(1:years <= historical_end, "Historical", "Projection")
    )
    plot_data <- rbind(plot_data, iter_data)
  }

  #calculate stats
  summary_data <- plot_data %>%
    group_by(user_year, period) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  #generate title
  if(is.null(title)) {
    title <- "Total Spawning Biomass"
  }

  if(is.null(color_palette)) {
    main_color <- "steelblue"
  } else {
    main_color <- color_palette[1]
  }

  p <- ggplot()

  #individual iterations?
  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = user_year, y = value, group = iteration),
                       color = main_color, alpha = 0.3, size = 0.5)
  }

  #quantile
  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75),
                         fill = main_color, alpha = 0.3)
  }

  #median line
  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value),
                       color = main_color, size = 1.5)
  }

  #formatting
  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed",
               color = "red", alpha = 0.7) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(title = title, x = "Year",
         y = get_metric_ylabel("SB", simulation_result)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  #save plot
  if(save_plot) {
    if(is.null(filename)) {
      filename <- "SB_total.jpeg"
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}

#' Plot total catch in numbers (sum across areas, by fleet if multifleet)
#'
#' @param simulation_result Output object from \code{runProjection()}
#' @param show_fleets Logical. Whether to show fleet breakdown in multifleet. Default is TRUE
#' @param show_median Logical. Whether to show median line. Default is TRUE
#' @param show_quantiles Logical. Whether to show quantile ribbon. Default is TRUE
#' @param show_individual Logical. Whether to show individual iteration lines. Default is FALSE
#' @param color_palette Character vector of colors. If NULL, uses default
#' @param title Character. Custom plot title. If NULL, generates automatic title
#' @param save_plot Logical. Whether to save plot. Default is FALSE
#' @param filename Character. Filename for saved plot
#' @param width Numeric. Plot width in inches. Default is 12
#' @param height Numeric. Plot height in inches. Default is 8
#'
#' @return A ggplot object
#' @export

plot_catchN_total <- function(simulation_result,
                              show_fleets = TRUE,
                              show_median = TRUE,
                              show_quantiles = TRUE,
                              show_individual = FALSE,
                              color_palette = NULL,
                              title = NULL,
                              save_plot = FALSE,
                              filename = NULL,
                              width = 12,
                              height = 8) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  is_multifleet <- !is.null(dynamics$multifleet)

  if(is_multifleet && show_fleets) {
    #using fleet-specific data: [year, iter, area, fleet]
    catchN_array <- dynamics$multifleet$catchN_by_fleet

    #sum across areas (dimension 3) for each fleet
    years <- dim(catchN_array)[1]
    total_iterations <- dim(catchN_array)[2]
    nfleets <- dim(catchN_array)[4]

    plot_data <- data.frame()

    for(fleet in 1:nfleets) {
      catchN_fleet <- apply(catchN_array[, , , fleet], c(1, 2), sum, na.rm = TRUE)

      for(iter in 1:total_iterations) {
        iter_data <- data.frame(
          year = 1:years,
          user_year = 0:(years-1),
          value = catchN_fleet[, iter],
          iteration = iter,
          fleet = paste("Fleet", fleet),
          period = ifelse(1:years <= historical_end, "Historical", "Projection")
        )
        plot_data <- rbind(plot_data, iter_data)
      }
    }

    #statistics by fleet
    summary_data <- plot_data %>%
      group_by(user_year, fleet, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )


    if(is.null(color_palette)) {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
    } else {
      colors <- color_palette[1:nfleets]
    }

    if(is.null(title)) {
      title <- "Total Catch in Numbers by Fleet"
    }

    p <- ggplot()

    #individual iterations?
    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, color = fleet,
                             group = interaction(fleet, iteration)),
                         alpha = 0.3, size = 0.5)
    }

    #quantile
    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75, fill = fleet),
                           alpha = 0.3)
    }

    #median lines
    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value, color = fleet),
                         size = 1.5)
    }

    #formatting
    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year",
           y = get_metric_ylabel("catchN", simulation_result),
           color = "Fleet", fill = "Fleet") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom"
      )

  } else {
    #single fleet or aggregated
    catchN_array <- dynamics$catchN

    #sum across areas (dimension 3)
    catchN_total <- apply(catchN_array, c(1, 2), sum, na.rm = TRUE)

    years <- dim(catchN_total)[1]
    total_iterations <- dim(catchN_total)[2]

    plot_data <- data.frame()
    for(iter in 1:total_iterations) {
      iter_data <- data.frame(
        year = 1:years,
        user_year = 0:(years-1),
        value = catchN_total[, iter],
        iteration = iter,
        period = ifelse(1:years <= historical_end, "Historical", "Projection")
      )
      plot_data <- rbind(plot_data, iter_data)
    }

    # statistics
    summary_data <- plot_data %>%
      group_by(user_year, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(title)) {
      title <- "Total Catch in Numbers"
    }

    if(is.null(color_palette)) {
      main_color <- "chocolate"
    } else {
      main_color <- color_palette[1]
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, group = iteration),
                         color = main_color, alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75),
                           fill = main_color, alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value),
                         color = main_color, size = 1.5)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year",
           y = get_metric_ylabel("catchN", simulation_result)) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
  }

  #save plot
  if(save_plot) {
    if(is.null(filename)) {
      fleet_suffix <- if(is_multifleet && show_fleets) "by_fleet" else "total"
      filename <- paste0("catchN_total_", fleet_suffix, ".jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


#' Plot total catch biomass (sum across areas, by fleet if multifleet)
#'
#' @param simulation_result Output object from \code{runProjection()}
#' @param show_fleets Logical. Whether to show fleet breakdown in multifleet. Default is TRUE
#' @param show_median Logical. Whether to show median line. Default is TRUE
#' @param show_quantiles Logical. Whether to show quantile ribbon. Default is TRUE
#' @param show_individual Logical. Whether to show individual iteration lines. Default is FALSE
#' @param color_palette Character vector of colors. If NULL, uses default
#' @param title Character. Custom plot title. If NULL, generates automatic title
#' @param save_plot Logical. Whether to save plot. Default is FALSE
#' @param filename Character. Filename for saved plot
#' @param width Numeric. Plot width in inches. Default is 12
#' @param height Numeric. Plot height in inches. Default is 8
#'
#' @return A ggplot object
#' @export

plot_catchB_total <- function(simulation_result,
                              show_fleets = TRUE,
                              show_median = TRUE,
                              show_quantiles = TRUE,
                              show_individual = FALSE,
                              color_palette = NULL,
                              title = NULL,
                              save_plot = FALSE,
                              filename = NULL,
                              width = 12,
                              height = 8) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  is_multifleet <- !is.null(dynamics$multifleet)

  if(is_multifleet && show_fleets) {
    #using fleet-specific data: [year, iter, area, fleet]
    catchB_array <- dynamics$multifleet$catchB_by_fleet

    #sum across areas (dimension 3) for each fleet
    years <- dim(catchB_array)[1]
    total_iterations <- dim(catchB_array)[2]
    nfleets <- dim(catchB_array)[4]

    plot_data <- data.frame()

    for(fleet in 1:nfleets) {
      catchB_fleet <- apply(catchB_array[, , , fleet], c(1, 2), sum, na.rm = TRUE)

      for(iter in 1:total_iterations) {
        iter_data <- data.frame(
          year = 1:years,
          user_year = 0:(years-1),
          value = catchB_fleet[, iter],
          iteration = iter,
          fleet = paste("Fleet", fleet),
          period = ifelse(1:years <= historical_end, "Historical", "Projection")
        )
        plot_data <- rbind(plot_data, iter_data)
      }
    }

    #statistics by fleet
    summary_data <- plot_data %>%
      group_by(user_year, fleet, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(color_palette)) {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
    } else {
      colors <- color_palette[1:nfleets]
    }

    if(is.null(title)) {
      title <- "Total Catch Biomass by Fleet"
    }

    p <- ggplot()

    #individual iterations?
    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, color = fleet,
                             group = interaction(fleet, iteration)),
                         alpha = 0.3, size = 0.5)
    }

    #quantile
    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75, fill = fleet),
                           alpha = 0.3)
    }

    #median lines
    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value, color = fleet),
                         size = 1.5)
    }

    #formatting
    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year",
           y = get_metric_ylabel("catchB", simulation_result),
           color = "Fleet", fill = "Fleet") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom"
      )

  } else {
    #single fleet or aggregated
    catchB_array <- dynamics$catchB

    #sum across areas (dimension 3)
    catchB_total <- apply(catchB_array, c(1, 2), sum, na.rm = TRUE)

    years <- dim(catchB_total)[1]
    total_iterations <- dim(catchB_total)[2]

    plot_data <- data.frame()
    for(iter in 1:total_iterations) {
      iter_data <- data.frame(
        year = 1:years,
        user_year = 0:(years-1),
        value = catchB_total[, iter],
        iteration = iter,
        period = ifelse(1:years <= historical_end, "Historical", "Projection")
      )
      plot_data <- rbind(plot_data, iter_data)
    }

    #statistics
    summary_data <- plot_data %>%
      group_by(user_year, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(title)) {
      title <- "Total Catch Biomass"
    }

    if(is.null(color_palette)) {
      main_color <- "brown"
    } else {
      main_color <- color_palette[1]
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, group = iteration),
                         color = main_color, alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75),
                           fill = main_color, alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value),
                         color = main_color, size = 1.5)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year",
           y = get_metric_ylabel("catchB", simulation_result)) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
  }

  #save plot
  if(save_plot) {
    if(is.null(filename)) {
      fleet_suffix <- if(is_multifleet && show_fleets) "by_fleet" else "total"
      filename <- paste0("catchB_total_", fleet_suffix, ".jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}

#' Plot total observed catch (sum across areas, by fleet if multifleet)
#'
#' @param simulation_result Output object from \code{runProjection()}
#' @param show_fleets Logical. Whether to show fleet breakdown in multifleet. Default is TRUE
#' @param show_median Logical. Whether to show median line. Default is TRUE
#' @param show_quantiles Logical. Whether to show quantile ribbon. Default is TRUE
#' @param show_individual Logical. Whether to show individual iteration points. Default is FALSE
#' @param point_size Numeric. Size of points. Default is 1.5
#' @param color_palette Character vector of colors. If NULL, uses default
#' @param title Character. Custom plot title. If NULL, generates automatic title
#' @param save_plot Logical. Whether to save plot. Default is FALSE
#' @param filename Character. Filename for saved plot
#' @param width Numeric. Plot width in inches. Default is 12
#' @param height Numeric. Plot height in inches. Default is 8
#'
#' @return A ggplot object
#' @export
plot_observed_catch_total <- function(simulation_result,
                                      show_fleets = TRUE,
                                      show_median = TRUE,
                                      show_quantiles = TRUE,
                                      show_individual = FALSE,
                                      point_size = 1.5,
                                      color_palette = NULL,
                                      title = NULL,
                                      save_plot = FALSE,
                                      filename = NULL,
                                      width = 12,
                                      height = 8) {

  obs_data <- simulation_result$HCR$decisionData
  if(is.null(obs_data)) {
    stop("no observation data found in simulation result")
  }

  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1

  #detect multifleet
  has_multifleet <- any(grepl("fleet_\\d+_observed_catch$", names(obs_data)))

  if(has_multifleet && show_fleets) {
    #find fleet observed catch columns
    obs_cols <- grep("fleet_\\d+_observed_catch$", names(obs_data), value = TRUE)
    fleet_nums <- unique(gsub("fleet_(\\d+)_observed_catch", "\\1", obs_cols))

    plot_data <- data.frame()

    for(fleet_num in fleet_nums) {
      obs_col <- paste0("fleet_", fleet_num, "_observed_catch")

      if(obs_col %in% names(obs_data)) {
        #extract data for this fleet
        fleet_data <- obs_data %>%
          filter(!is.na(!!sym(obs_col))) %>%
          mutate(
            user_year = j - 1,
            iteration = k,
            value = !!sym(obs_col),
            fleet = paste("Fleet", fleet_num),
            period = ifelse(j <= historical_end, "Historical", "Projection")
          ) %>%
          select(user_year, iteration, value, fleet, period)

        plot_data <- rbind(plot_data, fleet_data)
      }
    }

    if(nrow(plot_data) == 0) {
      stop("no multifleet observed catch data found")
    }

    #statistics by fleet
    summary_data <- plot_data %>%
      group_by(user_year, fleet, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
    nfleets <- length(fleet_nums)
    if(is.null(color_palette)) {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
    } else {
      colors <- color_palette[1:nfleets]
    }
    names(colors) <- paste("Fleet", fleet_nums)

    if(is.null(title)) {
      title <- "Total Observed Catch by Fleet"
    }

    p <- ggplot()

    #individual iterations?
    if(show_individual) {
      p <- p + geom_point(data = plot_data,
                          aes(x = user_year, y = value, color = fleet),
                          alpha = 0.4, size = point_size * 0.7)
    }

    #quantile ribbons
    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75, fill = fleet),
                           alpha = 0.3)
    }

    #median lines and points
    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value, color = fleet),
                         size = 1.5)
      p <- p + geom_point(data = summary_data,
                          aes(x = user_year, y = median_value, color = fleet),
                          size = point_size)
    }

    #formatting
    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Observed Catch",
           color = "Fleet", fill = "Fleet") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom"
      )

  } else {
    #single fleet or aggregated - use total observed catch
    if(!"observed_catch" %in% names(obs_data)) {
      stop("no observed_catch column found in observation data")
    }

    plot_data <- obs_data %>%
      filter(!is.na(observed_catch)) %>%
      mutate(
        user_year = j - 1,
        iteration = k,
        value = observed_catch,
        period = ifelse(j <= historical_end, "Historical", "Projection")
      ) %>%
      select(user_year, iteration, value, period)

    # Summary statistics
    summary_data <- plot_data %>%
      group_by(user_year, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(title)) {
      title <- "Total Observed Catch"
    }

    if(is.null(color_palette)) {
      main_color <- "orange"
    } else {
      main_color <- color_palette[1]
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_point(data = plot_data,
                          aes(x = user_year, y = value),
                          color = main_color, alpha = 0.4, size = point_size * 0.7)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75),
                           fill = main_color, alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value),
                         color = main_color, size = 1.5)
      p <- p + geom_point(data = summary_data,
                          aes(x = user_year, y = median_value),
                          color = main_color, size = point_size)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Observed Catch") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
      )
  }

  #save plot
  if(save_plot) {
    if(is.null(filename)) {
      fleet_suffix <- if(has_multifleet && show_fleets) "by_fleet" else "total"
      filename <- paste0("observed_catch_total_", fleet_suffix, ".jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


# ============================================================================
# MODIFIED PLOTTING FUNCTIONS WITH ENHANCED FONTS AND STANDARDIZATION
# Added: Oct 22, 2025
# These functions improve the original plotting functions above
# Changes: Larger fonts, standardization by last historical year
# ============================================================================

# Helper function for standardization
get_last_historical_value_mod <- function(data, historical_end, metric_col = "value") {
  last_hist_data <- data %>%
    filter(user_year == (historical_end - 1)) %>%
    pull(!!sym(metric_col))

  if(length(last_hist_data) == 0) return(1)
  return(median(last_hist_data, na.rm = TRUE))
}


#' Modified Total Catch Biomass by Fleet
#'
#' Enhanced version with larger fonts and standardization
#' @inheritParams plot_catchB_total
#' @export
plot_catchB_total_modified <- function(simulation_result,
                                       show_fleets = TRUE,
                                       show_median = TRUE,
                                       show_quantiles = TRUE,
                                       show_individual = FALSE,
                                       color_palette = NULL,
                                       title = NULL,
                                       save_plot = FALSE,
                                       filename = NULL,
                                       width = 14,
                                       height = 10) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  is_multifleet <- !is.null(dynamics$multifleet)

  if(is_multifleet && show_fleets) {
    catchB_array <- dynamics$multifleet$catchB_by_fleet
    years <- dim(catchB_array)[1]
    total_iterations <- dim(catchB_array)[2]
    nfleets <- dim(catchB_array)[4]

    plot_data <- data.frame()

    for(fleet in 1:nfleets) {
      catchB_fleet <- apply(catchB_array[, , , fleet], c(1, 2), sum, na.rm = TRUE)

      for(iter in 1:total_iterations) {
        iter_data <- data.frame(
          year = 1:years,
          user_year = 0:(years-1),
          value = catchB_fleet[, iter],
          iteration = iter,
          fleet = paste("Fleet", fleet),
          period = ifelse(1:years <= historical_end, "Historical", "Projection")
        )
        plot_data <- rbind(plot_data, iter_data)
      }
    }

    std_factor <- get_last_historical_value_mod(plot_data, historical_end)
    plot_data$value <- plot_data$value / std_factor

    summary_data <- plot_data %>%
      group_by(user_year, fleet, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(color_palette)) {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
    } else {
      colors <- color_palette[1:nfleets]
    }

    if(is.null(title)) {
      title <- "Total Catch Biomass by Fleet\n(Standardized by Last Historical Year)"
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, color = fleet,
                             group = interaction(fleet, iteration)),
                         alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75, fill = fleet),
                           alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value, color = fleet),
                         size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative Catch Biomass",
           color = "Fleet", fill = "Fleet") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20)),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(1.5, "cm")
      )

  } else {
    catchB_array <- dynamics$catchB
    catchB_total <- apply(catchB_array, c(1, 2), sum, na.rm = TRUE)
    years <- dim(catchB_total)[1]
    total_iterations <- dim(catchB_total)[2]

    plot_data <- data.frame()
    for(iter in 1:total_iterations) {
      iter_data <- data.frame(
        year = 1:years,
        user_year = 0:(years-1),
        value = catchB_total[, iter],
        iteration = iter,
        period = ifelse(1:years <= historical_end, "Historical", "Projection")
      )
      plot_data <- rbind(plot_data, iter_data)
    }

    std_factor <- get_last_historical_value_mod(plot_data, historical_end)
    plot_data$value <- plot_data$value / std_factor

    summary_data <- plot_data %>%
      group_by(user_year, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(title)) {
      title <- "Total Catch Biomass\n(Standardized by Last Historical Year)"
    }

    if(is.null(color_palette)) {
      main_color <- "brown"
    } else {
      main_color <- color_palette[1]
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, group = iteration),
                         color = main_color, alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75),
                           fill = main_color, alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value),
                         color = main_color, size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative Catch Biomass") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20))
      )
  }

  if(save_plot) {
    if(is.null(filename)) {
      fleet_suffix <- if(is_multifleet && show_fleets) "by_fleet" else "total"
      filename <- paste0("catchB_total_", fleet_suffix, "_modified.jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


#' Modified Total Catch Numbers by Fleet
#'
#' Enhanced version with larger fonts and standardization
#' @inheritParams plot_catchN_total
#' @export
plot_catchN_total_modified <- function(simulation_result,
                                       show_fleets = TRUE,
                                       show_median = TRUE,
                                       show_quantiles = TRUE,
                                       show_individual = FALSE,
                                       color_palette = NULL,
                                       title = NULL,
                                       save_plot = FALSE,
                                       filename = NULL,
                                       width = 14,
                                       height = 10) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  is_multifleet <- !is.null(dynamics$multifleet)

  if(is_multifleet && show_fleets) {
    catchN_array <- dynamics$multifleet$catchN_by_fleet
    years <- dim(catchN_array)[1]
    total_iterations <- dim(catchN_array)[2]
    nfleets <- dim(catchN_array)[4]

    plot_data <- data.frame()

    for(fleet in 1:nfleets) {
      catchN_fleet <- apply(catchN_array[, , , fleet], c(1, 2), sum, na.rm = TRUE)

      for(iter in 1:total_iterations) {
        iter_data <- data.frame(
          year = 1:years,
          user_year = 0:(years-1),
          value = catchN_fleet[, iter],
          iteration = iter,
          fleet = paste("Fleet", fleet),
          period = ifelse(1:years <= historical_end, "Historical", "Projection")
        )
        plot_data <- rbind(plot_data, iter_data)
      }
    }

    std_factor <- get_last_historical_value_mod(plot_data, historical_end)
    plot_data$value <- plot_data$value / std_factor

    summary_data <- plot_data %>%
      group_by(user_year, fleet, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(color_palette)) {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
    } else {
      colors <- color_palette[1:nfleets]
    }

    if(is.null(title)) {
      title <- "Total Catch Numbers by Fleet\n(Standardized by Last Historical Year)"
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, color = fleet,
                             group = interaction(fleet, iteration)),
                         alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75, fill = fleet),
                           alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value, color = fleet),
                         size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative Catch Numbers",
           color = "Fleet", fill = "Fleet") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20)),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(1.5, "cm")
      )

  } else {
    catchN_array <- dynamics$catchN
    catchN_total <- apply(catchN_array, c(1, 2), sum, na.rm = TRUE)
    years <- dim(catchN_total)[1]
    total_iterations <- dim(catchN_total)[2]

    plot_data <- data.frame()
    for(iter in 1:total_iterations) {
      iter_data <- data.frame(
        year = 1:years,
        user_year = 0:(years-1),
        value = catchN_total[, iter],
        iteration = iter,
        period = ifelse(1:years <= historical_end, "Historical", "Projection")
      )
      plot_data <- rbind(plot_data, iter_data)
    }

    std_factor <- get_last_historical_value_mod(plot_data, historical_end)
    plot_data$value <- plot_data$value / std_factor

    summary_data <- plot_data %>%
      group_by(user_year, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(title)) {
      title <- "Total Catch Numbers\n(Standardized by Last Historical Year)"
    }

    if(is.null(color_palette)) {
      main_color <- "brown"
    } else {
      main_color <- color_palette[1]
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, group = iteration),
                         color = main_color, alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75),
                           fill = main_color, alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value),
                         color = main_color, size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative Catch Numbers") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20))
      )
  }

  if(save_plot) {
    if(is.null(filename)) {
      fleet_suffix <- if(is_multifleet && show_fleets) "by_fleet" else "total"
      filename <- paste0("catchN_total_", fleet_suffix, "_modified.jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}

#' Modified Total Spawning Biomass
#'
#' Enhanced version with larger fonts and standardization
#' @inheritParams plot_SB_total
#' @export
plot_SB_total_modified <- function(simulation_result,
                                   show_median = TRUE,
                                   show_quantiles = TRUE,
                                   show_individual = FALSE,
                                   color_palette = NULL,
                                   title = NULL,
                                   save_plot = FALSE,
                                   filename = NULL,
                                   width = 14,
                                   height = 10) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1

  SB_array <- dynamics$SB
  SB_total <- apply(SB_array, c(1, 2), sum, na.rm = TRUE)

  years <- dim(SB_total)[1]
  total_iterations <- dim(SB_total)[2]

  plot_data <- data.frame()
  for(iter in 1:total_iterations) {
    iter_data <- data.frame(
      year = 1:years,
      user_year = 0:(years-1),
      value = SB_total[, iter],
      iteration = iter,
      period = ifelse(1:years <= historical_end, "Historical", "Projection")
    )
    plot_data <- rbind(plot_data, iter_data)
  }

  std_factor <- get_last_historical_value_mod(plot_data, historical_end)
  plot_data$value <- plot_data$value / std_factor

  summary_data <- plot_data %>%
    group_by(user_year, period) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  if(is.null(title)) {
    title <- "Total Spawning Biomass\n(Standardized by Last Historical Year)"
  }

  if(is.null(color_palette)) {
    main_color <- "steelblue"
  } else {
    main_color <- color_palette[1]
  }

  p <- ggplot()

  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = user_year, y = value, group = iteration),
                       color = main_color, alpha = 0.3, size = 0.5)
  }

  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75),
                         fill = main_color, alpha = 0.3)
  }

  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value),
                       color = main_color, size = 2)
  }

  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed",
               color = "red", alpha = 0.7, size = 1) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(title = title, x = "Year", y = "Relative Spawning Biomass") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20))
    )

  if(save_plot) {
    if(is.null(filename)) {
      filename <- "SB_total_modified.jpeg"
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


#' Modified Recruitment
#'
#' Enhanced version with larger fonts and standardization
#' @inheritParams plot_recN
#' @export
plot_recN_modified <- function(simulation_result,
                               show_median = TRUE,
                               show_quantiles = TRUE,
                               show_individual = FALSE,
                               color_palette = NULL,
                               title = NULL,
                               save_plot = FALSE,
                               filename = NULL,
                               width = 14,
                               height = 10) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1

  recN_array <- dynamics$recN
  years <- dim(recN_array)[1]
  total_iterations <- dim(recN_array)[2]

  plot_data <- data.frame()
  for(iter in 1:total_iterations) {
    iter_data <- data.frame(
      year = 1:years,
      user_year = 0:(years-1),
      value = recN_array[, iter],
      iteration = iter,
      period = ifelse(1:years <= historical_end, "Historical", "Projection")
    )
    plot_data <- rbind(plot_data, iter_data)
  }

  std_factor <- get_last_historical_value_mod(plot_data, historical_end)
  plot_data$value <- plot_data$value / std_factor

  summary_data <- plot_data %>%
    group_by(user_year, period) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  if(is.null(title)) {
    title <- "Total Recruitment\n(Standardized by Last Historical Year)"
  }

  if(is.null(color_palette)) {
    main_color <- "darkgreen"
  } else {
    main_color <- color_palette[1]
  }

  p <- ggplot()

  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = user_year, y = value, group = iteration),
                       color = main_color, alpha = 0.3, size = 0.5)
  }

  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75),
                         fill = main_color, alpha = 0.3)
  }

  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value),
                       color = main_color, size = 2)
  }

  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed",
               color = "red", alpha = 0.7, size = 1) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(title = title, x = "Year", y = "Relative Recruitment") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20))
    )

  if(save_plot) {
    if(is.null(filename)) {
      filename <- "recN_modified.jpeg"
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


#' Modified Total Discards by Fleet (numbers)
#'
#' Enhanced version with larger fonts and standardization
#' @inheritParams plot_catchN_total
#' @export
plot_discN_total_modified <- function(simulation_result,
                                      show_fleets = TRUE,
                                      show_median = TRUE,
                                      show_quantiles = TRUE,
                                      show_individual = FALSE,
                                      color_palette = NULL,
                                      title = NULL,
                                      save_plot = FALSE,
                                      filename = NULL,
                                      width = 14,
                                      height = 10) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  is_multifleet <- !is.null(dynamics$multifleet)

  if(is_multifleet && show_fleets) {
    discN_array <- dynamics$multifleet$discN_by_fleet
    years <- dim(discN_array)[1]
    total_iterations <- dim(discN_array)[2]
    nfleets <- dim(discN_array)[4]

    plot_data <- data.frame()

    for(fleet in 1:nfleets) {
      discN_fleet <- apply(discN_array[, , , fleet], c(1, 2), sum, na.rm = TRUE)

      for(iter in 1:total_iterations) {
        iter_data <- data.frame(
          year = 1:years,
          user_year = 0:(years-1),
          value = discN_fleet[, iter],
          iteration = iter,
          fleet = paste("Fleet", fleet),
          period = ifelse(1:years <= historical_end, "Historical", "Projection")
        )
        plot_data <- rbind(plot_data, iter_data)
      }
    }

    std_factor <- get_last_historical_value_mod(plot_data, historical_end)
    plot_data$value <- plot_data$value / std_factor

    summary_data <- plot_data %>%
      group_by(user_year, fleet, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(color_palette)) {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
    } else {
      colors <- color_palette[1:nfleets]
    }

    if(is.null(title)) {
      title <- "Total Discards by Fleet (Numbers)\n(Standardized by Last Historical Year)"
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, color = fleet,
                             group = interaction(fleet, iteration)),
                         alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75, fill = fleet),
                           alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value, color = fleet),
                         size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative Discard Numbers",
           color = "Fleet", fill = "Fleet") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20)),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(1.5, "cm")
      )

  } else {
    discN_array <- dynamics$discN
    discN_total <- apply(discN_array, c(1, 2), sum, na.rm = TRUE)
    years <- dim(discN_total)[1]
    total_iterations <- dim(discN_total)[2]

    plot_data <- data.frame()
    for(iter in 1:total_iterations) {
      iter_data <- data.frame(
        year = 1:years,
        user_year = 0:(years-1),
        value = discN_total[, iter],
        iteration = iter,
        period = ifelse(1:years <= historical_end, "Historical", "Projection")
      )
      plot_data <- rbind(plot_data, iter_data)
    }

    std_factor <- get_last_historical_value_mod(plot_data, historical_end)
    plot_data$value <- plot_data$value / std_factor

    summary_data <- plot_data %>%
      group_by(user_year, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(title)) {
      title <- "Total Discards (Numbers)\n(Standardized by Last Historical Year)"
    }

    if(is.null(color_palette)) {
      main_color <- "brown"
    } else {
      main_color <- color_palette[1]
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, group = iteration),
                         color = main_color, alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75),
                           fill = main_color, alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value),
                         color = main_color, size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative Discard Numbers") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20))
      )
  }

  if(save_plot) {
    if(is.null(filename)) {
      fleet_suffix <- if(is_multifleet && show_fleets) "by_fleet" else "total"
      filename <- paste0("discN_total_", fleet_suffix, "_modified.jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


#' Modified Total Fishing Mortality by Fleet (Area 1 only)
#'
#' Enhanced version with larger fonts and standardization
#' Uses only Area 1 as total F
#' @inheritParams plot_Ftotal_multi
#' @export
plot_Ftotal_modified <- function(simulation_result,
                                 show_fleets = TRUE,
                                 show_median = TRUE,
                                 show_quantiles = TRUE,
                                 show_individual = FALSE,
                                 color_palette = NULL,
                                 title = NULL,
                                 save_plot = FALSE,
                                 filename = NULL,
                                 width = 14,
                                 height = 10) {

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  is_multifleet <- !is.null(dynamics$multifleet)

  if(is_multifleet && show_fleets) {
    Ftotal_array <- dynamics$multifleet$Ftotal_by_fleet
    years <- dim(Ftotal_array)[1]
    total_iterations <- dim(Ftotal_array)[2]
    nfleets <- dim(Ftotal_array)[4]

    plot_data <- data.frame()

    # Extract Area 1 only (dimension 3, index 1)
    for(fleet in 1:nfleets) {
      Ftotal_fleet <- Ftotal_array[, , 1, fleet]  # Area 1 only

      for(iter in 1:total_iterations) {
        iter_data <- data.frame(
          year = 1:years,
          user_year = 0:(years-1),
          value = Ftotal_fleet[, iter],
          iteration = iter,
          fleet = paste("Fleet", fleet),
          period = ifelse(1:years <= historical_end, "Historical", "Projection")
        )
        plot_data <- rbind(plot_data, iter_data)
      }
    }

    std_factor <- get_last_historical_value_mod(plot_data, historical_end)
    plot_data$value <- plot_data$value / std_factor

    summary_data <- plot_data %>%
      group_by(user_year, fleet, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(color_palette)) {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
    } else {
      colors <- color_palette[1:nfleets]
    }

    if(is.null(title)) {
      title <- "Total Fishing Mortality by Fleet (Area 1)\n(Standardized by Last Historical Year)"
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, color = fleet,
                             group = interaction(fleet, iteration)),
                         alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75, fill = fleet),
                           alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value, color = fleet),
                         size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative Fishing Mortality",
           color = "Fleet", fill = "Fleet") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20)),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(1.5, "cm")
      )

  } else {
    Ftotal_array <- dynamics$Ftotal
    # Extract Area 1 only
    Ftotal_area1 <- Ftotal_array[, , 1]

    years <- dim(Ftotal_area1)[1]
    total_iterations <- dim(Ftotal_area1)[2]

    plot_data <- data.frame()
    for(iter in 1:total_iterations) {
      iter_data <- data.frame(
        year = 1:years,
        user_year = 0:(years-1),
        value = Ftotal_area1[, iter],
        iteration = iter,
        period = ifelse(1:years <= historical_end, "Historical", "Projection")
      )
      plot_data <- rbind(plot_data, iter_data)
    }

    std_factor <- get_last_historical_value_mod(plot_data, historical_end)
    plot_data$value <- plot_data$value / std_factor

    summary_data <- plot_data %>%
      group_by(user_year, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(title)) {
      title <- "Total Fishing Mortality (Area 1)\n(Standardized by Last Historical Year)"
    }

    if(is.null(color_palette)) {
      main_color <- "brown"
    } else {
      main_color <- color_palette[1]
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, group = iteration),
                         color = main_color, alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75),
                           fill = main_color, alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value),
                         color = main_color, size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative Fishing Mortality") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20))
      )
  }

  if(save_plot) {
    if(is.null(filename)) {
      fleet_suffix <- if(is_multifleet && show_fleets) "by_fleet" else "total"
      filename <- paste0("Ftotal_", fleet_suffix, "_modified.jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


#' Modified Indices Plots (Survey and CPUE)
#'
#' Enhanced version with larger fonts, keeps original flexible pattern matching
#' @param simulation_result Output from runProjection
#' @param index_pattern Pattern to match (e.g., "IDX_Survey_1", "IDX_CPUE_2")
#' @param show_median Show median line
#' @param show_quantiles Show quantile ribbon
#' @param show_individual Show individual iterations
#' @param point_size Size of points
#' @param line_alpha Transparency
#' @param color_palette Custom colors
#' @param title Custom title
#' @param save_plot Save to file
#' @param filename Custom filename
#' @param width Plot width
#' @param height Plot height
#' @export
plot_indices_modified <- function(simulation_result,
                                  index_pattern = "IDX_",
                                  show_median = TRUE,
                                  show_quantiles = TRUE,
                                  show_individual = FALSE,
                                  point_size = 2,
                                  line_alpha = 0.5,
                                  color_palette = NULL,
                                  title = NULL,
                                  save_plot = FALSE,
                                  filename = NULL,
                                  width = 14,
                                  height = 10) {

  obs_data <- simulation_result$HCR$decisionData
  if(is.null(obs_data)) {
    stop("no observation data found in simulation result")
  }

  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1

  #find index columns matching the pattern
  index_cols <- grep(paste0("^", index_pattern), names(obs_data), value = TRUE)
  index_cols <- index_cols[!grepl("_(areas|indexYears|fleet_id|indextype|selectivity_|survey_timing)", index_cols)]

  if(length(index_cols) == 0) {
    stop("No index columns found matching pattern: ", index_pattern)
  }

  #prepare data
  plot_data <- data.frame()

  for(idx_col in index_cols) {
    valid_rows <- !is.na(obs_data[[idx_col]])
    if(any(valid_rows)) {
      temp_df <- obs_data[valid_rows, ] %>%
        mutate(
          user_year = j - 1,
          value = !!sym(idx_col),
          index_name = idx_col,
          index_clean = gsub("^IDX_(Survey_|CPUE_)", "", idx_col),
          index_clean = gsub("_", " ", index_clean),
          iteration = k,
          period = ifelse(j <= historical_end, "Historical", "Projection")
        ) %>%
        select(user_year, value, index_name, index_clean, iteration, period)

      plot_data <- rbind(plot_data, temp_df)
    }
  }

  if(nrow(plot_data) == 0) {
    stop("no valid data for plotting")
  }

  summary_data <- plot_data %>%
    group_by(user_year, index_name, index_clean, period) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE),
      q75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  if(is.null(title)) {
    title <- paste("Index:", index_pattern)
  }

  if(is.null(color_palette)) {
    main_color <- "steelblue"
  } else {
    main_color <- color_palette[1]
  }

  p <- ggplot()

  if(show_individual) {
    p <- p + geom_point(data = plot_data,
                        aes(x = user_year, y = value),
                        color = main_color,
                        alpha = line_alpha, size = point_size * 0.7)
  }

  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = user_year, ymin = q25, ymax = q75),
                         fill = main_color, alpha = 0.3)
  }

  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = user_year, y = median_value),
                       color = main_color, size = 2)
    p <- p + geom_point(data = summary_data,
                        aes(x = user_year, y = median_value),
                        color = main_color, size = point_size)
  }

  p <- p +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed",
               color = "red", alpha = 0.7, size = 1) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(title = title, x = "Year", y = "Index Value") +
    theme_minimal() +
    theme(
      # INCREASED FONT SIZES
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20))
    )

  if(save_plot) {
    if(is.null(filename)) {
      pattern_clean <- gsub("_$", "", gsub("^IDX_", "", index_pattern))
      filename <- paste0("indices_", pattern_clean, "_modified.jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


#' Modified Length Composition
#'
#' Enhanced version with larger fonts and proportions, uses program_pattern like original
#' @param result Simulation result
#' @param program_pattern Pattern to match (e.g., "LC_Fishery", "LC_Survey")
#' @param area_filter Areas to include
#' @param fleet_filter Fleets to include
#' @param years Years to plot or "all"
#' @param show_median Show median line
#' @param show_quantiles Show quantile ribbon
#' @param show_individual Show individual iterations
#' @param color_palette Custom colors
#' @param title Custom title
#' @param save_plot Save to file
#' @param filename Custom filename
#' @param width Plot width
#' @param height Plot height
#' @export
plot_length_comp_modified <- function(result,
                                      program_pattern = "LC_Fishery",
                                      area_filter = "all",
                                      fleet_filter = "all",
                                      years = "all",
                                      show_median = TRUE,
                                      show_quantiles = TRUE,
                                      show_individual = FALSE,
                                      color_palette = NULL,
                                      title = NULL,
                                      save_plot = FALSE,
                                      filename = NULL,
                                      width = 14,
                                      height = 10) {

  #ensure matrix format
  if(!is.null(result$HCR$decisionData) && is.null(result$observation_matrices)) {
    obs_data <- result$HCR$decisionData
    has_lc <- any(grepl("^LC_.*_count_bin_\\d+$", names(obs_data)))

    if(has_lc && exists("adapt_observation_data")) {
      result <- adapt_observation_data(result, verbose = FALSE)
    }
  }

  if(is.null(result$observation_matrices)) {
    stop("No observation matrices found. Ensure observation data is available.")
  }

  obs_matrices <- result$observation_matrices
  historical_end <- result$TimeAreaObj@historicalYears + 1

  #find length comp columns matching the pattern
  lc_pattern <- paste0("^", program_pattern, ".*_count_bin_\\d+$")
  lc_cols <- grep(lc_pattern, names(obs_matrices), value = TRUE)

  if(length(lc_cols) == 0) {
    stop(paste("No length composition data found for pattern:", program_pattern))
  }

  #apply area filter if specified
  if(!identical(area_filter, "all")) {
    area_pattern <- paste0("_area", area_filter, collapse = "|")
    lc_cols <- grep(area_pattern, lc_cols, value = TRUE)
  }

  #apply fleet filter if specified
  if(!identical(fleet_filter, "all")) {
    fleet_pattern <- paste0("_fleet", fleet_filter, collapse = "|")
    lc_cols <- grep(fleet_pattern, lc_cols, value = TRUE)
  }

  if(length(lc_cols) == 0) {
    stop("No length composition data after filtering")
  }

  # get unique program/area/fleet combinations
  programs <- unique(gsub("_count_bin_\\d+$", "", lc_cols))

  # extract bin numbers from first program
  first_program_cols <- grep(paste0("^", programs[1], "_count_bin"), lc_cols, value = TRUE)
  bin_nums <- as.numeric(gsub(".*_bin_(\\d+)$", "\\1", first_program_cols))
  bin_nums <- sort(unique(bin_nums))

  #get length bins
  lh <- result$LifeHistoryObj
  length_bins <- seq(0, lh@Linf * 1.5, length.out = max(bin_nums) + 1)
  bin_mids <- (length_bins[-1] + length_bins[-length(length_bins)]) / 2

  #prepare data - CONVERT TO PROPORTIONS
  plot_data <- data.frame()

  for(program in programs) {
    program_cols <- grep(paste0("^", program, "_count_bin"), lc_cols, value = TRUE)
    program_cols <- program_cols[order(as.numeric(gsub(".*_bin_(\\d+)$", "\\1", program_cols)))]

    first_matrix <- obs_matrices[[program_cols[1]]]
    all_years_available <- 0:(nrow(first_matrix) - 1)

    years_to_plot <- if(identical(years, "all")) all_years_available else intersect(all_years_available, years)

    for(year in years_to_plot) {
      year_idx <- year + 1

      for(iter in 1:ncol(first_matrix)) {
        #extract counts for all bins
        counts <- sapply(program_cols, function(col) {
          obs_matrices[[col]][year_idx, iter]
        })

        #convert to proportions
        total_count <- sum(counts, na.rm = TRUE)
        if(total_count > 0 && !all(is.na(counts))) {
          proportions <- counts / total_count

          temp_df <- data.frame(
            user_year = year,
            length_bin = bin_mids[bin_nums],
            proportion = proportions,
            iteration = iter,
            program = program,
            period = ifelse(year < (historical_end - 1), "Historical", "Projection")
          )
          plot_data <- rbind(plot_data, temp_df)
        }
      }
    }
  }

  if(nrow(plot_data) == 0) {
    stop("No valid length composition data to plot")
  }

  # statistics
  summary_data <- plot_data %>%
    group_by(user_year, length_bin, program, period) %>%
    summarise(
      median_prop = median(proportion, na.rm = TRUE),
      q25 = quantile(proportion, 0.25, na.rm = TRUE),
      q75 = quantile(proportion, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  if(is.null(title)) {
    length_units <- if(length(lh@L_units) > 0) lh@L_units else "cm"
    title <- paste0("Length Composition - ", program_pattern, " (Proportions)")
  }

  if(is.null(color_palette)) {
    main_color <- "steelblue"
  } else {
    main_color <- color_palette[1]
  }

  # plot
  p <- ggplot()

  if(show_individual) {
    p <- p + geom_line(data = plot_data,
                       aes(x = length_bin, y = proportion,
                           group = interaction(user_year, program, iteration)),
                       color = main_color, alpha = 0.2, size = 0.3)
  }

  if(show_quantiles) {
    p <- p + geom_ribbon(data = summary_data,
                         aes(x = length_bin, ymin = q25, ymax = q75,
                             group = interaction(user_year, program)),
                         fill = main_color, alpha = 0.3)
  }

  if(show_median) {
    p <- p + geom_line(data = summary_data,
                       aes(x = length_bin, y = median_prop,
                           group = interaction(user_year, program)),
                       color = main_color, size = 1.5)
  }

  p <- p +
    facet_wrap(~ user_year + program, scales = "free_y") +
    labs(title = title,
         x = paste0("Length (", if(length(lh@L_units) > 0) lh@L_units else "cm", ")"),
         y = "Proportion") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
      axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", margin = margin(b = 20)),
      strip.text = element_text(size = 14, face = "bold")
    )

  if(save_plot) {
    if(is.null(filename)) {
      pattern_clean <- gsub("^LC_", "", program_pattern)
      filename <- paste0("length_comp_", pattern_clean, "_modified.jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


#' Wrapper for fishery length compositions
#' @export
plot_fishery_length_comp_modified <- function(result, areas = "all", fleets = "all", ...) {
  plot_length_comp_modified(result, program_pattern = "LC_Fishery",
                            area_filter = areas, fleet_filter = fleets, ...)
}

#' Wrapper for survey length compositions
#' @export
plot_survey_length_comp_modified <- function(result, areas = "all", ...) {
  plot_length_comp_modified(result, program_pattern = "LC_Survey",
                            area_filter = areas, ...)
}

#' Modified Total TAC by Fleet
#'
#' Enhanced version with larger fonts and standardization by catch
#' Uses HCR$decisionAnnual like the original
#' @param simulation_result Output from runProjection
#' @param show_fleets Show fleet breakdown
#' @param show_median Show median line
#' @param show_quantiles Show quantile ribbon
#' @param show_individual Show individual iterations
#' @param color_palette Custom colors
#' @param title Custom title
#' @param save_plot Save to file
#' @param filename Custom filename
#' @param width Plot width
#' @param height Plot height
#' @export
plot_TAC_total_modified <- function(simulation_result,
                                    show_fleets = TRUE,
                                    show_median = TRUE,
                                    show_quantiles = TRUE,
                                    show_individual = FALSE,
                                    color_palette = NULL,
                                    title = NULL,
                                    save_plot = FALSE,
                                    filename = NULL,
                                    width = 14,
                                    height = 10) {

  #get TAC data from HCR
  tac_data <- simulation_result$HCR$decisionAnnual

  if(is.null(tac_data) || !"TAC" %in% names(tac_data)) {
    stop("no TAC data found in simulation result")
  }

  dynamics <- simulation_result$dynamics
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1
  is_multifleet <- !is.null(dynamics$multifleet)
  has_fleets <- "fleet" %in% names(tac_data)

  #get catch data for standardization
  if(is_multifleet && has_fleets && show_fleets) {
    catchB_array <- dynamics$multifleet$catchB_by_fleet
    years <- dim(catchB_array)[1]
    total_iterations <- dim(catchB_array)[2]
    nfleets <- dim(catchB_array)[4]

    #calculate standardization factors from catch (last historical year)
    catch_std_factors <- numeric(nfleets)

    for(fleet in 1:nfleets) {
      catchB_fleet <- apply(catchB_array[, , , fleet], c(1, 2), sum, na.rm = TRUE)
      last_hist_catch <- catchB_fleet[historical_end, ]
      catch_std_factors[fleet] <- median(last_hist_catch, na.rm = TRUE)
    }

    #prepare TAC plot data (sum across areas for each fleet)
    plot_data <- tac_data %>%
      group_by(year, fleet, iteration) %>%
      summarise(TAC = sum(TAC, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        user_year = year - 1,
        fleet_label = paste("Fleet", fleet),
        period = ifelse(year <= historical_end, "Historical", "Projection")
      ) %>%
      left_join(data.frame(fleet = 1:nfleets, std_factor = catch_std_factors), by = "fleet") %>%
      mutate(value = TAC / std_factor) %>%
      filter(period == "Projection")

    if(nrow(plot_data) == 0) {
      stop("No TAC data found in projection period")
    }

    summary_data <- plot_data %>%
      group_by(user_year, fleet_label, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(color_palette)) {
      colors <- c("steelblue", "darkgreen", "orange", "purple", "brown", "pink")[1:nfleets]
    } else {
      colors <- color_palette[1:nfleets]
    }
    names(colors) <- paste("Fleet", 1:nfleets)

    if(is.null(title)) {
      title <- "Total TAC by Fleet\n(Standardized by Last Historical Year Catch)"
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, color = fleet_label,
                             group = interaction(fleet_label, iteration)),
                         alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75, fill = fleet_label),
                           alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value, color = fleet_label),
                         size = 2)
      p <- p + geom_point(data = summary_data,
                          aes(x = user_year, y = median_value, color = fleet_label),
                          size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative TAC",
           color = "Fleet", fill = "Fleet") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20)),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(1.5, "cm")
      )

  } else {
    #single fleet or aggregated
    catchB_array <- dynamics$catchB
    catchB_total <- apply(catchB_array, c(1, 2), sum, na.rm = TRUE)

    #calculate standardization factor
    last_hist_catch <- catchB_total[historical_end, ]
    catch_std_factor <- median(last_hist_catch, na.rm = TRUE)

    #prepare TAC plot data (sum across areas)
    plot_data <- tac_data %>%
      group_by(year, iteration) %>%
      summarise(TAC = sum(TAC, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        user_year = year - 1,
        value = TAC / catch_std_factor,
        period = ifelse(year <= historical_end, "Historical", "Projection")
      ) %>%
      filter(period == "Projection")

    if(nrow(plot_data) == 0) {
      stop("No TAC data found in projection period")
    }

    summary_data <- plot_data %>%
      group_by(user_year, period) %>%
      summarise(
        median_value = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    if(is.null(title)) {
      title <- "Total TAC\n(Standardized by Last Historical Year Catch)"
    }

    if(is.null(color_palette)) {
      main_color <- "brown"
    } else {
      main_color <- color_palette[1]
    }

    p <- ggplot()

    if(show_individual) {
      p <- p + geom_line(data = plot_data,
                         aes(x = user_year, y = value, group = iteration),
                         color = main_color, alpha = 0.3, size = 0.5)
    }

    if(show_quantiles) {
      p <- p + geom_ribbon(data = summary_data,
                           aes(x = user_year, ymin = q25, ymax = q75),
                           fill = main_color, alpha = 0.3)
    }

    if(show_median) {
      p <- p + geom_line(data = summary_data,
                         aes(x = user_year, y = median_value),
                         color = main_color, size = 2)
      p <- p + geom_point(data = summary_data,
                          aes(x = user_year, y = median_value),
                          color = main_color, size = 2)
    }

    p <- p +
      geom_vline(xintercept = historical_end - 1, linetype = "dashed",
                 color = "red", alpha = 0.7, size = 1) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(title = title, x = "Year", y = "Relative TAC") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold", margin = margin(b = 20))
      )
  }

  if(save_plot) {
    if(is.null(filename)) {
      fleet_suffix <- if(has_fleets && show_fleets) "by_fleet" else "total"
      filename <- paste0("TAC_total_", fleet_suffix, "_modified.jpeg")
    }
    ggsave(filename, p, width = width, height = height, dpi = 300)
    cat("Plot saved as:", filename, "\n")
  }

  return(p)
}


# ============================================================================
# END OF MODIFIED PLOTTING FUNCTIONS
# ============================================================================

cat("\n")
cat("=========================================================================\n")
cat("Modified plotting functions loaded successfully!\n")
cat("=========================================================================\n")
cat("\nAvailable modified functions:\n")
cat("  1. plot_catchB_total_modified()      - Catch biomass by fleet\n")
cat("  2. plot_catchN_total_modified()      - Catch numbers by fleet\n")
cat("  3. plot_SB_total_modified()          - Spawning biomass\n")
cat("  4. plot_recN_modified()              - Recruitment\n")
cat("  5. plot_discN_total_modified()       - Discards by fleet\n")
cat("  6. plot_Ftotal_modified()            - Fishing mortality (Area 1)\n")
cat("  7. plot_indices_modified()           - Survey and CPUE indices\n")
cat("  8. plot_length_comp_modified()       - Length compositions (proportions)\n")
cat("  9. plot_TAC_total_modified()         - TAC by fleet (standardized by catch)\n")
cat("\nAll modified functions include:\n")
cat("  - Larger fonts (axis: 16pt, titles: 20pt, plot title: 22pt)\n")
cat("  - Standardization by last historical year (where applicable)\n")
cat("  - Enhanced legend sizes\n")
cat("=========================================================================\n")
cat("\n")



