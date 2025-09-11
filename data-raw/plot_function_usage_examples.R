#Plot usage examples


# # ============================================================================
# # TEST SCENARIO: SINGLE FLEET WITH NUMBER-BASED FD INDICES (BUG TEST)
# # ============================================================================
#
# cat("============================================\n")
# cat("RUNNING BUG TEST: SINGLE FLEET NUMBERS FD\n")
# cat("============================================\n")
#
# # Create number-based FD indices that will trigger manual calculation (to see if bug was fixed)
# cpue_numbers_single <- new("Index")
# cpue_numbers_single@indexID <- "CPUE_Numbers_BugTest"
# cpue_numbers_single@title <- "Single Fleet Numbers CPUE Test"
# cpue_numbers_single@useWeight <- FALSE  # forces manual calculation (numbers-based)
#
# cpue_numbers_single@survey_design <- list(
#   # Single area FD index (numbers-based)
#   list(
#     indextype = "FD",
#     areas = c(1),  # Single area
#     indexYears = c(1:15), #1:5,  # Historical years only to test the bug
#     q_hist_bounds = c(0.0001, 0.0003),
#     q_proj_bounds = c(0.0001, 0.0003),
#     hyperstability_hist_bounds = c(0.9, 1.1),
#     hyperstability_proj_bounds = c(0.9, 1.1),
#     obsError_CV_hist_bounds = c(0.2, 0.3),
#     obsError_CV_proj_bounds = c(0.2, 0.3)
#   ),
#
#   # Multi-area FD index (forces manual calculation also for biomass)
#   list(
#     indextype = "FD",
#     areas = c(1, 2),  # Multi-area forces manual calculation path
#     indexYears = c(1:15), #c(2, 4, 6),  # Subset of years
#     q_hist_bounds = c(0.0002, 0.0005),
#     q_proj_bounds = c(0.0002, 0.0005),
#     hyperstability_hist_bounds = c(0.8, 1.2),
#     hyperstability_proj_bounds = c(0.8, 1.2),
#     obsError_CV_hist_bounds = c(0.15, 0.25),
#     obsError_CV_proj_bounds = c(0.15, 0.25)
#   )
# )
#
# cpue_numbers_single@selectivity_hist_list <- list()
# cpue_numbers_single@selectivity_proj_list <- list()
#
# # Strategy for this test (same as before but different name)
# simpleMP_single_bugtest <- function(phase, dataObject) {
#   for(r in 1:NROW(dataObject)) assign(names(dataObject)[r], dataObject[[r]])
#
#   if(phase==1) {
#     combined_data <- list()
#
#     if(!is.null(IndexObj)) {
#       index_result <- calculate_single_Index(dataObject)
#       for(col_name in names(index_result)) {
#         combined_data[[col_name]] <- index_result[[col_name]]
#       }
#     }
#
#     return(combined_data)
#   }
#
#   if(phase==2) return(list())
#   if(phase==3) {
#     year = rep(j, areas)
#     iteration = rep(k, areas)
#     area = 1:areas
#     Flocal = rep(0.15, areas)
#
#     return(list(year=year, iteration=iteration, area=area, Flocal=Flocal))
#   }
# }
#
# strategy_single_bugtest <- new("Strategy")
# strategy_single_bugtest@title <- "Single Fleet Bug Test"
# strategy_single_bugtest@projectionYears <- 5
# strategy_single_bugtest@projectionName <- "simpleMP_single_bugtest"
# strategy_single_bugtest@projectionParams <- list()
#
#
# # ============================================================================
# # RUN THE BUG TEST
# # ============================================================================
#
# cat("Testing single fleet with numbers-based FD indices...\n")
# cat("This should trigger the bug if the fix was not applied correctly.\n\n")
#
# # This test should expose the bug if the selectivity fix did not work
# tryCatch({
#   result_single_bugtest <- runProjection(
#     LifeHistoryObj = lh_obj,
#     TimeAreaObj = ta,
#     HistFisheryObj = hist_fishery,
#     ProFisheryObj_list = proj_fishery_list,
#     StrategyObj = strategy_single_bugtest,
#     StochasticObj = stochastic_obj,
#     MultifleetObj = NULL,  # Single fleet mode
#     IndexObj = cpue_numbers_single,  # Numbers-based FD indices
#     wd = getwd(),
#     fileName = "single_fleet_bugtest",
#     seed = validation_seed,
#     doPlot = FALSE,
#     doDiagnostic = FALSE
#   )
#
#
#   cat("OKEY SUCCESS: Single fleet numbers FD test completed without errors\n")
#   cat("The selectivity bug fix is working correctly.\n")
#
#   #examine results
#   result_single_bugtest <- readProjection(getwd(), "single_fleet_bugtest")
#
#   #check that indices were generated
#   index_cols <- grep("^CPUE_\\d+$", names(result_single_bugtest$HCR$decisionData), value = TRUE)
#   cat("Generated indices:", paste(index_cols, collapse = ", "), "\n")
#
#   #check some values
#   for(idx_col in index_cols) {
#     values <- result_single_bugtest$HCR$decisionData[[idx_col]]
#     valid_values <- values[!is.na(values)]
#
#     if(length(valid_values) > 0) {
#       cat(sprintf("  %s: %d valid values, range [%.2e, %.2e]\n",
#                   idx_col, length(valid_values), min(valid_values), max(valid_values)))
#     } else {
#       cat(sprintf("  %s: No valid values (potential issue)\n", idx_col))
#     }
#   }
#
#   }, error = function(e) {
#     cat("ERROR: Single fleet numbers FD test failed\n")
#     cat("Error message:", e$message, "\n")
#     cat("This confirms the selectivity bug exists and needs to be fixed.\n")
#
#     # Check if it's the specific selectivity error
#     if(grepl("object.*selectivity.*not found", e$message, ignore.case = TRUE)) {
#       cat("\nThis appears to be the selectivity assignment bug.\n")
#       cat("I added the missing 'else' for single fleet FD indices.\n")
#     }
#   })
#
# # ============================================================================
# # VERIFICATION OF BUG FIX
# # ============================================================================
#
# # Only run verification if the test succeeded
# if(exists("result_single_bugtest")) {
#   cat("\n=== SINGLE FLEET BUG FIX VERIFICATION ===\n")
#
#   # Compare single fleet vs multifleet(1-fleet) indices
#   # Both should give similar results if the fix is correct
#
#   decision_data_single <- result_single_bugtest$HCR$decisionData
#
#   # Check index structure
#   single_indices <- grep("^CPUE_\\d+$", names(decision_data_single), value = TRUE)
#   cat("Single fleet indices found:", length(single_indices), "\n")
#
#   # Check that they have correct metadata
#   for(idx in single_indices) {
#     indextype_col <- paste0(idx, "_indextype")
#     fleet_id_col <- paste0(idx, "_fleet_id")
#
#     if(indextype_col %in% names(decision_data_single)) {
#       indextype <- unique(decision_data_single[[indextype_col]])[1]
#       fleet_id <- if(fleet_id_col %in% names(decision_data_single)) {
#         unique(decision_data_single[[fleet_id_col]])[1]
#       } else "Missing"
#
#       cat(sprintf("  %s: Type=%s, Fleet_ID=%s\n", idx, indextype, fleet_id))
#     }
#   }
#
#   cat("\nBug fix verification complete.\n")
# } else {
#   cat("\nCannot verify bug fix - test failed to run.\n")
# }
#






# ============================================================================
# IMPROVED PLOTTING FUNCTIONS FOR SINGLE FLEET AND MULTIFLEET
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# ============================================================================
# MAIN PLOTTING FUNCTION
# ============================================================================

plot_fishery_dynamics <- function(simulation_result,
                                  save_plots = FALSE,
                                  output_dir = getwd(),
                                  plot_prefix = "fishery_dynamics") {

  #extract data
  dynamics <- simulation_result$dynamics
  is_multifleet <- !is.null(dynamics$multifleet)

  #get dimensions
  years <- dim(dynamics$SB)[1]
  iterations <- dim(dynamics$SB)[2]
  areas <- dim(dynamics$SB)[3]

  #time axis setup
  sim_years <- 1:years
  user_years <- sim_years - 1  # convert to user-friendly years
  historical_end <- simulation_result$TimeAreaObj@historicalYears + 1

  cat("creating plots for:", ifelse(is_multifleet, "Multifleet", "Single fleet"), "simulation\n")
  cat("dimensions: Years =", years, ", Iterations =", iterations, ", Areas =", areas, "\n")

  if(is_multifleet) {
    nfleets <- dynamics$multifleet$nfleets
    cat("Number of fleets:", nfleets, "\n")
  }

  # Create individual plots
  plots <- list()

  # 1. Spawning Biomass by Area
  plots$SB <- create_SB_plot(dynamics, sim_years, user_years, historical_end, areas, is_multifleet)

  # 2. Vulnerable Biomass by Area
  plots$VB <- create_VB_plot(dynamics, sim_years, user_years, historical_end, areas, is_multifleet)

  # 3. Fishing Mortality by Area and Fleet
  plots$F <- create_F_plot(dynamics, sim_years, user_years, historical_end, areas, is_multifleet)

  # 4. SPR (population-level)
  plots$SPR <- create_SPR_plot(dynamics, sim_years, user_years, historical_end, is_multifleet)

  # 5. Catch in Weight by Area and Fleet
  plots$Catch <- create_Catch_plot(dynamics, sim_years, user_years, historical_end, areas, is_multifleet)

  # 6. Recruitment by Area
  plots$RecN <- create_RecN_plot(dynamics, sim_years, user_years, historical_end, areas, is_multifleet)

  # 7. Index observations (if available)
  if(!is.null(simulation_result$HCR$decisionData)) {
    plots$Indices <- plotIndex_tibble_enhanced(simulation_result$HCR$decisionData)
  }

  # 8. NEW: Catch Observations (if available)
  if(!is.null(simulation_result$HCR$decisionData)) {
    catch_obs_plot <- create_catch_observations_plot(simulation_result$HCR$decisionData,
                                                     user_years, historical_end)
    if(!is.null(catch_obs_plot)) {
      plots$CatchObs <- catch_obs_plot
    }
  }

  # 9. NEW: True vs Observed Catch Comparison (if available)
  if(!is.null(simulation_result$HCR$decisionData)) {
    catch_comparison_plot <- create_catch_comparison_plot(simulation_result$HCR$decisionData,
                                                          user_years, historical_end)
    if(!is.null(catch_comparison_plot)) {
      plots$CatchComparison <- catch_comparison_plot
    }
  }



  # Save plots if requested
  if(save_plots) {
    for(plot_name in names(plots)) {
      filename <- file.path(output_dir, paste0(plot_prefix, "_", plot_name, ".jpeg"))
      ggsave(filename, plots[[plot_name]], width = 14, height = 10, dpi = 300)
      cat("Saved:", filename, "\n")
    }
  }

  return(plots)
}

# ============================================================================
# INDIVIDUAL PLOT CREATION FUNCTIONS
# ============================================================================

create_SB_plot <- function(dynamics, sim_years, user_years, historical_end, areas, is_multifleet) {

  plot_data <- data.frame()

  for(area in 1:areas) {
    #calculate median across iterations
    sb_median <- apply(dynamics$SB[, , area], 1, median, na.rm = TRUE)
    sb_q25 <- apply(dynamics$SB[, , area], 1, quantile, 0.25, na.rm = TRUE)
    sb_q75 <- apply(dynamics$SB[, , area], 1, quantile, 0.75, na.rm = TRUE)

    area_data <- data.frame(
      year = sim_years,
      user_year = user_years,
      median = sb_median,
      q25 = sb_q25,
      q75 = sb_q75,
      area = paste("Area", area)
    )

    plot_data <- rbind(plot_data, area_data)
  }

  p <- ggplot(plot_data, aes(x = year)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3, fill = "steelblue") +
    geom_line(aes(y = median), color = "steelblue", size = 1.2) +
    geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
    facet_wrap(~ area, scales = "free_y") +
    scale_x_continuous(breaks = sim_years, labels = user_years) +
    labs(title = "Spawning Biomass by Area",
         x = "Year",
         y = "Spawning Biomass",
         subtitle = "Median with 25th-75th percentile range") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}

create_VB_plot <- function(dynamics, sim_years, user_years, historical_end, areas, is_multifleet) {

  plot_data <- data.frame()

  for(area in 1:areas) {
    vb_median <- apply(dynamics$VB[, , area], 1, median, na.rm = TRUE)
    vb_q25 <- apply(dynamics$VB[, , area], 1, quantile, 0.25, na.rm = TRUE)
    vb_q75 <- apply(dynamics$VB[, , area], 1, quantile, 0.75, na.rm = TRUE)

    area_data <- data.frame(
      year = sim_years,
      user_year = user_years,
      median = vb_median,
      q25 = vb_q25,
      q75 = vb_q75,
      area = paste("Area", area)
    )

    plot_data <- rbind(plot_data, area_data)
  }

  p <- ggplot(plot_data, aes(x = year)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3, fill = "darkgreen") +
    geom_line(aes(y = median), color = "darkgreen", size = 1.2) +
    geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
    facet_wrap(~ area, scales = "free_y") +
    scale_x_continuous(breaks = sim_years, labels = user_years) +
    labs(title = "Vulnerable Biomass by Area",
         x = "Year",
         y = "Vulnerable Biomass",
         subtitle = "Median with 25th-75th percentile range") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}

create_F_plot <- function(dynamics, sim_years, user_years, historical_end, areas, is_multifleet) {

  plot_data <- data.frame()

  if(is_multifleet) {
    #plot fleet-specific F values
    nfleets <- dynamics$multifleet$nfleets

    for(area in 1:areas) {
      for(fleet in 1:nfleets) {
        f_median <- apply(dynamics$multifleet$Ftotal_by_fleet[, , area, fleet], 1, median, na.rm = TRUE)
        f_q25 <- apply(dynamics$multifleet$Ftotal_by_fleet[, , area, fleet], 1, quantile, 0.25, na.rm = TRUE)
        f_q75 <- apply(dynamics$multifleet$Ftotal_by_fleet[, , area, fleet], 1, quantile, 0.75, na.rm = TRUE)

        fleet_data <- data.frame(
          year = sim_years,
          user_year = user_years,
          median = f_median,
          q25 = f_q25,
          q75 = f_q75,
          panel = paste("Area", area, "- Fleet", fleet),
          fleet = paste("Fleet", fleet)
        )

        plot_data <- rbind(plot_data, fleet_data)
      }
    }

    p <- ggplot(plot_data, aes(x = year, color = fleet, fill = fleet)) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2) +
      geom_line(aes(y = median), size = 1.2) +
      geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
      facet_wrap(~ panel, scales = "free_y") +
      scale_x_continuous(breaks = sim_years, labels = user_years) +
      labs(title = "Fishing Mortality by Area and Fleet",
           x = "Year",
           y = "Fishing Mortality (F)",
           subtitle = "Median with 25th-75th percentile range") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  } else {
    #single fleet - total F by area
    for(area in 1:areas) {
      f_median <- apply(dynamics$Ftotal[, , area], 1, median, na.rm = TRUE)
      f_q25 <- apply(dynamics$Ftotal[, , area], 1, quantile, 0.25, na.rm = TRUE)
      f_q75 <- apply(dynamics$Ftotal[, , area], 1, quantile, 0.75, na.rm = TRUE)

      area_data <- data.frame(
        year = sim_years,
        user_year = user_years,
        median = f_median,
        q25 = f_q25,
        q75 = f_q75,
        area = paste("Area", area)
      )

      plot_data <- rbind(plot_data, area_data)
    }

    p <- ggplot(plot_data, aes(x = year)) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3, fill = "orange") +
      geom_line(aes(y = median), color = "orange", size = 1.2) +
      geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
      facet_wrap(~ area, scales = "free_y") +
      scale_x_continuous(breaks = sim_years, labels = user_years) +
      labs(title = "Fishing Mortality by Area",
           x = "Year",
           y = "Fishing Mortality (F)",
           subtitle = "Median with 25th-75th percentile range") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  return(p)
}

create_SPR_plot <- function(dynamics, sim_years, user_years, historical_end, is_multifleet) {

  #SPR is population-level (not area-specific)
  spr_median <- apply(dynamics$SPR, 1, median, na.rm = TRUE)
  spr_q25 <- apply(dynamics$SPR, 1, quantile, 0.25, na.rm = TRUE)
  spr_q75 <- apply(dynamics$SPR, 1, quantile, 0.75, na.rm = TRUE)

  plot_data <- data.frame(
    year = sim_years,
    user_year = user_years,
    median = spr_median,
    q25 = spr_q25,
    q75 = spr_q75
  )

  p <- ggplot(plot_data, aes(x = year)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3, fill = "purple") +
    geom_line(aes(y = median), color = "purple", size = 1.2) +
    geom_hline(yintercept = 0.3, linetype = "dotted", color = "red", alpha = 0.7) +
    geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = sim_years, labels = user_years) +
    labs(title = "Spawning Potential Ratio (SPR)",
         x = "Year",
         y = "SPR",
         subtitle = "Median with 25th-75th percentile range (dotted line = SPR 30%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p)
}

create_Catch_plot <- function(dynamics, sim_years, user_years, historical_end, areas, is_multifleet) {

  plot_data <- data.frame()

  if(is_multifleet) {
    #plot fleet-specific catches
    nfleets <- dynamics$multifleet$nfleets

    for(area in 1:areas) {
      for(fleet in 1:nfleets) {
        catch_median <- apply(dynamics$multifleet$catchB_by_fleet[, , area, fleet], 1, median, na.rm = TRUE)
        catch_q25 <- apply(dynamics$multifleet$catchB_by_fleet[, , area, fleet], 1, quantile, 0.25, na.rm = TRUE)
        catch_q75 <- apply(dynamics$multifleet$catchB_by_fleet[, , area, fleet], 1, quantile, 0.75, na.rm = TRUE)

        fleet_data <- data.frame(
          year = sim_years,
          user_year = user_years,
          median = catch_median,
          q25 = catch_q25,
          q75 = catch_q75,
          panel = paste("Area", area, "- Fleet", fleet),
          fleet = paste("Fleet", fleet)
        )

        plot_data <- rbind(plot_data, fleet_data)
      }
    }

    p <- ggplot(plot_data, aes(x = year, color = fleet, fill = fleet)) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2) +
      geom_line(aes(y = median), size = 1.2) +
      geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
      facet_wrap(~ panel, scales = "free_y") +
      scale_x_continuous(breaks = sim_years, labels = user_years) +
      labs(title = "Catch Biomass by Area and Fleet",
           x = "Year",
           y = "Catch Biomass",
           subtitle = "Median with 25th-75th percentile range") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

  } else {
    #single fleet - total catch by area
    for(area in 1:areas) {
      catch_median <- apply(dynamics$catchB[, , area], 1, median, na.rm = TRUE)
      catch_q25 <- apply(dynamics$catchB[, , area], 1, quantile, 0.25, na.rm = TRUE)
      catch_q75 <- apply(dynamics$catchB[, , area], 1, quantile, 0.75, na.rm = TRUE)

      area_data <- data.frame(
        year = sim_years,
        user_year = user_years,
        median = catch_median,
        q25 = catch_q25,
        q75 = catch_q75,
        area = paste("Area", area)
      )

      plot_data <- rbind(plot_data, area_data)
    }

    p <- ggplot(plot_data, aes(x = year)) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3, fill = "brown") +
      geom_line(aes(y = median), color = "brown", size = 1.2) +
      geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
      facet_wrap(~ area, scales = "free_y") +
      scale_x_continuous(breaks = sim_years, labels = user_years) +
      labs(title = "Catch Biomass by Area",
           x = "Year",
           y = "Catch Biomass",
           subtitle = "Median with 25th-75th percentile range") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  return(p)
}

create_RecN_plot <- function(dynamics, sim_years, user_years, historical_end, areas, is_multifleet) {

  #recruitment is population-level but can be split by area based on recArea
  if("recN" %in% names(dynamics)) {
    recn_median <- apply(dynamics$recN, 1, median, na.rm = TRUE)
    recn_q25 <- apply(dynamics$recN, 1, quantile, 0.25, na.rm = TRUE)
    recn_q75 <- apply(dynamics$recN, 1, quantile, 0.75, na.rm = TRUE)

    plot_data <- data.frame(
      year = sim_years,
      user_year = user_years,
      median = recn_median,
      q25 = recn_q25,
      q75 = recn_q75
    )

    p <- ggplot(plot_data, aes(x = year)) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.3, fill = "cyan") +
      geom_line(aes(y = median), color = "cyan4", size = 1.2) +
      geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = sim_years, labels = user_years) +
      labs(title = "Recruitment (Total)",
           x = "Year",
           y = "Recruitment (Numbers)",
           subtitle = "Median with 25th-75th percentile range") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    #create empty plot if recruitment data not available
    p <- ggplot() +
      geom_text(aes(x = 0.5, y = 0.5, label = "Recruitment data not available"), size = 6) +
      theme_void() +
      labs(title = "Recruitment Data Not Available")
  }

  return(p)
}

create_catch_observations_plot <- function(decision_data, user_years, historical_end) {

  #detect if we have catch observations
  catch_obs_cols <- grep("(^observed_catch$|fleet_\\d+_observed_catch$)", names(decision_data), value = TRUE)

  if(length(catch_obs_cols) == 0) {
    return(NULL)  # No catch observations found
  }

  #determine if multifleet or single fleet
  is_multifleet <- any(grepl("fleet_\\d+_", catch_obs_cols))

  if(is_multifleet) {
    return(create_multifleet_catch_obs_plot(decision_data, user_years, historical_end))
  } else {
    return(create_single_fleet_catch_obs_plot(decision_data, user_years, historical_end))
  }
}

create_single_fleet_catch_obs_plot <- function(decision_data, user_years, historical_end) {

  #check if we have the required columns
  if(!"observed_catch" %in% names(decision_data)) {
    return(NULL)
  }

  #prepare data
  plot_data <- decision_data %>%
    filter(!is.na(observed_catch)) %>%
    mutate(user_year = j - 1,
           iteration_label = paste("Iter", k))

  if(nrow(plot_data) == 0) {
    return(NULL)
  }

  #calculate median
  median_data <- plot_data %>%
    group_by(user_year) %>%
    summarise(median_catch = median(observed_catch, na.rm = TRUE), .groups = "drop")

  #create plot
  p <- ggplot(plot_data, aes(x = user_year, y = observed_catch)) +
    geom_line(aes(group = iteration_label, color = iteration_label), alpha = 0.7, size = 0.8) +
    geom_point(aes(color = iteration_label), alpha = 0.8, size = 1.5) +
    geom_line(data = median_data, aes(y = median_catch),
              color = "black", size = 1.5, linetype = "solid") +
    geom_point(data = median_data, aes(y = median_catch),
               color = "black", size = 2) +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_x_continuous(breaks = unique(plot_data$user_year)) +
    labs(title = "Single Fleet Catch Observations",
         subtitle = "Observed catch with reporting and observation error (black line = median)",
         x = "Year",
         y = "Observed Catch",
         color = "Iteration") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")

  return(p)
}

create_multifleet_catch_obs_plot <- function(decision_data, user_years, historical_end) {

  #find fleet-specific catch columns
  fleet_cols <- grep("fleet_\\d+_observed_catch$", names(decision_data), value = TRUE)

  if(length(fleet_cols) == 0) {
    return(NULL)
  }

  #prepare data for all fleets
  plot_data <- data.frame()

  for(col in fleet_cols) {
    #extract fleet number
    fleet_num <- gsub("fleet_(\\d+)_observed_catch", "\\1", col)

    #get data for this fleet
    fleet_data <- decision_data %>%
      filter(!is.na(!!sym(col))) %>%
      select(j, k, !!sym(col)) %>%
      mutate(user_year = j - 1,
             iteration = k,
             fleet = paste("Fleet", fleet_num),
             observed_catch = !!sym(col),
             iteration_label = paste("Iter", k))

    plot_data <- rbind(plot_data, fleet_data)
  }

  if(nrow(plot_data) == 0) {
    return(NULL)
  }

  #calculate median by fleet
  median_data <- plot_data %>%
    group_by(user_year, fleet) %>%
    summarise(median_catch = median(observed_catch, na.rm = TRUE), .groups = "drop")

  #create plot
  p <- ggplot(plot_data, aes(x = user_year, y = observed_catch, color = fleet)) +
    geom_line(aes(group = interaction(fleet, iteration_label)), alpha = 0.6, size = 0.7) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_line(data = median_data, aes(y = median_catch, color = fleet),
              size = 1.5, linetype = "solid") +
    geom_point(data = median_data, aes(y = median_catch, color = fleet),
               size = 2.5, shape = 15) +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_wrap(~ fleet, scales = "free_y") +
    scale_x_continuous(breaks = function(x) pretty(x, n = 6)) +
    labs(title = "Multifleet Catch Observations",
         subtitle = "Fleet-specific observed catches (thick lines = median)",
         x = "Year",
         y = "Observed Catch",
         color = "Fleet") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  return(p)
}

create_catch_comparison_plot <- function(decision_data, user_years, historical_end) {

  # Check what type of catch data we have
  has_single_fleet <- all(c("true_catch", "observed_catch") %in% names(decision_data))
  has_multifleet <- any(grepl("fleet_\\d+_true_catch", names(decision_data))) &&
    any(grepl("fleet_\\d+_observed_catch", names(decision_data)))

  if(!has_single_fleet && !has_multifleet) {
    return(NULL)  # No suitable data for comparison
  }

  if(has_multifleet) {
    return(create_multifleet_catch_comparison_plot(decision_data, user_years, historical_end))
  } else {
    return(create_single_fleet_catch_comparison_plot(decision_data, user_years, historical_end))
  }
}

create_single_fleet_catch_comparison_plot <- function(decision_data, user_years, historical_end) {

  # Prepare comparison data
  plot_data <- decision_data %>%
    filter(!is.na(true_catch) & !is.na(observed_catch)) %>%
    mutate(user_year = j - 1,
           iteration_label = paste("Iter", k)) %>%
    select(user_year, iteration_label, k, true_catch, observed_catch) %>%
    pivot_longer(cols = c(true_catch, observed_catch),
                 names_to = "catch_type", values_to = "catch_value") %>%
    mutate(catch_type = case_when(
      catch_type == "true_catch" ~ "True Catch",
      catch_type == "observed_catch" ~ "Observed Catch"
    ))

  if(nrow(plot_data) == 0) {
    return(NULL)
  }

  # Calculate median by type
  median_data <- plot_data %>%
    group_by(user_year, catch_type) %>%
    summarise(median_catch = median(catch_value, na.rm = TRUE), .groups = "drop")

  # Create comparison plot
  p <- ggplot(plot_data, aes(x = user_year, y = catch_value, color = catch_type)) +
    geom_line(aes(group = interaction(catch_type, iteration_label)), alpha = 0.5, size = 0.6) +
    geom_point(alpha = 0.6, size = 1) +
    geom_line(data = median_data, aes(y = median_catch, color = catch_type),
              size = 1.5, linetype = "solid") +
    geom_point(data = median_data, aes(y = median_catch, color = catch_type),
               size = 2.5, shape = 15) +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_color_manual(values = c("True Catch" = "darkblue", "Observed Catch" = "orange")) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
    labs(title = "True vs Observed Catch Comparison",
         subtitle = "Single fleet - showing effect of reporting rates and observation error",
         x = "Year",
         y = "Catch",
         color = "Catch Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

  return(p)
}

create_multifleet_catch_comparison_plot <- function(decision_data, user_years, historical_end) {

  #find available fleet data
  true_cols <- grep("fleet_\\d+_true_catch$", names(decision_data), value = TRUE)
  obs_cols <- grep("fleet_\\d+_observed_catch$", names(decision_data), value = TRUE)

  #match fleet numbers
  fleet_nums <- unique(c(
    gsub("fleet_(\\d+)_true_catch", "\\1", true_cols),
    gsub("fleet_(\\d+)_observed_catch", "\\1", obs_cols)
  ))

  #prepare data for all available fleets
  plot_data <- data.frame()

  for(fleet_num in fleet_nums) {
    true_col <- paste0("fleet_", fleet_num, "_true_catch")
    obs_col <- paste0("fleet_", fleet_num, "_observed_catch")

    if(true_col %in% names(decision_data) && obs_col %in% names(decision_data)) {

      fleet_data <- decision_data %>%
        filter(!is.na(!!sym(true_col)) & !is.na(!!sym(obs_col))) %>%
        mutate(user_year = j - 1,
               iteration = k,
               fleet = paste("Fleet", fleet_num)) %>%
        select(user_year, iteration, fleet,
               true_catch = !!sym(true_col),
               observed_catch = !!sym(obs_col)) %>%
        pivot_longer(cols = c(true_catch, observed_catch),
                     names_to = "catch_type", values_to = "catch_value") %>%
        mutate(catch_type = case_when(
          catch_type == "true_catch" ~ "True Catch",
          catch_type == "observed_catch" ~ "Observed Catch"
        ))

      plot_data <- rbind(plot_data, fleet_data)
    }
  }

  if(nrow(plot_data) == 0) {
    return(NULL)
  }

  #calculate median by fleet and type
  median_data <- plot_data %>%
    group_by(user_year, fleet, catch_type) %>%
    summarise(median_catch = median(catch_value, na.rm = TRUE), .groups = "drop")

  #create comparison plot
  p <- ggplot(plot_data, aes(x = user_year, y = catch_value, color = catch_type)) +
    geom_line(aes(group = interaction(catch_type, iteration)), alpha = 0.4, size = 0.5) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_line(data = median_data, aes(y = median_catch, color = catch_type),
              size = 1.2, linetype = "solid") +
    geom_point(data = median_data, aes(y = median_catch, color = catch_type),
               size = 2, shape = 15) +
    geom_vline(xintercept = historical_end - 1, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_wrap(~ fleet, scales = "free_y") +
    scale_color_manual(values = c("True Catch" = "darkblue", "Observed Catch" = "orange")) +
    scale_x_continuous(breaks = function(x) pretty(x, n = 5)) +
    labs(title = "Multifleet True vs Observed Catch Comparison",
         subtitle = "Fleet-specific comparison showing observation model effects",
         x = "Year",
         y = "Catch",
         color = "Catch Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

  return(p)
}


plot_multifleet_with_catch_obs <- function() {

  #assuming we have 'result_multifleet'
  if(exists("result_multifleet")) {

    #create plots including catch observations
    enhanced_plots <- plot_fishery_dynamics(result_multifleet,
                                            save_plots = FALSE,
                                            plot_prefix = "multifleet_enhanced")

    #dDisplay the new catch observation plots
    if("CatchObs" %in% names(enhanced_plots)) {
      print("Catch Observations Plot:")
      print(enhanced_plots$CatchObs)
    }

    if("CatchComparison" %in% names(enhanced_plots)) {
      print("True vs Observed Catch Comparison:")
      print(enhanced_plots$CatchComparison)
    }

    # We can also view other plots
    # enhanced_plots$SB        # Spawning biomass
    # enhanced_plots$F         # Fishing mortality
    # enhanced_plots$Indices   # Index observations

    return(enhanced_plots)
  } else {
    cat("result_multifleet not found. Please run your multifleet simulation first.\n")
    return(NULL)
  }
}


#quick function to create a grid of key plots including catch obs
create_summary_plot_grid <- function(simulation_result, save_file = NULL) {

  plots <- plot_fishery_dynamics(simulation_result, save_plots = FALSE)

  #select key plots for summary
  key_plots <- list()

  if("SB" %in% names(plots)) key_plots$SB <- plots$SB
  if("F" %in% names(plots)) key_plots$F <- plots$F
  if("CatchObs" %in% names(plots)) key_plots$CatchObs <- plots$CatchObs
  if("CatchComparison" %in% names(plots)) key_plots$CatchComparison <- plots$CatchComparison

  #create grid
  if(length(key_plots) > 0) {
    grid_plot <- gridExtra::grid.arrange(grobs = key_plots, ncol = 2)

    if(!is.null(save_file)) {
      ggsave(save_file, grid_plot, width = 16, height = 12, dpi = 300)
      cat("Saved summary plot grid to:", save_file, "\n")
    }

    return(grid_plot)
  }

  return(NULL)
}



# ============================================================================
# ENHANCED INDEX PLOTTING (FROM PREVIOUS FUNCTION)
# ============================================================================

plotIndex_tibble_enhanced <- function(tibble_data, save_plot = FALSE,
                                      filename = "index_plot.jpeg") {

  #extract metadata
  title <- unique(tibble_data$title)[1]
  historical_end <- unique(tibble_data$historical_end)[1]
  is_multifleet <- unique(tibble_data$is_multifleet)[1]

  #enhanced pattern to detect both single fleet and multifleet indices
  if(is_multifleet) {
    index_columns <- grep("^(CPUE_\\d+(_Fleet_\\d+)?|Survey_\\d+)$", names(tibble_data), value = TRUE)
  } else {
    index_columns <- grep("^(CPUE_|Survey_)\\d+$", names(tibble_data), value = TRUE)
  }

  if(length(index_columns) == 0) {
    return(ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "No index data available"), size = 6) + theme_void())
  }

  all_data <- data.frame()

  for(index_col in index_columns) {
    indextype_col <- paste0(index_col, "_indextype")
    areas_col <- paste0(index_col, "_areas")
    indexyears_col <- paste0(index_col, "_indexYears")
    fleet_id_col <- paste0(index_col, "_fleet_id")

    if(!all(c(indextype_col, areas_col, indexyears_col) %in% names(tibble_data))) {
      next
    }

    index_years_str <- unique(tibble_data[[indexyears_col]])[1]
    if(is.na(index_years_str)) next

    index_years <- as.numeric(unlist(strsplit(as.character(index_years_str), "_")))
    areas_str <- unique(tibble_data[[areas_col]])[1]

    #enhanced panel naming
    if(fleet_id_col %in% names(tibble_data)) {
      fleet_id <- unique(tibble_data[[fleet_id_col]])[1]
      if(!is.na(fleet_id)) {
        panel_name <- paste(index_col, "- Fleet", fleet_id, "- Area(s)", areas_str)
      } else {
        panel_name <- paste(index_col, "- Area(s)", areas_str)
      }
    } else {
      panel_name <- paste(index_col, "- Area(s)", areas_str)
    }

    index_data <- tibble_data %>%
      select(k, j, all_of(index_col)) %>%
      rename(iteration = k, year = j, value = !!sym(index_col)) %>%
      mutate(panel = panel_name, type = "Iterations", iteration_label = paste0("Iter_", iteration))

    median_data <- index_data %>%
      group_by(year, panel) %>%
      summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
      mutate(type = "Median", iteration = "Median", iteration_label = "Median")

    combined_data <- rbind(
      index_data %>% select(year, value, panel, type, iteration_label),
      median_data %>% select(year, value, panel, type, iteration_label)
    )

    all_data <- rbind(all_data, combined_data)
  }

  #Y-axis labeling
  use_weight <- unique(tibble_data$useWeight)[1]
  final_indextype <- unique(tibble_data$final_indextype)[1]

  if(final_indextype == "FD") {
    y_label <- if(use_weight) "CPUE (biomass)" else "CPUE (numbers)"
  } else if(final_indextype == "FI") {
    y_label <- if(use_weight) "Survey (biomass)" else "Survey (numbers)"
  } else {
    y_label <- if(use_weight) "Index (biomass)" else "Index (numbers)"
  }



  #enhanced title
  if(is_multifleet) {
    nfleets <- unique(tibble_data$nfleets)[1]
    plot_title <- paste(title, "(", nfleets, "fleets )")
  } else {
    plot_title <- title
  }

  sim_years <- sort(unique(all_data$year))
  user_years <- sim_years - 1

  p <- ggplot(all_data, aes(x = year, y = value)) +
    geom_line(data = subset(all_data, type == "Iterations"),
              aes(group = iteration_label), color = "steelblue", alpha = 0.6, size = 0.5) +
    geom_point(data = subset(all_data, type == "Iterations" & !is.na(value)),
               color = "steelblue", alpha = 0.7, size = 1) +
    geom_line(data = subset(all_data, type == "Median"),
              color = "black", size = 1.2) +
    geom_point(data = subset(all_data, type == "Median" & !is.na(value)),
               color = "black", size = 1.5) +
    facet_wrap(~ panel, scales = "free_y") +
    geom_vline(xintercept = historical_end, linetype = "dashed", color = "red") +
    scale_x_continuous(breaks = sim_years, labels = user_years) +
    labs(title = plot_title, x = "Year", y = y_label) +
    theme_minimal() +
    theme(strip.text = element_text(size = 9), plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))

  if(save_plot) {
    ggsave(filename, plot = p, width = 14, height = 10, dpi = 300)
  }

  return(p)
}

# ============================================================================
# HOW THE USER CAN USE THESE FUCNTIONS
# ============================================================================

#for single fleet:
plots_single <- plot_fishery_dynamics(result_single,
                                      save_plots = FALSE,
                                      plot_prefix = "single_fleet")

#for multifleet:
plots_multi <- plot_fishery_dynamics(result_multifleet,
                                     save_plots = FALSE,
                                     plot_prefix = "multifleet")

#view individual plots:
plots_single$SB      # Spawning biomass
plots_single$F       # Fishing mortality
plots_multi$Catch    # Fleet-specific catches
plots_multi$Indices  # Index observations



stop() #(testing with more years again - the example: result_single_bugtest )

#for single fleet (NOTE Need to fix plot scale when using numbers)
# plots_single_bugtest  <- plot_fishery_dynamics(result_single_bugtest,
#                                       save_plots = FALSE,
#                                       plot_prefix = "single_fleet_bugtest")
# plots_single_bugtest$SB      # Spawning biomass
# plots_single_bugtest$F       # Fishing mortality


# Save all plots including catch observations
enhanced_plots <- plot_fishery_dynamics(result_multifleet,
                                        save_plots = TRUE,
                                        output_dir = getwd(),
                                        plot_prefix = "multifleet_complete")

