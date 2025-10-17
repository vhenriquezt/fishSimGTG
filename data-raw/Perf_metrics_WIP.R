#Performance metrics

#Arguments

# 1) result_list: of simulation results from runProjection() or readProjection().
# Can be a single result or named list of results for comparison.


# 2) time_horizons: list with 'short', 'medium', 'long' as vectors of year indices.

#3) iav_threshold= default 15%

#4) by_fleet = f TRUE and multifleet data exists, calculate fleet-specific metrics

#5) by_area =  If TRUE, calculate area-specific metrics

#6 output = "table" returns data.frame

calculate_performance_metrics <- function(result_list,
                                          mp_names = NULL,
                                          time_horizons = NULL,
                                          iav_threshold = 0.15,
                                          by_fleet = TRUE,
                                          by_area = FALSE,
                                          output_format = "table") {


  #1. single result vs list of results
  if(!is.list(result_list) || "dynamics" %in% names(result_list)) {
    # Single result provided
    result_list <- list(MP1 = result_list)
    if(!is.null(mp_names) && length(mp_names) == 1) {
      names(result_list) <- mp_names
    }
  }

  #validate MP names
  if(is.null(names(result_list))) {
    names(result_list) <- paste0("MP", seq_along(result_list))
  }

  n_mps <- length(result_list)
  mp_names <- names(result_list)

  cat(sprintf("Calculating performance metrics for %d management procedure(s):\n", n_mps))
  for(i in seq_along(mp_names)) {
    cat(sprintf("  %d. %s\n", i, mp_names[i]))
  }
  cat("\n")

  #2. extract dims

  first_result <- result_list[[1]]

  #detect multifleet mode
  is_multifleet <- !is.null(first_result$dynamics$multifleet) &&
    !is.null(first_result$dynamics$multifleet$catchB_by_fleet)

  if(is_multifleet) {
    years <- dim(first_result$dynamics$multifleet$catchB_by_fleet)[1]
    iterations <- dim(first_result$dynamics$multifleet$catchB_by_fleet)[2]
    areas <- dim(first_result$dynamics$multifleet$catchB_by_fleet)[3]
    nfleets <- dim(first_result$dynamics$multifleet$catchB_by_fleet)[4]
    cat(sprintf("Multifleet mode detected: %d fleets, %d areas\n", nfleets, areas))
  } else {
    years <- dim(first_result$dynamics$catchB)[1]
    iterations <- dim(first_result$dynamics$catchB)[2]
    areas <- dim(first_result$dynamics$catchB)[3]
    nfleets <- 1
    cat("Single fleet mode detected\n")
  }

  historical_years <- first_result$TimeAreaObj@historicalYears
  projection_start <- historical_years + 2  # +1 for equilibrium, +1 for indexing
  projection_years <- years - projection_start + 1

  cat(sprintf("Dimensions: %d years total (%d historical + %d projection), %d iterations, %d areas\n\n",
              years, historical_years, projection_years, iterations, areas))


  #3. time horizon

  if(is.null(time_horizons)) {
    #default time horizons (relative to projection start)
    time_horizons <- list(
      short = projection_start:(projection_start + 4),      # Years 1-5 of projection
      medium = (projection_start + 5):(projection_start + 9), # Years 6-10
      long = (projection_start + 10):years                    # Years 11+
    )

    #adjust if projection period is shorter
    time_horizons$short <- time_horizons$short[time_horizons$short <= years]
    time_horizons$medium <- time_horizons$medium[time_horizons$medium <= years]
    time_horizons$long <- time_horizons$long[time_horizons$long <= years]

    #remove empty horizons
    time_horizons <- time_horizons[sapply(time_horizons, length) > 0]
  }

  cat("Time horizons defined:\n")
  for(horizon in names(time_horizons)) {
    cat(sprintf("  %s: years %d-%d (%d years)\n",
                horizon,
                min(time_horizons[[horizon]]) - projection_start + 1,
                max(time_horizons[[horizon]]) - projection_start + 1,
                length(time_horizons[[horizon]])))
  }
  cat("\n")

  #4. main calcs

  results_list <- list()

  for(mp_idx in seq_along(result_list)) {
    mp_name <- mp_names[mp_idx]
    result <- result_list[[mp_idx]]

    cat(sprintf("Processing MP: %s\n", mp_name))
    cat(strrep("=", 60), "\n")

    #extract data
    SB <- result$dynamics$SB

    if(is_multifleet) {
      catchB_by_fleet <- result$dynamics$multifleet$catchB_by_fleet
      #total catch (sum across fleets)
      catchB <- apply(catchB_by_fleet, c(1, 2, 3), sum)
    } else {
      catchB <- result$dynamics$catchB
      catchB_by_fleet <- NULL
    }

    #calculate SB at end of historical period (reference point)
    SB_final_historical <- SB[historical_years + 1, , ]  # [iterations, areas]

    #Metric 1: P(SB > SB_historical)

    cat("\nMetric 1: P(SB > SB_historical)\n")
    cat(strrep("-", 40), "\n")

    metric1_results <- list()

    for(horizon_name in names(time_horizons)) {
      horizon_years <- time_horizons[[horizon_name]]

      if(by_area) {
        # Area-specific probabilities
        for(area in 1:areas) {
          prob_above <- sapply(horizon_years, function(yr) {
            mean(SB[yr, , area] > SB_final_historical[, area])
          })

          metric1_results[[paste0(horizon_name, "_area", area)]] <- list(
            mean_prob = mean(prob_above),
            min_prob = min(prob_above),
            max_prob = max(prob_above),
            prob_by_year = prob_above
          )

          cat(sprintf("  %s (Area %d): Mean P = %.3f (range: %.3f - %.3f)\n",
                      horizon_name, area, mean(prob_above), min(prob_above), max(prob_above)))
        }
      } else {
        # Total across areas
        SB_total <- apply(SB[horizon_years, , , drop=FALSE], c(1, 2), sum)
        SB_historical_total <- apply(SB_final_historical, 1, sum)

        prob_above <- sapply(1:length(horizon_years), function(i) {
          mean(SB_total[i, ] > SB_historical_total)
        })

        metric1_results[[horizon_name]] <- list(
          mean_prob = mean(prob_above),
          min_prob = min(prob_above),
          max_prob = max(prob_above),
          prob_by_year = prob_above
        )

        cat(sprintf("  %s (Total): Mean P = %.3f (range: %.3f - %.3f)\n",
                    horizon_name, mean(prob_above), min(prob_above), max(prob_above)))
      }
    }

    #metric 2: mean catch
    cat("\nMetric 2: Mean Catch Over Projection Period\n")
    cat(strrep("-", 40), "\n")

    metric2_results <- list()

    for(horizon_name in names(time_horizons)) {
      horizon_years <- time_horizons[[horizon_name]]

      if(is_multifleet && by_fleet) {
        # Fleet-specific means
        for(fleet in 1:nfleets) {
          if(by_area) {
            for(area in 1:areas) {
              catches <- catchB_by_fleet[horizon_years, , area, fleet]
              mean_catch <- mean(catches)
              sd_catch <- sd(catches)
              cv_catch <- sd_catch / mean_catch

              metric2_results[[paste0(horizon_name, "_fleet", fleet, "_area", area)]] <- list(
                mean = mean_catch,
                sd = sd_catch,
                cv = cv_catch,
                median = median(catches),
                q25 = quantile(catches, 0.25),
                q75 = quantile(catches, 0.75)
              )

              cat(sprintf("  %s (Fleet %d, Area %d): Mean = %.2f (SD = %.2f, CV = %.3f)\n",
                          horizon_name, fleet, area, mean_catch, sd_catch, cv_catch))
            }
          } else {
            # Fleet total across areas
            catches <- apply(catchB_by_fleet[horizon_years, , , fleet, drop=FALSE], c(1, 2), sum)
            mean_catch <- mean(catches)
            sd_catch <- sd(catches)
            cv_catch <- sd_catch / mean_catch

            metric2_results[[paste0(horizon_name, "_fleet", fleet)]] <- list(
              mean = mean_catch,
              sd = sd_catch,
              cv = cv_catch,
              median = median(catches),
              q25 = quantile(catches, 0.25),
              q75 = quantile(catches, 0.75)
            )

            cat(sprintf("  %s (Fleet %d): Mean = %.2f (SD = %.2f, CV = %.3f)\n",
                        horizon_name, fleet, mean_catch, sd_catch, cv_catch))
          }
        }
      }

      # Total catch (across fleets and areas)
      if(by_area) {
        for(area in 1:areas) {
          catches <- catchB[horizon_years, , area]
          mean_catch <- mean(catches)
          sd_catch <- sd(catches)
          cv_catch <- sd_catch / mean_catch

          metric2_results[[paste0(horizon_name, "_total_area", area)]] <- list(
            mean = mean_catch,
            sd = sd_catch,
            cv = cv_catch,
            median = median(catches),
            q25 = quantile(catches, 0.25),
            q75 = quantile(catches, 0.75)
          )

          cat(sprintf("  %s (Total, Area %d): Mean = %.2f (SD = %.2f, CV = %.3f)\n",
                      horizon_name, area, mean_catch, sd_catch, cv_catch))
        }
      } else {
        catches <- apply(catchB[horizon_years, , , drop=FALSE], c(1, 2), sum)
        mean_catch <- mean(catches)
        sd_catch <- sd(catches)
        cv_catch <- sd_catch / mean_catch

        metric2_results[[paste0(horizon_name, "_total")]] <- list(
          mean = mean_catch,
          sd = sd_catch,
          cv = cv_catch,
          median = median(catches),
          q25 = quantile(catches, 0.25),
          q75 = quantile(catches, 0.75)
        )

        cat(sprintf("  %s (Total): Mean = %.2f (SD = %.2f, CV = %.3f)\n",
                    horizon_name, mean_catch, sd_catch, cv_catch))
      }
    }

    #metric 3: Cacth stability (IAV < threshold)

    cat("\nMetric 3: Catch Stability (P(IAV < 15%))\n")
    cat(strrep("-", 40), "\n")

    metric3_results <- list()

    # Function to calculate IAV
    calculate_iav <- function(catch_series) {
      # catch_series: matrix [years, iterations]
      n_years <- nrow(catch_series)
      if(n_years < 2) return(NA)

      # Calculate year-to-year changes for each iteration
      iav_by_iter <- sapply(1:ncol(catch_series), function(iter) {
        catches <- catch_series[, iter]
        changes <- abs(diff(catches)) / catches[-length(catches)]
        mean(changes, na.rm = TRUE)
      })

      return(iav_by_iter)
    }

    for(horizon_name in names(time_horizons)) {
      horizon_years <- time_horizons[[horizon_name]]

      if(length(horizon_years) < 2) {
        cat(sprintf("  %s: Skipped (need at least 2 years)\n", horizon_name))
        next
      }

      if(is_multifleet && by_fleet) {
        for(fleet in 1:nfleets) {
          if(by_area) {
            for(area in 1:areas) {
              catch_series <- catchB_by_fleet[horizon_years, , area, fleet]
              iav_values <- calculate_iav(catch_series)
              prob_stable <- mean(iav_values < iav_threshold, na.rm = TRUE)

              metric3_results[[paste0(horizon_name, "_fleet", fleet, "_area", area)]] <- list(
                prob_stable = prob_stable,
                mean_iav = mean(iav_values, na.rm = TRUE),
                median_iav = median(iav_values, na.rm = TRUE),
                threshold = iav_threshold
              )

              cat(sprintf("  %s (Fleet %d, Area %d): P(IAV < %.0f%%) = %.3f (Mean IAV = %.3f)\n",
                          horizon_name, fleet, area, iav_threshold*100, prob_stable,
                          mean(iav_values, na.rm = TRUE)))
            }
          } else {
            catch_series <- apply(catchB_by_fleet[horizon_years, , , fleet, drop=FALSE], c(1, 2), sum)
            iav_values <- calculate_iav(catch_series)
            prob_stable <- mean(iav_values < iav_threshold, na.rm = TRUE)

            metric3_results[[paste0(horizon_name, "_fleet", fleet)]] <- list(
              prob_stable = prob_stable,
              mean_iav = mean(iav_values, na.rm = TRUE),
              median_iav = median(iav_values, na.rm = TRUE),
              threshold = iav_threshold
            )

            cat(sprintf("  %s (Fleet %d): P(IAV < %.0f%%) = %.3f (Mean IAV = %.3f)\n",
                        horizon_name, fleet, iav_threshold*100, prob_stable,
                        mean(iav_values, na.rm = TRUE)))
          }
        }
      }

      #total catch stability
      if(by_area) {
        for(area in 1:areas) {
          catch_series <- catchB[horizon_years, , area]
          iav_values <- calculate_iav(catch_series)
          prob_stable <- mean(iav_values < iav_threshold, na.rm = TRUE)

          metric3_results[[paste0(horizon_name, "_total_area", area)]] <- list(
            prob_stable = prob_stable,
            mean_iav = mean(iav_values, na.rm = TRUE),
            median_iav = median(iav_values, na.rm = TRUE),
            threshold = iav_threshold
          )

          cat(sprintf("  %s (Total, Area %d): P(IAV < %.0f%%) = %.3f (Mean IAV = %.3f)\n",
                      horizon_name, area, iav_threshold*100, prob_stable,
                      mean(iav_values, na.rm = TRUE)))
        }
      } else {
        catch_series <- apply(catchB[horizon_years, , , drop=FALSE], c(1, 2), sum)
        iav_values <- calculate_iav(catch_series)
        prob_stable <- mean(iav_values < iav_threshold, na.rm = TRUE)

        metric3_results[[paste0(horizon_name, "_total")]] <- list(
          prob_stable = prob_stable,
          mean_iav = mean(iav_values, na.rm = TRUE),
          median_iav = median(iav_values, na.rm = TRUE),
          threshold = iav_threshold
        )

        cat(sprintf("  %s (Total): P(IAV < %.0f%%) = %.3f (Mean IAV = %.3f)\n",
                    horizon_name, iav_threshold*100, prob_stable,
                    mean(iav_values, na.rm = TRUE)))
      }
    }

    cat("\n")

    #save results

    results_list[[mp_name]] <- list(
      biomass_prob = metric1_results,
      mean_catch = metric2_results,
      catch_stability = metric3_results
    )
  }

  #output as a table

  if(output_format == "table") {
    # Convert nested list to data frame for easy comparison
    metrics_df <- data.frame()

    for(mp_name in names(results_list)) {
      mp_results <- results_list[[mp_name]]

      #biomass probability metrics
      for(metric_name in names(mp_results$biomass_prob)) {
        metrics_df <- rbind(metrics_df, data.frame(
          MP = mp_name,
          Metric = "P(SB > SB_historical)",
          Horizon = gsub("_area.*", "", metric_name),
          Fleet_Area = gsub(".*_(fleet.*|area.*|total)", "\\1", metric_name),
          Value = mp_results$biomass_prob[[metric_name]]$mean_prob,
          Min = mp_results$biomass_prob[[metric_name]]$min_prob,
          Max = mp_results$biomass_prob[[metric_name]]$max_prob,
          stringsAsFactors = FALSE
        ))
      }

      #mean catch metrics
      for(metric_name in names(mp_results$mean_catch)) {
        metrics_df <- rbind(metrics_df, data.frame(
          MP = mp_name,
          Metric = "Mean Catch",
          Horizon = gsub("_fleet.*|_area.*|_total.*", "", metric_name),
          Fleet_Area = gsub(".*_(fleet.*|area.*|total.*)", "\\1", metric_name),
          Value = mp_results$mean_catch[[metric_name]]$mean,
          Min = mp_results$mean_catch[[metric_name]]$q25,
          Max = mp_results$mean_catch[[metric_name]]$q75,
          stringsAsFactors = FALSE
        ))
      }

      #catch stability metrics
      for(metric_name in names(mp_results$catch_stability)) {
        metrics_df <- rbind(metrics_df, data.frame(
          MP = mp_name,
          Metric = "P(IAV < 15%)",
          Horizon = gsub("_fleet.*|_area.*|_total.*", "", metric_name),
          Fleet_Area = gsub(".*_(fleet.*|area.*|total.*)", "\\1", metric_name),
          Value = mp_results$catch_stability[[metric_name]]$prob_stable,
          Min = NA,
          Max = NA,
          stringsAsFactors = FALSE
        ))
      }
    }

    return(metrics_df)

  } else {
    #return nested list structure
    return(results_list)
  }
}




