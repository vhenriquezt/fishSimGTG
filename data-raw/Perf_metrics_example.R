
#==============================================================================
# Performance Metrics
#==============================================================================

# 1. P(SB > SB_historical): Probability spawning biomass exceeds last historical SB
# 2. P(Catch > Catch_historical): Probability catch exceeds last historical catch
# 3. Interannual Variability (IAV): Probability catch cariability < 15%


calculate_performance_metrics <- function(sim_result, mp_name) {

  # extract information
  n_hist <- sim_result$TimeAreaObj@historicalYears
  n_proj <- if(!is.null(sim_result$StrategyObj)) {
    sim_result$StrategyObj@projectionYears
  } else {
    0
  }
  n_areas <- sim_result$TimeAreaObj@areas
  n_iter <- dim(sim_result$dynamics$SB)[2]

  # year indices
  hist_end <- n_hist + 1  # last historical year (in j-index)
  proj_start <- hist_end + 1
  proj_end <- hist_end + n_proj

  cat("\n=== Processing MP:", mp_name, "===\n")
  cat("Historical years:", n_hist, "(j=2 to j=", hist_end, ")\n")
  cat("Projection years:", n_proj, "(j=", proj_start, "to j=", proj_end, ")\n")
  cat("Areas:", n_areas, "\n")
  cat("Iterations:", n_iter, "\n")

  # detect multifleet mode
  is_multifleet <- !is.null(sim_result$dynamics$multifleet)

  if(is_multifleet) {
    n_fleets <- dim(sim_result$dynamics$multifleet$catchB_by_fleet)[4] #[years, iter, area, fleet]
    cat("Multifleet mode: YES (", n_fleets, "fleets)\n")
  } else {
    cat("Multifleet mode: NO (single fleet)\n")
  }

  #============================================================================
  # METRIC 1: P(SB_year_proj / SB_final_historical > 1.0)
  #============================================================================

  # for each projection year, calculate the probability (across iterations)
  # that spawning biomass exceeds the final historical year

  #get spawning biomass (aggregate across all areas)
  SB_total <- apply(sim_result$dynamics$SB, c(1, 2), sum)  # [year, iter]

  #reference: last historical year spawning biomass for each iteration
  SB_hist_final <- SB_total[hist_end, ]  #vector of length n_iter

  # calculate ratio for each projection year
  sb_metrics <- data.frame()

  for(proj_year in proj_start:proj_end) {
    #if proj_start = 12, proj_end = 16, loop for years 12, 13, 14, 15, 16
    real_year <- proj_year - 1  # convert j-index to real year (e.g., j=12 -> year 11)

    # ratio for each iteration
    sb_ratio <- SB_total[proj_year, ] / SB_hist_final

    #SB_total[proj_year, ] extracts one row from SB_total
    #this gives SB for all iterations in this projection year
    #example: [110.2, 90.5, 120.3, 105.1, 95.8, 103.2] for proj_year = 12
    #SB_hist_final : one value for each ietrartion


    # calculate probability
    # mean (proportion of true values)
    # TRUE=1, FALSE=0, then calculate the mean
    # eg., [TRUE, FALSE, TRUE, TRUE, FALSE, TRUE]
    #      [1,0,1,1,0,1]/6=4/6=0.667
    prob_above_1 <- mean(sb_ratio > 1.0) #sb_ratio > 1.0 creates logical vector [TRUE, FALSE, TRUE]

    mean_ratio_val <- mean(sb_ratio)      # average ratio across all iterations
    median_ratio_val <- median(sb_ratio)  # median ratio across all iterations
    sd_ratio_val <- sd(sb_ratio)          # standard deviation across iterations
    min_ratio_val <- min(sb_ratio)        # minimum ratio across all iterations
    max_ratio_val <- max(sb_ratio)        # maximum ratio across all iterations



    # store results
    sb_metrics <- rbind(sb_metrics, data.frame(
      mp_name = mp_name,                                # mame of MP
      year = real_year,                                 # real year number (e.g., 11)
      proj_year_index = proj_year - proj_start + 1,    # projection year index (1, 2, 3, ...)
      metric = "P(SB > SB_hist)",                      # metric name
      probability = prob_above_1,                       # probability across iterations
      mean_ratio = mean_ratio_val,                      # mean ratio across iterations
      median_ratio = median_ratio_val,                  # median ratio across iterations
      sd_ratio = sd_ratio_val,                          # SD across iterations
      min_ratio = min_ratio_val,                        # min across iterations
      max_ratio = max_ratio_val,                        # max across iterations
      n_iterations = n_iter                             # total number of iterations used
    ))
  }
  # example output final:
  # Row 1: Year 11, P(SB > SB_hist) = 0.6 (calculated from 6 iterations)
  # Row 2: Year 12, P(SB > SB_hist) = 0.8 (calculated from 6 iterations)
  # Row 3: Year 13, P(SB > SB_hist) = 0.5 (calculated from 6 iterations)


  #============================================================================
  # METRIC 2: P(Catch_year / Catch_final_historical > 1.0)
  #============================================================================
  # for each projection year, we calculate the probability (across iterations)
  # that total catch exceeds the final historical year


  # get total catch (aggregate across all areas and fleets)
  if(is_multifleet) {
    # sum across areas and fleets: [year, iter, area, fleet] -> [year, iter]
    catchB_total <- apply(sim_result$dynamics$multifleet$catchB_by_fleet,
                          c(1, 2), sum)
  } else {
    # sum across areas: [year, iter, area] -> [year, iter]
    catchB_total <- apply(sim_result$dynamics$catchB, c(1, 2), sum)
  }

  # reference: last historical year catch for each iteration
  catch_hist_final <- catchB_total[hist_end, ]  # vector of length n_iter

  # calculate ratio for each projection year
  catch_metrics <- data.frame()

  for(proj_year in proj_start:proj_end) {
    real_year <- proj_year - 1

    # ratio for each iteration
    catch_ratio <- catchB_total[proj_year, ] / catch_hist_final

    # calculate probability
    # mean (proportion of true values)
    prob_above_1 <- mean(catch_ratio > 1.0) #catch_ratio > 1.0 creates logical vector [TRUE, FALSE, TRUE]

    mean_ratio_val <- mean(catch_ratio)      # mean of [1.07, 0.96, 1.11, 1.05, 0.95, 1.06]
    median_ratio_val <- median(catch_ratio)  # median of the vector
    sd_ratio_val <- sd(catch_ratio)          # standard deviation
    min_ratio_val <- min(catch_ratio)        # minimum value = 0.95
    max_ratio_val <- max(catch_ratio)        # maximum value = 1.11


    # store results
    catch_metrics <- rbind(catch_metrics, data.frame(
      mp_name = mp_name,
      year = real_year,
      proj_year_index = proj_year - proj_start + 1,
      metric = "P(Catch > Catch_hist)",
      probability = prob_above_1,                    #probability across iterations
      mean_ratio = mean_ratio_val,
      median_ratio = median_ratio_val,
      sd_ratio = sd_ratio_val,
      min_ratio = min_ratio_val,
      max_ratio = max_ratio_val,
      n_iterations = n_iter
    ))
  }

  #============================================================================
  # METRIC 3: Interannual Variability (IAV) - P(|change| < 15%)
  #============================================================================

  # for each projection year, calculate the probability (across iterations)
  # that the year-to-year catch change is less than 15%

  #example 6 iter
  # year 11 catch: [150, 130, 160, 145, 135, 148]
  # year 12 catch: [165, 125, 170, 150, 140, 145]
  # percent of change: |[(165-150)/150*100, (125-130)/130*100, ...]|
  #               = |[10%, -3.8%, 6.25%, 3.4%, 3.7%, -2.0%]|
  #               = [10%, 3.8%, 6.25%, 3.4%, 3.7%, 2.0%] (abs values)
  # changes < 15%: [TRUE, TRUE, TRUE, TRUE, TRUE, TRUE] = 6 out of 6
  #  P(IAV < 15%) = 6/6 = 1.0 = 100%

  #calculate year-to-year percent change in catch
  iav_metrics <- data.frame()
  #need 2 consecutive years to calculate change
  #proj_start is the first projection year
  #proj_start + 1 is the second projection year

  for(proj_year in (proj_start + 1):proj_end) {  # start from 2nd projection year
    real_year <- proj_year - 1                  #convert j-index to real year
    prev_year <- proj_year - 1                  #previous year in j-index

    # Calculate percent change for each iteration
    # eg:
    # catchB_total[proj_year, ] extracts one row
    # catch_current is a vector of length n_iter

    catch_current <- catchB_total[proj_year, ]

    #get catch for previous year (all iterations)
    #catchB_total[prev_year, ] extracts one row
    #catch_previous is a vector of length n_iter
    catch_previous <- catchB_total[prev_year, ]

    # calc percent change: (current - previous) / previous * 100 (abs values)
    # vector of length n_iter (one value per iteration)
    pct_change <- abs((catch_current - catch_previous) / catch_previous * 100)

    # calculate prob that change is less than 15%
    # mean of TRUE values ( 5TRUE out 6 total)
    prob_stable <- mean(pct_change < 15) #logical

    mean_pct <- mean(pct_change)        # mean of change  eg, [9.8, 3.8, 6.2, 3.3, 3.4, 1.7] = 4.7%
    median_pct <- median(pct_change)    # median of the vector = 3.6%
    sd_pct <- sd(pct_change)            # standard deviation = 3.0%
    min_pct <- min(pct_change)          # minimum = 1.7%
    max_pct <- max(pct_change)          # maximum = 9.8%


    #store results
    iav_metrics <- rbind(iav_metrics, data.frame(
      mp_name = mp_name,
      year = real_year,                               #real year number
      proj_year_index = proj_year - proj_start + 1,   #projection year index
      metric = "P(IAV < 15%)",
      probability = prob_stable,                     #probability across iterations
      mean_pct_change = mean_pct,                    #mean % change across iterations
      median_pct_change = median_pct,                #median % change across iterations
      sd_pct_change = sd_pct,                        #SD across iterations
      min_pct_change = min_pct,                      #min across iterations
      max_pct_change = max_pct,                      #max across iterations
      n_iterations = n_iter                          #total iterations used
    ))
  }

  #example of the results:
  # Row 1: Year 13, P(IAV < 15%) = 1.00 (100% of iterations changed < 15%)
  # Row 2: Year 14, P(IAV < 15%) = 0.83 (83% of iterations changed < 15%)
  # Row 3: Year 15, P(IAV < 15%) = 0.67 (67% of iterations changed < 15%)

  # Return all metrics
  return(list(
    sb_metrics = sb_metrics,
    catch_metrics = catch_metrics,
    iav_metrics = iav_metrics,
    summary = data.frame(
      mp_name = mp_name,
      n_hist = n_hist,
      n_proj = n_proj,
      n_areas = n_areas,
      n_iter = n_iter,
      is_multifleet = is_multifleet,
      n_fleets = if(is_multifleet) n_fleets else 1
    )
  ))
}

#Now for multiple MPs
calculate_multi_mp_metrics <- function(mp_list) {

  all_sb_metrics <- data.frame()
  all_catch_metrics <- data.frame()
  all_iav_metrics <- data.frame()
  all_summaries <- data.frame()

  for(mp_name in names(mp_list)) {

    metrics <- calculate_performance_metrics(mp_list[[mp_name]], mp_name)

    all_sb_metrics <- rbind(all_sb_metrics, metrics$sb_metrics)
    all_catch_metrics <- rbind(all_catch_metrics, metrics$catch_metrics)
    all_iav_metrics <- rbind(all_iav_metrics, metrics$iav_metrics)
    all_summaries <- rbind(all_summaries, metrics$summary)
  }

  cat("\n======================\n")
  cat("COMPLETED: Processed", length(mp_list), "management procedures\n")
  cat("\n======================\n")

  return(list(
    sb_metrics = all_sb_metrics,
    catch_metrics = all_catch_metrics,
    iav_metrics = all_iav_metrics,
    summaries = all_summaries
  ))
}


plot_performance_metrics <- function(metrics_combined) {

  # colorfor MPs
  n_mps <- length(unique(c(metrics_combined$sb_metrics$mp_name)))
  mp_colors <- scales::hue_pal()(n_mps)

  #short-term: years 1-5, Medium-term: years 6-10, Long-term: years 10+
  short_medium_boundary <- 5   # After year 5
  medium_long_boundary <- 10   # After year 10

  #==========================================================================
  # SB Performance
  #==========================================================================

  p_sb <- ggplot(metrics_combined$sb_metrics,
                 aes(x = proj_year_index, y = probability,
                     color = mp_name, shape = mp_name)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = short_medium_boundary, linetype = "dotted",
               color = "blue", alpha = 1) +
    geom_vline(xintercept = medium_long_boundary, linetype = "dotted",
               color = "blue", alpha = 1) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = seq(1, max(metrics_combined$sb_metrics$proj_year_index))) +
    labs(
      title = "Spawning Biomass Performance",
      subtitle = "Probability SB exceeds final historical year",
      x = "Projection Year",
      y = "P(SB > SB_hist)",
      color = "Management Procedure",
      shape = "Management Procedure"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray70")
    )

  #==========================================================================
  # Mean catch Performance
  #==========================================================================

  p_catch <- ggplot(metrics_combined$catch_metrics,
                    aes(x = proj_year_index, y = probability,
                        color = mp_name, shape = mp_name)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = short_medium_boundary, linetype = "dotted",
               color = "blue", alpha = 1) +
    geom_vline(xintercept = medium_long_boundary, linetype = "dotted",
               color = "blue", alpha = 1) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = seq(1, max(metrics_combined$catch_metrics$proj_year_index))) +
    labs(
      title = "Catch Performance",
      subtitle = "Probability catch exceeds final historical year",
      x = "Projection Year",
      y = "P(Catch > Catch_hist)",
      color = "Management Procedure",
      shape = "Management Procedure"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray70")
    )

  #==========================================================================
  # Interannual Variability (IAV)
  #==========================================================================

  p_iav <- ggplot(metrics_combined$iav_metrics,
                  aes(x = proj_year_index, y = probability,
                      color = mp_name, shape = mp_name)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "darkgreen",
               alpha = 0.5) +
    geom_vline(xintercept = short_medium_boundary, linetype = "dotted",
               color = "blue", alpha = 1) +
    geom_vline(xintercept = medium_long_boundary, linetype = "dotted",
               color = "blue", alpha = 1) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = seq(1, max(metrics_combined$iav_metrics$proj_year_index))) +
    labs(
      title = "Catch Stability (Interannual Variability)",
      subtitle = "Probability year-to-year catch change < 15%",
      x = "Projection Year",
      y = "P(IAV < 15%)",
      color = "Management Procedure",
      shape = "Management Procedure"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray70")
    )

  #==========================================================================
  # all three metrics in one figure
  #==========================================================================

  p_combined <- grid.arrange(p_sb, p_catch, p_iav, ncol = 1)

  return(list(
    p_sb = p_sb,
    p_catch = p_catch,
    p_iav = p_iav,
    p_combined = p_combined
  ))
}


create_summary_tables <- function(metrics_combined) {

  #==========================================================================
  # TABLE 1: Performance Summary
  #==========================================================================

  # average probability across all projection years
  sb_summary <- metrics_combined$sb_metrics %>%
    group_by(mp_name) %>%
    summarize(
      avg_prob_sb = mean(probability),
      avg_mean_ratio = mean(mean_ratio),
      .groups = "drop"
    )

  catch_summary <- metrics_combined$catch_metrics %>%
    group_by(mp_name) %>%
    summarize(
      avg_prob_catch = mean(probability),
      avg_mean_ratio = mean(mean_ratio),
      .groups = "drop"
    )

  iav_summary <- metrics_combined$iav_metrics %>%
    group_by(mp_name) %>%
    summarize(
      avg_prob_stable = mean(probability),
      avg_pct_change = mean(mean_pct_change),
      .groups = "drop"
    )

  # combine into overall summary
  overall_summary <- sb_summary %>%
    left_join(catch_summary, by = "mp_name") %>%
    left_join(iav_summary, by = "mp_name") %>%
    left_join(metrics_combined$summaries, by = "mp_name")

  #==========================================================================
  # TABLE 2: year to year comparison
  #==========================================================================

  # Combine all metrics into long format
  year_comparison <- bind_rows(
    metrics_combined$sb_metrics %>%
      select(mp_name, proj_year_index, metric, probability),
    metrics_combined$catch_metrics %>%
      select(mp_name, proj_year_index, metric, probability),
    metrics_combined$iav_metrics %>%
      select(mp_name, proj_year_index, metric, probability)
  ) %>%
    pivot_wider(names_from = metric, values_from = probability)

  #==========================================================================
  # TABLE 3: Final year performance
  #==========================================================================

  final_year_sb <- metrics_combined$sb_metrics %>%
    group_by(mp_name) %>%
    filter(proj_year_index == max(proj_year_index)) %>%
    select(mp_name, prob_sb_final = probability,
           mean_ratio_sb_final = mean_ratio)

  final_year_catch <- metrics_combined$catch_metrics %>%
    group_by(mp_name) %>%
    filter(proj_year_index == max(proj_year_index)) %>%
    select(mp_name, prob_catch_final = probability,
           mean_ratio_catch_final = mean_ratio)

  final_year_iav <- metrics_combined$iav_metrics %>%
    group_by(mp_name) %>%
    filter(proj_year_index == max(proj_year_index)) %>%
    select(mp_name, prob_stable_final = probability,
           avg_change_final = mean_pct_change)

  final_year_summary <- final_year_sb %>%
    left_join(final_year_catch, by = "mp_name") %>%
    left_join(final_year_iav, by = "mp_name")

  return(list(
    overall = overall_summary,
    year_by_year = year_comparison,
    final_year = final_year_summary
  ))
}



save_performance_outputs <- function(metrics_combined,
                                     output_dir = getwd(),
                                     file_prefix = "mp_performance") {

  # Create plots
  plots <- plot_performance_metrics(metrics_combined)

  # Create tables
  tables <- create_summary_tables(metrics_combined)

  # Save plots
  ggsave(filename = file.path(output_dir, paste0(file_prefix, "_sb.jpeg")),
         plot = plots$p_sb, width = 10, height = 6, dpi = 300)

  ggsave(filename = file.path(output_dir, paste0(file_prefix, "_catch.jpeg")),
         plot = plots$p_catch, width = 10, height = 6, dpi = 300)

  ggsave(filename = file.path(output_dir, paste0(file_prefix, "_iav.jpeg")),
         plot = plots$p_iav, width = 10, height = 6, dpi = 300)

  ggsave(filename = file.path(output_dir, paste0(file_prefix, "_combined.jpeg")),
         plot = plots$p_combined, width = 10, height = 16, dpi = 300)

  # Save tables as CSV
  write.csv(tables$overall,
            file = file.path(output_dir, paste0(file_prefix, "_summary.csv")),
            row.names = FALSE)

  write.csv(tables$year_by_year,
            file = file.path(output_dir, paste0(file_prefix, "_year_comparison.csv")),
            row.names = FALSE)

  write.csv(tables$final_year,
            file = file.path(output_dir, paste0(file_prefix, "_final_year.csv")),
            row.names = FALSE)

  #detailed metrics
  write.csv(metrics_combined$sb_metrics,
            file = file.path(output_dir, paste0(file_prefix, "_sb_detailed.csv")),
            row.names = FALSE)

  write.csv(metrics_combined$catch_metrics,
            file = file.path(output_dir, paste0(file_prefix, "_catch_detailed.csv")),
            row.names = FALSE)

  write.csv(metrics_combined$iav_metrics,
            file = file.path(output_dir, paste0(file_prefix, "_iav_detailed.csv")),
            row.names = FALSE)

  return(list(plots = plots, tables = tables))
}

# usage
result_mp1 <- readProjection(getwd(), "test_multi_indexratio_length_MP")
result_mp2 <- readProjection(getwd(), "test_multifleet_islope_length_MP")

# create named list of MPs to compare
mp_list <- list(
  "Indexratio MP" = result_mp1,
  "slopeMP" = result_mp2
)

# calculate all metrics
metrics <- calculate_multi_mp_metrics(mp_list)

# create and save all outputs
outputs <- save_performance_outputs(metrics,
                                   output_dir = getwd(),
                                   file_prefix = "mp_comparison")


#individual plot
print(outputs$plots$p_sb)
print(outputs$plots$p_catch)
print(outputs$plots$p_iav)

# view tables in console
print(outputs$tables$overall)
print(outputs$tables$final_year)



