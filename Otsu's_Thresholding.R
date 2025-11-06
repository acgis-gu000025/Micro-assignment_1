# Otsu's Thresholding Method for Ontario CD Population Data

# Ontario CD Data
ontario_data <- data.frame(
  Geos = c(3557000, 3529000, 3541000, 3536000),
  CD = c("Ontario - Algoma District", "Ontario - Brant County", 
         "Ontario - Bruce County", "Ontario - Chatham-Kent Division"),
  AvInc2000 = c(26153, 28705, 28097, 28757),
  AvIncMale = c(33108, 35750, 35435, 35687),
  AvIncFem = c(19451, 21947, 20692, 22060),
  Pop2001 = c(118567, 118485, 63892, 107709),
  AreaSqKm = c(48737.2, 1092.9, 4155.5, 2470.8)
)

# Function to calculate Otsu's threshold
otsu_threshold <- function(data_vector) {
  # Calculate histogram with appropriate breaks for population data
  hist_data <- hist(data_vector, breaks=length(data_vector), plot=FALSE)
  counts <- hist_data$counts
  bin_centers <- hist_data$mids
  
  # Remove zero counts for efficiency
  non_zero <- counts > 0
  counts <- counts[non_zero]
  bin_centers <- bin_centers[non_zero]
  
  # Normalize histogram to get probability distribution
  total_points <- sum(counts)
  probabilities <- counts / total_points
  
  # Initialize variables for optimal threshold search
  optimal_threshold <- 0
  max_variance <- 0
  best_within_variance <- Inf
  
  # Pre-calculate cumulative sums for efficiency
  cumulative_sum <- cumsum(probabilities)
  cumulative_mean <- cumsum(bin_centers * probabilities)
  global_mean <- cumulative_mean[length(cumulative_mean)]
  
  # Store results for all thresholds
  results <- data.frame(
    threshold = numeric(),
    weight_bg = numeric(),
    mean_bg = numeric(),
    var_bg = numeric(),
    weight_fg = numeric(),
    mean_fg = numeric(),
    var_fg = numeric(),
    within_var = numeric(),
    between_var = numeric()
  )
  
  # Iterate through all possible thresholds
  for (threshold in 1:(length(probabilities)-1)) {
    
    # b. Calculate weight, mean, variance for background
    weight_bg <- cumulative_sum[threshold]
    if (weight_bg > 0) {
      mean_bg <- cumulative_mean[threshold] / weight_bg
    } else {
      mean_bg <- 0
    }
    
    # Variance calculation for background
    var_bg <- 0
    if (weight_bg > 0) {
      for (i in 1:threshold) {
        var_bg <- var_bg + probabilities[i] * (bin_centers[i] - mean_bg)^2
      }
      var_bg <- var_bg / weight_bg
    }
    
    # c. Calculate weight, mean, variance for foreground
    weight_fg <- 1 - weight_bg
    if (weight_fg > 0) {
      mean_fg <- (global_mean - cumulative_mean[threshold]) / weight_fg
    } else {
      mean_fg <- 0
    }
    
    # Variance calculation for foreground
    var_fg <- 0
    if (weight_fg > 0) {
      for (i in (threshold+1):length(probabilities)) {
        var_fg <- var_fg + probabilities[i] * (bin_centers[i] - mean_fg)^2
      }
      var_fg <- var_fg / weight_fg
    }
    
    # d. Calculate within-class variance
    within_variance <- weight_bg * var_bg + weight_fg * var_fg
    
    # e. Alternative: Calculate between-class variance (faster approach)
    between_variance <- weight_bg * (mean_bg - global_mean)^2 + 
      weight_fg * (mean_fg - global_mean)^2
    
    # Store results
    results <- rbind(results, data.frame(
      threshold = bin_centers[threshold],
      weight_bg = weight_bg,
      mean_bg = mean_bg,
      var_bg = var_bg,
      weight_fg = weight_fg,
      mean_fg = mean_fg,
      var_fg = var_fg,
      within_var = within_variance,
      between_var = between_variance
    ))
    
    # Update optimal threshold (using between-class variance)
    if (between_variance > max_variance && weight_bg > 0 && weight_fg > 0) {
      max_variance <- between_variance
      optimal_threshold <- bin_centers[threshold]
    }
  }
  
  return(list(
    optimal_threshold = optimal_threshold,
    all_results = results,
    global_mean = global_mean,
    data_vector = data_vector
  ))
}

# Apply Otsu's method to Population data
population_data <- ontario_data$Pop2001
otsu_result <- otsu_threshold(population_data)

# Display detailed results
cat("=== ONTARIO CD POPULATION DATA ANALYSIS ===\n")
cat("Population Data:", population_data, "\n")
cat("Global mean population:", otsu_result$global_mean, "\n")
cat("Optimal threshold:", otsu_result$optimal_threshold, "\n\n")

# Find the optimal threshold details
optimal_idx <- which.min(otsu_result$all_results$within_var)
optimal_details <- otsu_result$all_results[optimal_idx, ]

cat("=== DETAILED RESULTS FOR OPTIMAL THRESHOLD ===\n")
cat("Threshold:", round(optimal_details$threshold, 2), "\n")
cat("Background - Weight:", round(optimal_details$weight_bg, 4), 
    "Mean:", round(optimal_details$mean_bg, 2), 
    "Variance:", round(optimal_details$var_bg, 2), "\n")
cat("Foreground - Weight:", round(optimal_details$weight_fg, 4), 
    "Mean:", round(optimal_details$mean_fg, 2), 
    "Variance:", round(optimal_details$var_fg, 2), "\n")
cat("Within-class variance:", round(optimal_details$within_var, 2), "\n")
cat("Between-class variance:", round(optimal_details$between_var, 2), "\n\n")

# Apply the threshold and classify the CDs
cat("=== CLASSIFICATION RESULTS ===\n")
for (i in 1:nrow(ontario_data)) {
  classification <- ifelse(population_data[i] >= otsu_result$optimal_threshold, 
                           "High Population", "Low Population")
  cat(sprintf("%-30s: %8d people -> %s\n", 
              ontario_data$CD[i], population_data[i], classification))
}

# Visualization
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# Histogram with threshold
hist(population_data, breaks = 10, main = "Population Distribution with Otsu Threshold",
     xlab = "Population", col = "lightblue", border = "black")
abline(v = otsu_result$optimal_threshold, col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Threshold =", round(otsu_result$optimal_threshold, 2)), 
       col = "red", lty = 2)

# Within-class variance
if (nrow(otsu_result$all_results) > 1) {
  plot(otsu_result$all_results$threshold, otsu_result$all_results$within_var, 
       type = "l", main = "Within-class Variance vs Threshold",
       xlab = "Threshold", ylab = "Within-class Variance", col = "blue")
  points(optimal_details$threshold, optimal_details$within_var, 
         col = "red", pch = 19, cex = 1.5)
}

# Between-class variance
if (nrow(otsu_result$all_results) > 1) {
  plot(otsu_result$all_results$threshold, otsu_result$all_results$between_var, 
       type = "l", main = "Between-class Variance vs Threshold",
       xlab = "Threshold", ylab = "Between-class Variance", col = "green")
  points(optimal_details$threshold, optimal_details$between_var, 
         col = "red", pch = 19, cex = 1.5)
}

# Bar plot of populations with classification
colors <- ifelse(population_data >= otsu_result$optimal_threshold, "red", "blue")
barplot(population_data, names.arg = abbreviate(ontario_data$CD, 15),
        main = "Population by CD with Classification",
        ylab = "Population", col = colors, las = 2)
abline(h = otsu_result$optimal_threshold, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("High Population", "Low Population"), 
       fill = c("red", "blue"))

# Statistical summary
cat("\n=== STATISTICAL SUMMARY ===\n")
cat("Minimum population:", min(population_data), "\n")
cat("Maximum population:", max(population_data), "\n")
cat("Mean population:", mean(population_data), "\n")
cat("Standard deviation:", sd(population_data), "\n")
cat("Optimal threshold separates data into:\n")
low_pop <- population_data[population_data < otsu_result$optimal_threshold]
high_pop <- population_data[population_data >= otsu_result$optimal_threshold]
cat("  Low population group:", length(low_pop), "CDs\n")
cat("  High population group:", length(high_pop), "CDs\n")