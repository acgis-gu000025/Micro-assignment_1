# Otsu's Thresholding Method for Ontario CD Population Data

# Ontario CD Data - Complete dataset
ontario_data <- data.frame(
  Geos = c(3557000, 3529000, 3541000, 3536000, 3556000, 3522000, 3518000, 3534000, 
           3537000, 3510000, 3553000, 3542000, 3528000, 3546000, 3524000, 3525000, 
           3512000, 3540000, 3516000, 3560000, 3538000, 3509000, 3507000, 3511000, 
           3551000, 3539000, 3544000, 3526000, 3548000, 3514000, 3506000, 3532000, 
           3549000, 3521000, 3531000, 3515000, 3502000, 3513000, 3559000, 3547000, 
           3543000, 3501000, 3552000, 3558000, 3554000, 3520000, 3530000, 3523000, 3519000),
  CD = c("Ontario - Algoma District", "Ontario - Brant County", "Ontario - Bruce County", 
         "Ontario - Chatham-Kent Division", "Ontario - Cochrane District", "Ontario - Dufferin County", 
         "Ontario - Durham Regional Municipality", "Ontario - Elgin County", "Ontario - Essex County", 
         "Ontario - Frontenac County", "Ontario - Greater Sudbury Division", "Ontario - Grey County", 
         "Ontario - Haldimand-Norfolk Regional Municipality", "Ontario - Haliburton County", 
         "Ontario - Halton Regional Municipality", "Ontario - Hamilton Division", "Ontario - Hastings County", 
         "Ontario - Huron County", "Ontario - Kawartha Lakes Division", "Ontario - Kenora District", 
         "Ontario - Lambton County", "Ontario - Lanark County", "Ontario - Leeds and Grenville United Counties", 
         "Ontario - Lennox and Addington County", "Ontario - Manitoulin District", "Ontario - Middlesex County", 
         "Ontario - Muskoka District Municipality", "Ontario - Niagara Regional Municipality", 
         "Ontario - Nipissing District", "Ontario - Northumberland County", "Ontario - Ottawa Division", 
         "Ontario - Oxford County", "Ontario - Parry Sound District", "Ontario - Peel Regional Municipality", 
         "Ontario - Perth County", "Ontario - Peterborough County", "Ontario - Prescott and Russell United Counties", 
         "Ontario - Prince Edward Division", "Ontario - Rainy River District", "Ontario - Renfrew County", 
         "Ontario - Simcoe County", "Ontario - Stormont, Dundas and Glengarry United Counties", 
         "Ontario - Sudbury District", "Ontario - Thunder Bay District", "Ontario - Timiskaming District", 
         "Ontario - Toronto Division", "Ontario - Waterloo Regional Municipality", "Ontario - Wellington County", 
         "Ontario - York Regional Municipality"),
  AvInc2000 = c(26153, 28705, 28097, 28757, 27939, 33509, 34961, 28069, 33241, 30091, 
                29012, 27228, 27800, 22418, 43691, 29489, 26259, 27101, 27148, 27337, 
                30149, 29905, 29026, 26526, 21351, 31082, 26382, 29030, 26184, 28596, 
                38584, 29587, 24323, 34436, 29877, 27668, 29373, 27237, 28185, 26286, 
                30321, 26376, 26094, 29998, 24978, 34406, 32483, 32802, 38130),
  AvIncMale = c(33108, 35750, 35435, 35687, 35676, 43417, 44277, 35301, 42861, 36191, 
                36231, 32911, 34829, 26735, 57585, 36704, 32202, 33239, 33982, 33550, 
                39332, 36470, 35339, 32722, 24211, 38426, 31975, 36694, 32110, 35971, 
                47328, 37363, 29535, 42550, 36856, 34123, 35110, 32638, 35333, 32735, 
                37937, 32405, 33974, 37887, 30936, 42020, 41223, 41266, 47858),
  AvIncFem = c(19451, 21947, 20692, 22060, 19938, 23496, 25949, 20963, 23595, 24355, 
               22060, 21758, 20768, 17978, 30273, 22550, 20545, 21001, 20438, 21002, 
               21324, 23608, 22942, 20476, 18595, 24242, 20957, 21701, 20515, 21512, 
               30230, 22019, 19067, 26447, 23149, 21662, 23608, 22064, 20998, 19939, 
               22882, 20521, 17351, 22172, 19159, 27306, 23904, 24480, 28621),
  Pop2001 = c(118567, 118485, 63892, 107709, 85247, 51013, 506901, 81553, 374975, 
              138606, 155268, 89073, 104670, 15085, 375229, 490268, 125915, 59701, 
              69179, 61802, 126971, 62495, 96606, 39461, 12679, 403185, 53106, 410574, 
              82910, 77497, 774072, 99270, 39665, 988948, 73675, 125856, 76446, 24901, 
              22109, 95138, 377050, 109522, 22894, 150860, 34442, 2481494, 438515, 
              187313, 729254),
  AreaSqKm = c(48737.2, 1092.9, 4155.5, 2470.8, 141244, 1485.7, 2523.5, 1880.8, 
               1851.6, 3672.8, 3365, 4508.2, 2894.2, 4025.3, 967, 1117.1, 5978.3, 
               3407.6, 3059.2, 407167.3, 3001.7, 2979.4, 3350.6, 2777, 4759.6, 3317.1, 
               3890.4, 1863.1, 17064.6, 1903.2, 2778.6, 2039.4, 9222.3, 1242, 2218.4, 
               3806, 2001.2, 1050, 15473.7, 7403.6, 4840.5, 3307.2, 38350.5, 103714.4, 
               13280.1, 629.9, 1368.5, 2656.6, 1761.6)
)

# Function to calculate Otsu's threshold
otsu_threshold <- function(data_vector) {
  # Calculate histogram with appropriate breaks for population data
  hist_data <- hist(data_vector, breaks=50, plot=FALSE)
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
cat("Number of Census Divisions:", length(population_data), "\n")
cat("Population range:", min(population_data), "-", max(population_data), "\n")
cat("Global mean population:", round(otsu_result$global_mean, 2), "\n")
cat("Optimal threshold:", round(otsu_result$optimal_threshold, 2), "\n\n")

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
low_pop_count <- 0
high_pop_count <- 0

for (i in 1:nrow(ontario_data)) {
  classification <- ifelse(population_data[i] >= otsu_result$optimal_threshold, 
                           "High Population", "Low Population")
  if (classification == "Low Population") low_pop_count <- low_pop_count + 1
  if (classification == "High Population") high_pop_count <- high_pop_count + 1
  
  cat(sprintf("%-40s: %8d people -> %s\n", 
              ontario_data$CD[i], population_data[i], classification))
}

cat("\nClassification Summary:\n")
cat("Low population CDs:", low_pop_count, "\n")
cat("High population CDs:", high_pop_count, "\n")

# Visualization
par(mfrow = c(2, 2), mar = c(5, 4, 4, 2))

# Histogram with threshold
hist(population_data, breaks = 20, main = "Population Distribution with Otsu Threshold",
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
  legend("topleft", legend = "Optimal Threshold", col = "red", pch = 19)
}

# Between-class variance
if (nrow(otsu_result$all_results) > 1) {
  plot(otsu_result$all_results$threshold, otsu_result$all_results$between_var, 
       type = "l", main = "Between-class Variance vs Threshold",
       xlab = "Threshold", ylab = "Between-class Variance", col = "green")
  points(optimal_details$threshold, optimal_details$between_var, 
         col = "red", pch = 19, cex = 1.5)
  legend("topleft", legend = "Optimal Threshold", col = "red", pch = 19)
}

# Bar plot of populations with classification (top 20 for readability)
colors <- ifelse(population_data >= otsu_result$optimal_threshold, "red", "blue")
top_n <- min(20, length(population_data))
barplot(sort(population_data, decreasing = TRUE)[1:top_n], 
        names.arg = abbreviate(ontario_data$CD[order(population_data, decreasing = TRUE)][1:top_n], 20),
        main = "Top 20 Populations by CD with Classification",
        ylab = "Population", col = colors[order(population_data, decreasing = TRUE)][1:top_n], 
        las = 2, cex.names = 0.7)
abline(h = otsu_result$optimal_threshold, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("High Population", "Low Population"), 
       fill = c("red", "blue"), cex = 0.8)

# Statistical summary
cat("\n=== STATISTICAL SUMMARY ===\n")
cat("Minimum population:", min(population_data), "\n")
cat("Maximum population:", max(population_data), "\n")
cat("Mean population:", round(mean(population_data), 2), "\n")
cat("Median population:", median(population_data), "\n")
cat("Standard deviation:", round(sd(population_data), 2), "\n")

low_pop_data <- population_data[population_data < otsu_result$optimal_threshold]
high_pop_data <- population_data[population_data >= otsu_result$optimal_threshold]

cat("\nLow population group statistics:\n")
cat("  Count:", length(low_pop_data), "\n")
cat("  Range:", min(low_pop_data), "-", max(low_pop_data), "\n")
cat("  Mean:", round(mean(low_pop_data), 2), "\n")

cat("\nHigh population group statistics:\n")
cat("  Count:", length(high_pop_data), "\n")
cat("  Range:", min(high_pop_data), "-", max(high_pop_data), "\n")
cat("  Mean:", round(mean(high_pop_data), 2), "\n")

# Show the top 5 highest and lowest population CDs
cat("\n=== EXTREME POPULATION CDs ===\n")
cat("Top 5 Highest Population:\n")
high_pop_sorted <- ontario_data[order(-ontario_data$Pop2001), ][1:5, c("CD", "Pop2001")]
print(high_pop_sorted)

cat("\nTop 5 Lowest Population:\n")
low_pop_sorted <- ontario_data[order(ontario_data$Pop2001), ][1:5, c("CD", "Pop2001")]
print(low_pop_sorted)