# Otsu's Thresholding for Grayscale Satellite Images
# This script implements Otsu's automatic thresholding method for image segmentation

# Load required libraries
library(jpeg)

# Function to download and process grayscale satellite image
process_satellite_image <- function(image_url) {
  temp_file <- tempfile(fileext = ".jpg")
  download.file(image_url, temp_file, mode = "wb", quiet = TRUE)
  img_gray <- readJPEG(temp_file)
  return(img_gray)
}

# Create a sample grayscale image for demonstration if download fails
create_sample_image <- function() {
  # Create a synthetic grayscale satellite-like image
  set.seed(123)
  img <- matrix(0, nrow = 300, ncol = 300)
  
  # Add different regions with varying grayscale intensities
  img[1:100, 1:100] <- 0.2      # Dark region (water)
  img[1:100, 201:300] <- 0.5    # Medium region (vegetation)
  img[201:300, 1:100] <- 0.7    # Light region (urban)
  img[101:200, 101:200] <- runif(100*100, 0.3, 0.8)  # Mixed areas
  
  return(img)
}

# Enhanced Otsu's threshold function for grayscale images
# Otsu's method automatically determines the optimal threshold value
otsu_threshold_image <- function(image) {
  
  # ========== STEP a: Calculate Histogram ==========
  # Convert grayscale image to vector
  data_vector <- as.vector(image)
  
  # Normalize to 0-255 range for 8-bit image processing
  data_vector <- round(data_vector * 255)
  data_vector <- data_vector[!is.na(data_vector)]
  
  # a. Apply a threshold (histogram can be used in determining this)
  # Create histogram with 256 bins (one for each possible intensity value)
  hist_data <- hist(data_vector, breaks = 0:256, plot = FALSE)
  counts <- hist_data$counts
  probabilities <- counts / sum(counts)
  
  cat("=== Otsu's Thresholding Algorithm ===\n")
  cat("a. Histogram calculated with", length(counts), "bins\n\n")
  
  # ========== STEP b, c, d: Calculate Statistics for All Thresholds ==========
  optimal_threshold <- 0
  max_between_variance <- 0
  
  # Pre-calculate cumulative sums for efficiency
  cumsum_prob <- cumsum(probabilities)
  cumsum_weighted <- cumsum((0:255) * probabilities)
  global_mean <- cumsum_weighted[256]
  
  # Create results dataframe to store all threshold statistics
  results_list <- list()
  
  # e. Repeat the calculation for all possible thresholds
  # (Using between-class-variance as it's computationally faster)
  cat("e. Calculating between-class-variance for all possible thresholds:\n")
  cat("   (Scanning", 255, "threshold values...)\n\n")
  
  for (t in 1:255) {
    # b. Calculate weight, mean, variance for background (pixels <= threshold)
    weight_bg <- cumsum_prob[t]
    if (weight_bg > 0) {
      mean_bg <- cumsum_weighted[t] / weight_bg
    } else {
      mean_bg <- 0
    }
    
    # c. Calculate weight, mean, variance for foreground (pixels > threshold)
    weight_fg <- 1 - weight_bg
    if (weight_fg > 0) {
      mean_fg <- (global_mean - cumsum_weighted[t]) / weight_fg
    } else {
      mean_fg <- 0
    }
    
    # d. Calculate between-class-variance
    # (More efficient than within-class-variance for finding optimal threshold)
    between_variance <- weight_bg * weight_fg * (mean_bg - mean_fg)^2
    
    # Store results
    results_list[[t]] <- data.frame(
      threshold = t,
      weight_bg = weight_bg,
      weight_fg = weight_fg,
      mean_bg = mean_bg,
      mean_fg = mean_fg,
      between_variance = between_variance
    )
    
    # f. Specify the best threshold (maximum between-class-variance)
    if (between_variance > max_between_variance && weight_bg > 0 && weight_fg > 0) {
      max_between_variance <- between_variance
      optimal_threshold <- t
    }
  }
  
  # Convert results list to dataframe
  results_df <- do.call(rbind, results_list)
  
  cat("b-c. Weight and mean calculated for background and foreground\n")
  cat("d.   Within-class-variance used in threshold evaluation\n")
  cat("f. OPTIMAL THRESHOLD FOUND:", optimal_threshold, "\n")
  cat("   Maximum between-class-variance:", round(max_between_variance, 6), "\n\n")
  
  return(list(
    optimal_threshold = optimal_threshold / 255,  # Normalize back to 0-1
    threshold_int = optimal_threshold,
    all_results = results_df,
    histogram = hist_data,
    max_variance = max_between_variance
  ))
}

# Function to apply threshold to grayscale image
apply_image_threshold <- function(image, threshold) {
  # Create binary image where pixels >= threshold are foreground (1)
  # and pixels < threshold are background (0)
  binary_image <- matrix(0, nrow = nrow(image), ncol = ncol(image))
  binary_image[image >= threshold] <- 1
  return(binary_image)
}

# Main analysis function
analyze_satellite_image <- function(image_url) {
  cat("Processing GRAYSCALE satellite image...\n")
  cat("=" %*% rep("=", 40), "\n\n")
  
  # Download and load grayscale image
  img_gray <- process_satellite_image(image_url)
  
  # Perform Otsu's thresholding analysis
  otsu_result <- otsu_threshold_image(img_gray)
  
  # Display results
  cat("Threshold in 0-1 range:", round(otsu_result$optimal_threshold, 4), "\n")
  cat("Threshold in 0-255 range:", otsu_result$threshold_int, "\n\n")
  
  # Apply threshold
  binary_img <- apply_image_threshold(img_gray, otsu_result$optimal_threshold)
  
  # Calculate coverage statistics
  foreground_pixels <- sum(binary_img)
  total_pixels <- length(binary_img)
  coverage <- foreground_pixels / total_pixels
  
  cat("=== SEGMENTATION RESULTS ===\n")
  cat("Foreground coverage:", round(coverage * 100, 2), "%\n")
  cat("Foreground pixels:", foreground_pixels, "\n")
  cat("Background pixels:", total_pixels - foreground_pixels, "\n\n")
  
  # Visualization
  par(mfrow = c(2, 2), mar = c(2, 2, 2, 1))
  
  # Original grayscale image
  plot(as.raster(img_gray), main = "Original Grayscale Image")
  
  # Histogram with optimal threshold
  hist(as.vector(round(img_gray * 255)), breaks = 50, 
       main = "Intensity Histogram",
       xlab = "Intensity (0-255)", col = "lightblue", border = "black")
  abline(v = otsu_result$threshold_int, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("Optimal threshold =", otsu_result$threshold_int), 
         col = "red", lty = 2, lwd = 2)
  
  # Binary segmentation result
  plot(as.raster(binary_img), main = "Binary Segmentation Result")
  
  # Between-class-variance plot
  plot(otsu_result$all_results$threshold, otsu_result$all_results$between_variance,
       type = "l", main = "Between-class-variance vs Threshold",
       xlab = "Threshold Value", ylab = "Between-class-variance", 
       col = "darkblue", lwd = 2)
  points(otsu_result$threshold_int, otsu_result$max_variance, 
         col = "red", pch = 16, cex = 2)
  legend("topleft", legend = "Optimal threshold", col = "red", pch = 16, cex = 1.2)
  
  return(list(
    image = img_gray,
    threshold = otsu_result$optimal_threshold,
    binary = binary_img,
    results = otsu_result
  ))
}

# Apply to the satellite image URL
satellite_url <- "https://www.beneaththewaves.net/Photography/Images/01-Greyscale-Basic-02-OnEarth_Example-GS-SRTM.JPG"

# Run the analysis
final_results <- analyze_satellite_image(satellite_url)

cat("\n========================================\n")
cat("=== ANALYSIS COMPLETE ===\n")
cat("========================================\n")
cat("Otsu's method has successfully segmented the grayscale satellite image\n")
cat("The optimal threshold separates foreground objects from background\n")