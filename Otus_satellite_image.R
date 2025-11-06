# Otsu's Thresholding for Satellite Images

# Load required libraries
library(jpeg)
library(png)
library(httr)
library(imager)

# Function to download and process satellite image
process_satellite_image <- function(image_url) {
  # Create temporary file
  temp_file <- tempfile(fileext = ".jpg")
  
  # Try to download the image
  tryCatch({
    download.file(image_url, temp_file, mode = "wb", quiet = TRUE)
    cat("Image downloaded successfully\n")
    
    # Check file type and load image
    if (grepl("\\.jpg$|\\.jpeg$", image_url, ignore.case = TRUE)) {
      img <- readJPEG(temp_file)
    } else if (grepl("\\.png$", image_url, ignore.case = TRUE)) {
      img <- readPNG(temp_file)
    } else {
      # Try with imager for other formats
      img <- load.image(temp_file)
      img <- as.array(img)[,,1:3]  # Convert to array format
    }
    
    return(img)
  }, error = function(e) {
    cat("Error downloading image:", e$message, "\n")
    cat("Using built-in sample image instead\n")
    return(create_sample_image())
  })
}

# Create a sample image for demonstration if download fails
create_sample_image <- function() {
  # Create a synthetic satellite-like image
  set.seed(123)
  img <- array(0, dim = c(300, 300, 3))
  
  # Add different regions (water, vegetation, urban)
  # Water (dark blue)
  img[1:100, 1:100, 1] <- 0.1  # R
  img[1:100, 1:100, 2] <- 0.2  # G
  img[1:100, 1:100, 3] <- 0.4  # B
  
  # Vegetation (green)
  img[1:100, 201:300, 1] <- 0.2
  img[1:100, 201:300, 2] <- 0.5
  img[1:100, 201:300, 3] <- 0.1
  
  # Urban (gray)
  img[201:300, 1:100, 1] <- 0.6
  img[201:300, 1:100, 2] <- 0.6
  img[201:300, 1:100, 3] <- 0.6
  
  # Mixed areas
  img[101:200, 101:200, ] <- runif(100*100*3, 0.3, 0.7)
  
  return(img)
}

# Enhanced Otsu's threshold function for images
otsu_threshold_image <- function(image, channel = "grayscale") {
  
  # Convert image to appropriate channel
  if (length(dim(image)) == 3) {
    if (channel == "red") {
      data_vector <- as.vector(image[,,1])
    } else if (channel == "green") {
      data_vector <- as.vector(image[,,2])
    } else if (channel == "blue") {
      data_vector <- as.vector(image[,,3])
    } else {
      # Convert to grayscale using standard weights
      data_vector <- as.vector(0.299 * image[,,1] + 0.587 * image[,,2] + 0.114 * image[,,3])
    }
  } else {
    data_vector <- as.vector(image)
  }
  
  # Normalize to 0-255 for 8-bit image processing
  data_vector <- round(data_vector * 255)
  data_vector <- data_vector[!is.na(data_vector)]
  
  # Calculate histogram
  hist_data <- hist(data_vector, breaks=0:256, plot=FALSE)
  counts <- hist_data$counts
  bin_centers <- hist_data$mids
  
  # Remove zero counts for efficiency
  non_zero <- counts > 0
  counts <- counts[non_zero]
  bin_centers <- bin_centers[non_zero]
  
  # Normalize histogram to get probability distribution
  total_points <- sum(counts)
  probabilities <- counts / total_points
  
  # Initialize variables
  optimal_threshold <- 0
  max_variance <- 0
  
  # Pre-calculate cumulative sums for efficiency
  cumulative_sum <- cumsum(probabilities)
  cumulative_mean <- cumsum(bin_centers * probabilities)
  global_mean <- cumulative_mean[length(cumulative_mean)]
  
  # Store results
  results <- data.frame(
    threshold = numeric(),
    weight_bg = numeric(),
    mean_bg = numeric(),
    weight_fg = numeric(),
    mean_fg = numeric(),
    within_var = numeric(),
    between_var = numeric()
  )
  
  # Iterate through all possible thresholds
  for (threshold in 1:(length(probabilities)-1)) {
    
    # Calculate background statistics
    weight_bg <- cumulative_sum[threshold]
    mean_bg <- ifelse(weight_bg > 0, cumulative_mean[threshold] / weight_bg, 0)
    
    # Calculate foreground statistics
    weight_fg <- 1 - weight_bg
    mean_fg <- ifelse(weight_fg > 0, (global_mean - cumulative_mean[threshold]) / weight_fg, 0)
    
    # Calculate variances
    between_variance <- weight_bg * (mean_bg - global_mean)^2 + 
                       weight_fg * (mean_fg - global_mean)^2
    
    # Store results
    results <- rbind(results, data.frame(
      threshold = bin_centers[threshold],
      weight_bg = weight_bg,
      mean_bg = mean_bg,
      weight_fg = weight_fg,
      mean_fg = mean_fg,
      between_var = between_variance
    ))
    
    # Update optimal threshold
    if (between_variance > max_variance && weight_bg > 0 && weight_fg > 0) {
      max_variance <- between_variance
      optimal_threshold <- bin_centers[threshold]
    }
  }
  
  return(list(
    optimal_threshold = optimal_threshold / 255,  # Normalize back to 0-1
    all_results = results,
    global_mean = global_mean / 255,
    histogram = hist_data,
    data_vector = data_vector
  ))
}

# Function to apply threshold to image
apply_image_threshold <- function(image, threshold, channel = "grayscale") {
  if (length(dim(image)) == 3) {
    if (channel == "grayscale") {
      gray_image <- 0.299 * image[,,1] + 0.587 * image[,,2] + 0.114 * image[,,3]
      binary_image <- matrix(0, nrow = nrow(image), ncol = ncol(image))
      binary_image[gray_image >= threshold] <- 1
    } else {
      channel_idx <- switch(channel,
                           "red" = 1,
                           "green" = 2, 
                           "blue" = 3)
      binary_image <- matrix(0, nrow = nrow(image), ncol = ncol(image))
      binary_image[image[,,channel_idx] >= threshold] <- 1
    }
  } else {
    binary_image <- matrix(0, nrow = nrow(image), ncol = ncol(image))
    binary_image[image >= threshold] <- 1
  }
  return(binary_image)
}

# Main analysis function
analyze_satellite_image <- function(image_url) {
  cat("Processing satellite image...\n")
  
  # Download and load image
  img <- process_satellite_image(image_url)
  
  # Analyze different channels
  channels <- c("grayscale", "red", "green", "blue")
  results <- list()
  
  for (channel in channels) {
    cat("\n=== Analyzing", channel, "channel ===\n")
    otsu_result <- otsu_threshold_image(img, channel)
    results[[channel]] <- otsu_result
    
    cat("Optimal threshold:", round(otsu_result$optimal_threshold, 3), "\n")
    cat("Global mean:", round(otsu_result$global_mean, 3), "\n")
    
    # Apply threshold
    binary_img <- apply_image_threshold(img, otsu_result$optimal_threshold, channel)
    
    # Calculate coverage statistics
    foreground_pixels <- sum(binary_img)
    total_pixels <- length(binary_img)
    coverage <- foreground_pixels / total_pixels
    
    cat("Foreground coverage:", round(coverage * 100, 2), "%\n")
    cat("Foreground pixels:", foreground_pixels, "\n")
    cat("Background pixels:", total_pixels - foreground_pixels, "\n")
  }
  
  # Visualization
  par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))
  
  # Original image
  plot(as.raster(img), main = "Original Satellite Image")
  
  # Grayscale version
  gray_img <- 0.299 * img[,,1] + 0.587 * img[,,2] + 0.114 * img[,,3]
  plot(as.raster(gray_img), main = "Grayscale Image")
  
  # Binary results for each channel
  for (channel in channels) {
    binary_img <- apply_image_threshold(img, results[[channel]]$optimal_threshold, channel)
    plot(as.raster(binary_img), main = paste("Binary -", channel))
  }
  
  # Histogram with thresholds
  gray_vector <- as.vector(gray_img)
  hist(gray_vector, breaks = 50, main = "Intensity Histogram with Thresholds",
       xlab = "Intensity", col = "lightblue", border = "black")
  
  colors <- c("red", "green", "blue", "black")
  for (i in 1:length(channels)) {
    abline(v = results[[channels[i]]]$optimal_threshold, 
           col = colors[i], lwd = 2, lty = i)
  }
  legend("topright", legend = channels, col = colors, lty = 1:4, lwd = 2)
  
  # Detailed threshold comparison
  cat("\n=== THRESHOLD COMPARISON ACROSS CHANNELS ===\n")
  comparison_df <- data.frame(
    Channel = channels,
    Threshold = sapply(channels, function(ch) round(results[[ch]]$optimal_threshold, 3)),
    Global_Mean = sapply(channels, function(ch) round(results[[ch]]$global_mean, 3)),
    Coverage_Percent = sapply(channels, function(ch) {
      binary_img <- apply_image_threshold(img, results[[ch]]$optimal_threshold, ch)
      round(sum(binary_img) / length(binary_img) * 100, 2)
    })
  )
  print(comparison_df)
  
  return(list(
    image = img,
    results = results,
    comparison = comparison_df
  ))
}

# Apply to the satellite image URL
satellite_url <- "https://wp-cdn.apollomapping.com/web_assets/user_uploads/2021/10/19084817/DoleMazeGarden_11_28_2018_WV4_30cmcolor_ENHANCE.jpg"

# Run the analysis
final_results <- analyze_satellite_image(satellite_url)

# Additional analysis for image segmentation
cat("\n=== IMAGE SEGMENTATION ANALYSIS ===\n")

# Use the best channel (usually grayscale for overall segmentation)
best_channel <- "grayscale"
best_threshold <- final_results$results[[best_channel]]$optimal_threshold
binary_segmentation <- apply_image_threshold(final_results$image, best_threshold, best_channel)

# Calculate region properties
cat("Segmentation using", best_channel, "channel:\n")
cat("Optimal threshold:", round(best_threshold, 3), "\n")

# Connected components analysis (simple version)
find_connected_components <- function(binary_image) {
  labeled <- matrix(0, nrow = nrow(binary_image), ncol = ncol(binary_image))
  current_label <- 1
  
  for (i in 1:nrow(binary_image)) {
    for (j in 1:ncol(binary_image)) {
      if (binary_image[i, j] == 1 && labeled[i, j] == 0) {
        # Simple region growing (for demonstration)
        stack <- list(c(i, j))
        region_size <- 0
        
        while (length(stack) > 0) {
          pixel <- stack[[1]]
          stack <- stack[-1]
          row <- pixel[1]
          col <- pixel[2]
          
          if (row >= 1 && row <= nrow(binary_image) && 
              col >= 1 && col <= ncol(binary_image) &&
              binary_image[row, col] == 1 && labeled[row, col] == 0) {
            
            labeled[row, col] <- current_label
            region_size <- region_size + 1
            
            # Add neighbors to stack
            neighbors <- list(
              c(row-1, col), c(row+1, col),
              c(row, col-1), c(row, col+1)
            )
            stack <- c(stack, neighbors)
          }
        }
        
        cat("Region", current_label, "size:", region_size, "pixels\n")
        current_label <- current_label + 1
      }
    }
  }
  
  return(labeled)
}

# Simple connected components (comment out for large images as it can be slow)
# cat("\nConnected components analysis:\n")
# components <- find_connected_components(binary_segmentation)
# cat("Total regions found:", max(components), "\n")

# Show final visualization
par(mfrow = c(1, 2), mar = c(2, 2, 2, 1))
plot(as.raster(final_results$image), main = "Original Satellite Image")
plot(as.raster(binary_segmentation), main = "Otsu Segmentation Result")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Otsu's method has successfully segmented the satellite image\n")
cat("The optimal threshold separates foreground objects from background\n")
cat("This can be used for feature extraction, land cover classification, etc.\n")