#Resets R environment
rm(list = ls())

library(rstudioapi)
# Set working directory to the folder containing the current script
setwd(dirname(getActiveDocumentContext()$path))
# Check
getwd()

otsu_threshold <- function(x, bins = 256) {
  # Build histogram
  h <- hist(x, breaks = bins, plot = FALSE)
  
  counts <- h$counts
  mids <- h$mids
  total <- sum(counts)
  
  # Probabilities
  p <- counts / total
  
  # Cumulative sums
  w_b <- cumsum(p)                # background weight
  w_f <- 1 - w_b                  # foreground weight
  
  # Cumulative means
  mu_b <- cumsum(p * mids) / w_b
  mu_f <- (sum(p * mids) - cumsum(p * mids)) / w_f
  
  # Between-class variance
  sigma_b_sq <- w_b * w_f * (mu_b - mu_f)^2
  
  # Find max variance index
  t_index <- which.max(sigma_b_sq)
  
  # Corresponding threshold
  threshold <- mids[t_index]
  
  return(threshold)
}

install.packages("readxl") # Install if not already installed
library(readxl)
OntarioCDs <- read_excel("OntarioCDs.xlsx")

t_pop <- otsu_threshold(OntarioCDs$Pop2001)
t_pop



#Reads and assigns image to img, which is converted to grayscale as gray
library(png)
img <- readPNG("Grayscale-satellite-image-of-the-Fremont-River-basin.png")
dim(img)
gray <- 0.299*img[,,1] + 0.587*img[,,2] + 0.114*img[,,3]

#Plot the rasters to see the difference
plot(as.raster(img))
plot(as.raster(gray))

#Applies the otsu threshold on the grayscale image
thr <- otsu_threshold(gray)
thr

#Converts the grayscale image to monochrome based on otsu_threshold
binary <- gray > thr
plot(as.raster(binary))
