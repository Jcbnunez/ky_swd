#KY Weather slicing for Genomic samples!!
#Feb 28 2025

### libraries
library(tidyverse)
library(data.table)
library(gdata)
library(foreach)
library(nasapower)
library(sp)
library(lubridate)
library(magrittr)
library(readr)

weather.data <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/Ellie_MS/Weather_data/KYWeatherData.mgarvin.Jan25.csv")
weather.data %<>%
  mutate(date = as.Date(paste(YEAR,MO,DY, sep = "-"),
                        format = "%Y-%m-%d"))

meta <- fread("/netfiles/nunezlab/D_suzukii_resources/Datasets/KY_2020_2023/MASTER_SWD_METADATA-METADATA_v2_MG.csv")
meta %>%
  filter(province == "Kentucky") ->
  meta.ky

calc_stat <- function(variable, window, sample, statistic) {
  #variable=combo$var[1];window=combo$win[1];sample=combo$samp[1];statistic=combo$stat[1]
  
  #a lil sanity check
  cat("variable:", variable, "\n")
  cat("window:", toString(window), "\n")  # Convert list to a string
  cat("sample:", sample, "\n")
  cat("statistic:", statistic, "\n")
  
  cat("Processing Samples", "\n")
  # Extract the sample row from KYsamples for the given sample id
  sample_row <- meta.ky[meta.ky$sampleId_orig == sample, ]
  if (nrow(sample_row) == 0) {
    stop(paste("Sample", sample, "not found in KYsamples"))
  }
  cat("Processing Dates", "\n")
  # Convert the Collection_date to a Date object (assuming format "YEAR-MO-DY")
  collection_date <- as.Date(sample_row$Collection_date, format="%Y-%m-%d")
  
  # Define window boundaries:
  # Assuming window[1] is the offset for the end of the window (0 days before = collection day)
  # and window[2] is the offset for the start of the window (e.g., 7 days before)
  window_vec <- window[[1]] 
  print(window_vec[1])
  print(window_vec[2])
  window_end <- collection_date - window_vec[1]  
  window_start <- collection_date - window_vec[2]
  #cat("testing2", "\n")
  # Create a Date column in weather.data from YEAR, MO, DY (assuming they are numeric)
  
  # Filter weather.data to include only the dates in the desired window
  window_data <- weather.data[weather.data$date >= window_start & weather.data$date <= window_end, ]
  
  # Calculate the statistic's value of the specified variable in this window (ignore NA values)
  cat("Processing Statistics", "\n")
  if (statistic == "minimum"){
    calc_value <- min(window_data[[variable]], na.rm = TRUE)
  }
  if (statistic == "maximum"){
    calc_value <- max(window_data[[variable]], na.rm = TRUE)
  }
  if (statistic == "mean"){
    calc_value <- mean(window_data[[variable]], na.rm = TRUE)
  }
  if (statistic == "variance"){
    calc_value <- var(window_data[[variable]], na.rm = TRUE)
  }
  
  if (variable == "T2M"){
    if (statistic == "prop. min") {
      # Group the data by date and compute the daily minimum for the specified variable
      daily_min <- aggregate(window_data[[variable]], 
                             by = list(date = window_data$date), 
                             FUN = min, na.rm = TRUE)
      # Count the number of days where the minimum is below 5 degrees and divide by total days
      calc_value <- (sum(daily_min$x < 5)) / abs(window_vec[1] - window_vec[2])
    }
    
    if (statistic == "prop. max") {
      # Group the data by date and compute the daily maximum for the specified variable
      daily_max <- aggregate(window_data[[variable]], 
                             by = list(date = window_data$date), 
                             FUN = max, na.rm = TRUE)
      # Count the number of days where the maximum is above 32 degrees, the divide by total days
      calc_value <- (sum(daily_max$x > 32)) / abs(window_vec[1] - window_vec[2])
    }
  }
  else {
    if (statistic == "prop. max") {
      return()
    }
    if (statistic == "prop. min") {
      return()
    }}
  
  # Create a one-row data frame with the desired output structure
  cat("Processing Return Variable", "\n")
  result <- data.frame(
    Sample = sample,
    win = as.character(window),
    windowstart = as.character(window_start),
    windowend = as.character(window_end),
    variable = variable,
    stat = statistic,
    value = calc_value,
    stringsAsFactors = FALSE
  )
  cat("finishing...", "\n")
  print(calc_value)
  print(result)
  final_results <- rbind(final_results, result)  # Append results
  return(result)
}
final_results <- data.frame()

###
# Define window sizes
sets <- data.table(mod = 1:11,
                   start = c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end   = c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))


window_list <- lapply(1:nrow(sets), function(i) { 
  c(sets$start[i], sets$end[i])  # Store as numeric vector
})
start = c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0)
end   = c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90)

#head(window_list)

# Example variables and statistics
variables <- c("T2M", "PRECTOTCORR", "RH2M")
statistics <- c("mean", "maximum", "minimum", "prop. max", "prop. min", "variance")

combos <- CJ(win = window_list,
             var = variables, 
             stat = statistics, 
             samp = meta.ky$sampleId_orig, 
             sorted = FALSE)


