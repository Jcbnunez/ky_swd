# Miles Garvin
# Jan 29 2025
## SWD Weather Manipulatiuon Script

### libraries
library(tidyverse)
library(data.table)
library(gdata)
library(foreach)
library(nasapower)
library(sp)
library(lubridate)
library(doMC)
library(magrittr)
library(readr)




# get files
samples <- read_csv("~/Desktop/MASTER_SWD_METADATA.csv")
#View(dest_v2_samps_24Aug2024)
data <- samples
view(data)


## Rounds columns 3 and 4 to 2 decimal places

data <- data %>%
  mutate(
    lat = round(lat, 2),
    long = round(long, 2))


#####

# Find Each Unique City

# Get unique cities from column 5
unique_cities <- unique(data$city)

# Convert to a list
city_list <- as.list(unique_cities)

# Print the list
print(city_list)
head(city_list)



### Weather analysis ###


### set data definitions
query_parameters(community = "ag",
                 temporal_api = "hourly")

####################
## WEATHER DATA LOOP
# Fetches weather data for each city in city_list and adds to weather.data

weather.data <- foreach(city_name = city_list, .combine = rbind) %do% {   
  # Subset the data for the current city
  city_data <- subset(data, city == city_name) 
  
  if (nrow(city_data) == 0) {
    # Skip cities with no data
    return(NULL)
  }
  
  # Find the row with the earliest date
  earliest_row <- city_data[order(city_data$Collection_date), ][1, ]
  #head(earliest_row)
  
  # Format the date CORRECTLY
  year <- sub("^(\\d{4}).*", "\\1", earliest_row$Collection_date)
  
  # Extract remaining date part, keeping the dash
  month <- sub("^\\d{4}", "", earliest_row$Collection_date)
  month <- ifelse(nchar(month) == 3, paste0(month, "-01"), month)
  
print(month)
  
  # Print 4 debugging
  print(earliest_row)
  
  # Get weather data for each city
  temp.data <- tryCatch({
    get_power(
      community = "ag",
      lonlat = c(earliest_row$long, earliest_row$lat),
      pars = c("RH2M", "T2M", "PRECTOTCORR"),
      dates = c(paste(year, month, sep = ""), 
                paste(2024, "-08-31", sep = "")), # gets data from earliest occurence - summer 2024
      temporal_api = "hourly",
      time_standard = "UTC"
    ) %>% mutate(city = city_name)
  }, error = function(e) {
    message(sprintf("Failed to fetch data for city: %s. Error: %s", city_name, e))
    return(NULL)
  })
  
  # Return the data + add to weather.data
  temp.data
}

view(weather.data)
# Print the results
#view(weather.data)
#view(data)


##Write results  to a file!
write.csv(weather.data, "KYWeatherData.mgarvin.Jan25.csv", row.names = FALSE)
write.csv(data, "KYSampleData.mgarvin.Jan25.csv", row.names = FALSE)




############################
#KY Weather slicing!!
#Feb 28 2025

#Packages?
library(data.table)
library(foreach)
library(doParallel)


## Functions
calc_stat <- function(variable, window, sample, statistic) {
  
  #a lil sanity check
  cat("variable:", variable, "\n")
  cat("window:", toString(window), "\n")  # Convert list to a string
  cat("sample:", sample, "\n")
  cat("statistic:", statistic, "\n")
  
  cat("Processing Samples", "\n")
  # Extract the sample row from KYsamples for the given sample id
  sample_row <- KYsamples[KYsamples$sampleId_orig == sample, ]
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
  cat("testing2", "\n")
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


### Define Variables:

#main data sets
weather.data <- read_csv("KYWeatherData.mgarvin.Jan25.csv")
sample.data <- read_csv("SampleData.mgarvin.Jan25.csv")
#isolate KY samples
KYsamples <- sample.data[sample.data$province == "Kentucky", ]
# Remove rows with any NA values
KYsamples <- KYsamples[complete.cases(KYsamples), ]

# Define window sizes
sets <- data.table(mod = 1:11,
                   start = c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end   = c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))


window_list <- lapply(1:nrow(sets), function(i) { 
  c(sets$start[i], sets$end[i])  # Store as numeric vector
})
start = c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0)
end   = c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90)






head(window_list)

# Example variables and statistics
variables <- c("T2M", "PRECTOTCORR", "RH2M")
statistics <- c("mean", "maximum", "minimum", "prop. max", "prop. min", "variance")

combos <- CJ(win = window_list,
             var = variables, 
             stat = statistics, 
             samp = KYsamples$sampleId_orig, sorted = FALSE)


view(weather.data)

weather.data$date <- as.Date(with(weather.data, paste(YEAR, MO, DY, sep = "-")), format = "%Y-%m-%d")


#### Run the thing!!!!!
results <- foreach(i = 1:nrow(combos), .combine = rbind, .packages = "data.table") %do% { 
  combo <- combos[i, ]  # Extract the i-th row properly
  calc_stat(combo$var, combo$win, combo$samp, combo$stat)
}

view(results)

write.csv(results, "SlicedWeatherData.mgarvin.Feb25.csv", row.names = FALSE)
