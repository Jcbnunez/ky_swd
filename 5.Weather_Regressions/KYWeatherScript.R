# Miles Garvin
# Feb 5 2025
## SWD Weather Manipulatiuon Script (KY Only)

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


##Unique City in kentucky!!!

# Get unique cities from column 5
unique_cities <- unique(data$city)

# Convert to a list
city_list <- as.list(unique_cities)

# Print the list
print(city_list)
head(city_list)
KY_city_list <- list("Lexington","Berea")


### Weather analysis ###


### set data definitions
query_parameters(community = "ag",
                 temporal_api = "hourly")

####################
## WEATHER DATA LOOP
# Fetches weather data for each city in city_list and adds to weather.data

weather.data <- foreach(city_name = KY_city_list, .combine = rbind) %do% {   
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
write.csv(weather.data, "KYWeatherData.mgarvin.csv", row.names = FALSE)


