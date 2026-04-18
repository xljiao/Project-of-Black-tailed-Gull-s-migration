
# Load required packages
library(geosphere)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggthemes)

# Read tracking data
df <- read.csv("tracking_data.csv")
colnames(df) <- c("pop", "bird", "id", "date", "stage", "lon", "lat")

# Unify date formats: some are formatted as "YYYY/MM/DD HH:MM", others as "YYYY/MM/DD"
# Append " 00:00" to dates lacking time information
df$datetime_clean <- ifelse(
  grepl("^\\d{4}/\\d{2}/\\d{2}$", df$date),
  paste0(df$date, " 00:00"),
  df$date
)

# Convert to POSIXct object
df$datetime_clean <- as.POSIXct(df$datetime_clean, format = "%Y/%m/%d %H:%M")

# Remove rows with invalid dates (NAs) generated during conversion
df <- df[!is.na(df$datetime_clean), ]
write.csv(df, file = "df.csv", row.names = FALSE)

# Function to calculate displacement from origin
calculate_nsd <- function(data) {
  # Ensure coordinates are numeric
  data$lon <- as.numeric(as.character(data$lon))
  data$lat <- as.numeric(as.character(data$lat))
  
  # Set the first point as the reference (origin)
  ref <- c(data$lon[1], data$lat[1])
  
  # Calculate the great-circle distance between each point and the origin

  distances_km <- distGeo(cbind(data$lon, data$lat), ref) / 1000^2
  data$nsd <- distances_km
  
  return(data)
}

df_nsd <- df %>%
  group_by(pop, bird) %>%
  group_modify(~ calculate_nsd(.x)) %>%
  ungroup()

# Save the calculated distance data
write.csv(df_nsd, file = "nsd_all_nosquare.csv", row.names = FALSE)

# Add a 'month_day' column for x-axis plotting.
# Convert dates to a "virtual year" (e.g., year 2000) using format() and as.Date() to align phenology
df_nsd$month_day <- format(df_nsd$datetime_clean, "%m-%d")
df_nsd$plot_date <- as.Date(paste0("2000-", df_nsd$month_day), format = "%Y-%m-%d")

# Create a "virtual date" column, retaining only month and day
df_nsd$month <- as.numeric(format(df_nsd$plot_date, "%m"))
df_nsd$day <- format(df_nsd$plot_date, "%d")

# Assign months May-Dec to year 2000, and Jan-Apr to year 2001 to ensure a continuous timeline
df_nsd$virtual_year <- ifelse(df_nsd$month >= 5, 2000, 2001)

# Construct the final virtual date column
df_nsd$virtual_date <- as.Date(paste0(df_nsd$virtual_year, "-", sprintf("%02d", df_nsd$month), "-", df_nsd$day))

# Generate plotting timelines for each population using a loop
populations <- unique(df_nsd$pop)

for (p in populations) {
  df_pop <- df_nsd %>% filter(pop == p)
  
  plot <- ggplot(df_pop, aes(x = virtual_date, y = nsd, color = bird, group = bird)) +
    geom_line(alpha = 0.8) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    labs(title = paste0("Annual displacement pattern - Population: ", p),
         x = "Month", y = "Distance from origin (km)") +
    theme_few() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  ggsave(filename = paste0("NSD_virtual_timeline_", p, ".pdf"), plot = plot, width = 8, height = 10)
}
