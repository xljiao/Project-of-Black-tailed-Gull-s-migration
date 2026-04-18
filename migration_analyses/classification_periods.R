# This script was modified from the reference: 
# Xia, H., Nilsson, C., Thorup, K., Jia, C. & Lei, F. Non-breeding movements of the Black-tailed Gull (Larus crassirostris). Avian Research 14, 100103, doi:10.1016/j.avrs.2023.100103 (2023)

# Load required packages

library(lubridate)
library(tidyverse)
library(ggpubr)

#Function to calculate haversine distances between points (taken from stack overflow...)
dt.haversine <- function(lat_from, lon_from, lat_to, lon_to, r = 6378137) { radians <- pi / 180
  lat_to <- lat_to * radians
  lat_from <- lat_from * radians
  lon_to <- lon_to * radians
  lon_from <- lon_from * radians
  dLat <- (lat_to - lat_from)
  dLon <- (lon_to - lon_from)
  a <- (sin(dLat / 2)^2) + (cos(lat_from) * cos(lat_to)) * (sin(dLon / 2)^2)
  return(2 * atan2(sqrt(a), sqrt(1 - a)) * r) }
# Load gull data
gulls <- read.csv(input_file, stringsAsFactors = FALSE, header = T)
gulls[is.na(gulls)] <- 0
gulls %>%
  mutate(time = as.POSIXct(Timestamp, format = "%Y/%m/%d %H:%M")) %>%
  mutate(time = format(as.POSIXct(time), format = "%M")) -> gulls
gulls[with(gulls, time < 30),] -> gulls
gulls %>%
  mutate(Date = as.POSIXct(Timestamp, format = "%Y/%m/%d %H:%M")) %>%
  mutate(Date = round_date(Date, "hour")) -> gulls1
gulls1[!duplicated(gulls1$Date),] -> gulls1
gulls1 %>%
  mutate(log_prev = NA) %>%
  mutate(lat_prev = NA) -> gulls1
for (i in c(2:nrow(gulls1))) {
  gulls1$log_prev[i] = gulls1$Longitude[i - 1]
  gulls1$lat_prev[i] = gulls1$Latitude[i - 1]
}

# Create a date vector with all days, to NA fill data
date_vector <- as.data.frame(seq(from = min(gulls1$Date, na.rm = TRUE), to = max(gulls1$Date, na.rm = TRUE), by = "hour"))
colnames(date_vector) <- "Date"

# Make sure that "gulls" are in chronological order
gulls1 %>%
  full_join(date_vector, by = "Date") %>%
  arrange(Date) -> gulls2

### SET NR OF TIMESTAMPS(NR OF TIMESTAMPS CORRESPONDING TO HOURS) AND BUFFER DISTANCE ###
t = 192
km = 15
dist_buffer = km * 1000

#Calculate New segment if moved #Correct segment for returns
gulls2 %>%
  rownames_to_column() %>%
  mutate(dist = dt.haversine(lat_prev, log_prev, Latitude, Longitude)) %>%
  mutate(Newsegment = ifelse(dist < dist_buffer, 0, rowname)) %>%
  mutate(Newsegment = ifelse(is.na(Newsegment), 0, Newsegment)) %>%
  mutate(Newsegment = cumsum(Newsegment)) %>%
  mutate(Rsegment = NA) -> gulls3
gulls3$Rsegment[c(1:t)] = 0
for (n in c(c(t + 1):nrow(gulls3))) {
  #
  if (is.na(gulls3$Latitude[n])) {
    gulls3$Rsegment[n] = gulls3$Rsegment[n - 1]
    next
  }
  gull <- gulls3[c(c(n - t):c(n - 1)),]
  gull <- na.omit(gull)
  #
  if (nrow(gull) == 0) {
    gulls3$Rsegment[n] = ifelse(gulls3$Newsegment[n] == gulls3$Newsegment[n - 1],
                                gulls3$Rsegment[n - 1], gulls3$Newsegment[n])
    next
  }
  for (j in c(1:nrow(gull))) {
    dis = dt.haversine(gulls3$Latitude[n], gulls3$Longitude[n], gull$Latitude[j], gull$Longitude[j])
    #
    if (dis < dist_buffer) {
      gulls3$Rsegment[n] = gull$Rsegment[j]
      break
    }
    #
    if (j == nrow(gull)) {
      gulls3$Rsegment[n] = gulls3$Newsegment[n]
    }
  }
}

#Calculate segment
gulls3 %>% mutate(segment = NA) -> gulls4
gulls4$segment[1] = 0

for (k in c(2:t)) {
  if (min(gulls4$Rsegment[c(c(k + 1):c(k + t))]) < gulls4$Rsegment[k]) {
    gulls4$segment[k] = min(gulls4$Rsegment[c(c(k + 1):c(k + t))])
    next
  }
  
  gulls4$segment[k] = gulls4$Rsegment[k]
}

for (m in c(c(t + 1):c(nrow(gulls4) - t))) {
  if (min(gulls4$Rsegment[c(c(m + 1):c(m + t))]) < gulls4$Rsegment[m]) {
    gulls4$segment[m] = gulls4$segment[c(m - 1)] 
    next
  }

  if (length(which(gulls4$Rsegment[c(c(m - t):c(m - 1))] == gulls4$Rsegment[m])) > 0) {
    gulls4$segment[m] = gulls4$segment[c(m - 1)]
    next
  }
  
  gulls4$segment[m] = gulls4$Rsegment[m]
}

for (b in c(c(nrow(gulls4) - t + 1):c(nrow(gulls4) - 1))) {
  if (min(gulls4$Rsegment[c(c(b + 1):nrow(gulls4))]) < gulls4$Rsegment[b]) {
    gulls4$segment[b] = gulls4$segment[c(b - 1)]
    next
  }

  if (length(which(gulls4$Rsegment[c(c(b - t):c(b - 1))] == gulls4$Rsegment[b])) > 0) {
    gulls4$segment[b] = gulls4$segment[c(b - 1)]
    next
  }
  
  gulls4$segment[b] = gulls4$Rsegment[b]
}

if (length(which(gulls4$Rsegment[c(c(nrow(gulls4) - t):c(nrow(gulls4) - 1))] == gulls4$Rsegment[nrow(gulls4)])) > 0) {
  gulls4$segment[nrow(gulls4)] = gulls4$segment[c(nrow(gulls4) - 1)] 
} else {
  gulls4$segment[nrow(gulls4)] = gulls4$Rsegment[nrow(gulls4)] 
}
# Identify migration periods(Borrow the idea of density clustering)
# Calculate Nr(density) of each cluster and plot line chart to find out kickpoint(threshold)

seg <- unique(gulls4$segment)
segdensity <- data.frame(segment = seg, density = NA)
for (h in c(1:length(seg))) { segdensity$density[h] = length(which(gulls4$segment == seg[h])) }
newseg <- segdensity[order(segdensity$density, decreasing = T),]
nsd <- data.frame(newseg$density)
colnames(nsd) = "density"
nsd %>% rownames_to_column() -> nsd
nsd$rowname <- as.integer(rownames(nsd))

###SET THRESHOLD DENSITY###
den_buffer = 65

# Join the clusters those densities are below the threshold(equal to noise point in density clustering) together
migration = newseg[newseg$density <= den_buffer,]$segment
gulls4 %>% mutate(segment_sp = segment) -> gull_segment_sp
for (f in c(1:length(migration))) { mig = migration[f]
  gull_segment_sp$segment_sp[gull_segment_sp$segment_sp == mig] = -1 }
gull_segment_sp <- gull_segment_sp %>% drop_na()
gull_segment_sp$stage <- paste0("stage", cumsum(c(TRUE, diff(gull_segment_sp$segment_sp) != 0)))

write.csv(gull_segment_sp, output_file, row.names = FALSE)

