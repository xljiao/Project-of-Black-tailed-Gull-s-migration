# Load required packages
library(ggplot2)
library(dplyr)
library(stringr)   
library(geosphere) 
library(ade4)      
library(maps)      
library(ggpubr)    


# Load the dataset
df <- read.csv("mig.csv")

# Data cleaning and coordinate extraction pipeline
clean_data <- df %>%
  # Select relevant columns
  select(Group, Bird_year, Breeding_coordinates, Wintering_coordinates) %>%
  
  # Remove empty rows
  filter(!is.na(Breeding_coordinates), !is.na(Wintering_coordinates),
         Breeding_coordinates != "", Wintering_coordinates != "") %>%
  
  # Parse coordinate strings (Expected format: "39.528°N,123.042°E")
  # Use regular expressions to extract numeric values
  mutate(
    # Extract breeding coordinates
    Breed_Lat = as.numeric(str_extract(Breeding_coordinates, "^[0-9.]+(?=°N)")),
    Breed_Lon = as.numeric(str_extract(Breeding_coordinates, "(?<=,)[0-9.]+(?=°E)")),
    
    # Extract wintering coordinates
    Winter_Lat = as.numeric(str_extract(Wintering_coordinates, "^[0-9.]+(?=°N)")),
    Winter_Lon = as.numeric(str_extract(Wintering_coordinates, "(?<=,)[0-9.]+(?=°E)"))
  ) %>%
  
  # Filter out rows where parsing failed (NAs)
  filter(!is.na(Breed_Lat), !is.na(Breed_Lon), 
         !is.na(Winter_Lat), !is.na(Winter_Lon))

cat("Effective Sample Size (N):", nrow(clean_data), "\n")


# Overall Migratory Connectivity (Mantel Test)

# Calculate geographic distance matrices for breeding and wintering grounds
dist_breed <- distm(cbind(clean_data$Breed_Lon, clean_data$Breed_Lat), fun = distVincentyEllipsoid)
dist_winter <- distm(cbind(clean_data$Winter_Lon, clean_data$Winter_Lat), fun = distVincentyEllipsoid)

# Perform Mantel Test
mantel_res <- mantel.rtest(as.dist(dist_breed), as.dist(dist_winter), nrepet = 999)


# Population Comparison: Wintering Dispersion

# Calculate the wintering geometric centroid for each population
centroids <- clean_data %>%
  group_by(Group) %>%
  summarise(
    Center_Lon = mean(Winter_Lon),
    Center_Lat = mean(Winter_Lat)
  )

# Calculate distance from each individual to its population centroid
df_disp <- clean_data %>%
  left_join(centroids, by = "Group") %>%
  rowwise() %>%
  mutate(
    # Calculate point-to-point distance (converted to km)
    Dist_to_Center = distVincentyEllipsoid(c(Winter_Lon, Winter_Lat), 
                                           c(Center_Lon, Center_Lat)) / 1000 
  )

# Statistical Comparison
compare_res <- kruskal.test(Dist_to_Center ~ Group, data = df_disp)

# Visualization: Dispersion Boxplot
df_disp$Group <- factor(df_disp$Group, levels = c("Japan", "Russia", "Dalian", "Qingdao"))

my_comparisons <- list(c("Japan", "Russia"), c("Japan", "Dalian"), 
                       c("Japan", "Qingdao"), c("Russia", "Dalian"), 
                       c("Russia", "Qingdao"), c("Qingdao", "Dalian"))

p_disp <- ggplot(df_disp, aes(x = Group, y = Dist_to_Center, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  stat_compare_means(label.y.npc = "top") + # Add overall P-value
  labs(title = "Comparison of Wintering Dispersion (Weak Connectivity Proxy)",
       subtitle = "Higher distance = More scattered wintering = Weaker connectivity",
       y = "Distance to Wintering Centroid (km)",
       x = "Population") +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif",
                     method.args = list(exact = FALSE)) + 
  theme_classic() +
  scale_fill_brewer(palette = "Set2")

ggsave("Wintering_Dispersion_Comparison.pdf", p_disp, width = 6, height = 5)

# Visualization: Migratory Connectivity Map
world_map <- map_data("world")

p_map <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "lightgrey", color = "white") +
  coord_fixed(xlim = c(110, 145), ylim = c(15, 47), ratio = 1.3) +
  geom_segment(data = clean_data, aes(x = Breed_Lon, y = Breed_Lat, 
                                      xend = Winter_Lon, yend = Winter_Lat, 
                                      color = Group), 
               alpha = 0.6, linewidth = 0.5) +
  geom_point(data = clean_data, aes(x = Breed_Lon, y = Breed_Lat, color = Group), shape = 17, size = 2) + # Triangle = Breeding
  geom_point(data = clean_data, aes(x = Winter_Lon, y = Winter_Lat, color = Group), size = 2) + # Circle = Wintering
  
  scale_color_brewer(palette = "Set2") +
  labs(title = "Migratory Connectivity Map", x = "Longitude", y = "Latitude") +
  theme_bw()

ggsave("Mig_Connectivity_Map.pdf", p_map, width = 8, height = 6)