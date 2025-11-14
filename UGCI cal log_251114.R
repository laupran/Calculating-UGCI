#install.packages("raster")
#install.packages("sf")
#install.packages("zoo")
#install.packages("dplyr")
#install.packages("terra")

library(raster)
library(sf)
library(zoo)
library(dplyr)
library(terra)
library(ggplot2)
library(tidyr)
library(broom)
library(ggpubr)

rm(list = ls())

VCS <- raster("C:/Users/user/Desktop/UGCI_raw/INPUT/vegetationCstorage_raw.tif") # tC/grid
VCU <- raster("C:/Users/user/Desktop/UGCI_raw/INPUT/netCuptake_raw.tif") # tC/ha/yr
SCS <- raster("C:/Users/user/Desktop/UGCI_raw/INPUT/soilCstorage_raw.tif") # tC/ha
ST <- raster("C:/Users/user/Desktop/UGCI_raw/INPUT/siltandclay_raw.tif") # %
LC <- raster("C:/Users/user/Desktop/UGCI_raw/INPUT/landcovermap_raw.tif") 

SCS_rp <- projectRaster(SCS, VCS, res=res(VCS), crs=crs(VCS), method="bilinear")
ST_rp <- projectRaster(ST, VCS, res=res(VCS), crs=crs(VCS), method="bilinear")
VCU_rp <- projectRaster(VCU, VCS, res=res(VCS), crs=crs(VCS), method="bilinear")
LC_rp <- projectRaster(LC, to = VCS, method = "ngb")

forest_codes <- c(1,2,3)
agri_codes <- c(4)
park_codes <- c(5)
roadside_codes <- c(6)
remove_codes <- c(7,8,9)

roadside_mask <- LC_rp %in% roadside_codes
agri_mask <- LC_rp %in% agri_codes
forest_mask <- LC_rp %in% forest_codes
park_mask <- LC_rp %in% park_codes
remove_mask <- LC_rp %in% remove_codes

#install.packages("igraph")
library(igraph)
forest_clumps <- clump(forest_mask, directions = 8)
freq_data <- freq(forest_clumps)
freq_data_df <- as.data.frame(freq_data)  

small_threshold <- 3
small_clumps <- freq_data_df[freq_data_df$count < small_threshold, "value"]

setwd("C:/Users/user/Desktop/UGCI_raw/Processed")
forest_filtered <- forest_mask
forest_filtered[forest_clumps %in% small_clumps] <- NA
plot(forest_filtered, main = "Filtered Forest Grid")
writeRaster(forest_filtered, "Forest_grid.tif", overwrite=TRUE)

forest_filtered[forest_filtered == 0] <- NA
agri_mask[agri_mask == 0] <- NA
park_mask[park_mask == 0] <- NA
roadside_mask[roadside_mask ==0] <- NA
remove_mask[remove_mask==0] <-NA

r_grid <- sum(roadside_mask, agri_mask, park_mask, forest_filtered, na.rm = TRUE)
r_grid[r_grid == 0] <- NA   
r_grid[r_grid > 0] <- 1     

plot(r_grid, main = "Combined Grid (Union of All Masks)")


writeRaster(r_grid, "Combined_mask_grid.tif", overwrite = TRUE)

#Preprocessing_masking
VCS_masked <- mask(VCS, r_grid)
VCU_masked <- mask(VCU_rp, r_grid)
SCS_masked <- mask(SCS_rp, r_grid)
ST_masked <- mask(ST_rp, r_grid)

VCS_masked_ha <- VCS_masked /900 *10000

writeRaster(VCS_masked_ha, filename = "VCS_masked_raw.tif", format = "GTiff", overwrite = TRUE)
writeRaster(VCU_masked, filename = "VCU_masked_raw.tif", format = "GTiff", overwrite = TRUE)
writeRaster(SCS_masked, filename = "SCS_masked_raw.tif", format = "GTiff", overwrite = TRUE)
writeRaster(ST_masked, filename = "ST_masked_raw.tif", format = "GTiff", overwrite = TRUE)



normal <- function(r) {
  name <- deparse(substitute(r))
  min_val <- min(values(r), na.rm = TRUE)
  max_val <- max(values(r), na.rm = TRUE)
  r_normalized <- (r - min_val) / (max_val - min_val)
  writeRaster(r_normalized, file=paste0("normalized_", name, ".tif"), format="GTiff", overwrite=TRUE)
}



VCS_masked_norm <- normal(VCS_masked)
VCU_masked_norm <- normal(VCU_masked)
SCS_masked_norm <- normal(SCS_masked)
ST_masked_norm <- normal(ST_masked)


data <- as.data.frame(stack(VCS_masked_norm, SCS_masked_norm, VCU_masked_norm, ST_masked_norm), xy = TRUE)
names(data) <- c("x", "y", "VCS", "SCS", "VCU", "ST")


data_clean <- data %>%
  filter(!is.na(VCS) & !is.na(SCS) & !is.na(VCU) & !is.na(ST)) %>%
  filter_all(all_vars(is.finite(.)))


decision_matrix <- data_clean %>% select(VCS, SCS, VCU, ST)
scale_minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
dm_scaled <- as.data.frame(lapply(decision_matrix, scale_minmax))

# p_{ij} calculation
col_sum <- colSums(dm_scaled)
p_ij <- sweep(dm_scaled, 2, col_sum, FUN = "/")  # 열 단위 나눗셈

# E_j, W_j calculation
m <- nrow(p_ij)
k <- 1 / log(m)

E_j <- sapply(p_ij, function(col_j) {
  col_j_non0 <- col_j[col_j > 0]
  -k * sum(col_j_non0 * log(col_j_non0))
})

d_j <- 1 - E_j
W_j <- d_j / sum(d_j)

W_j


UGCI_cal_weighted <- function(VCS_r, SCS_r, VCU_r, ST_r, W_j) {
  # W_j[1] = w_VCS, W_j[2] = w_SCS, W_j[3] = w_VCU, W_j[4] = w_ST
  w_VCS <- W_j[1]
  w_SCS <- W_j[2]
  w_vcu <- W_j[3]
  w_st  <- W_j[4]
  
  total_w <- sum(W_j)  
  
  weighted_sum <- ((w_VCS * VCS_r) + (w_SCS * SCS_r) + (w_vcu * VCU_r) + (w_st * ST_r))
  
  UGCI_w <- (weighted_sum / total_w)
  
  writeRaster(UGCI_w, filename = "C:/Users/user/Desktop/UGCI_raw/OUTPUT/UGCI_weighted.tif", format = "GTiff", overwrite = TRUE)
  
  plot(UGCI_w, main="UGCI (Weighted)")
  
  
  return(UGCI_w)
}

UGCI_masked_weighted <- UGCI_cal_weighted(
  VCS_r = VCS_masked_norm,
  SCS_r = SCS_masked_norm,
  VCU_r = VCU_masked_norm,
  ST_r  = ST_masked_norm,
  W_j   = W_j
)

vals <- values(UGCI_masked_weighted)

sum(!is.na(vals))


data_weighted <- as.data.frame(stack(VCS_masked_norm,
                                     SCS_masked_norm,
                                     VCU_masked_norm,
                                     ST_masked_norm,
                                     UGCI_masked_weighted),
                               xy = TRUE)
names(data_weighted) <- c("x", "y", "VCS", "SCS", "VCU", "ST", "UGCI")

data_processed <- data_weighted %>%
  filter_all(all_vars(is.finite(.))) %>%
  filter(!is.na(VCS) & !is.na(SCS) & !is.na(VCU) & !is.na(ST) & !is.na(UGCI))

data_processed <- data_processed %>%
  mutate(UGCI_level = cut(UGCI, 
                          breaks = quantile(UGCI, probs = c(0,0.25,0.5,0.75,1)), 
                          labels = c("Extremely low","Low","Moderate","High"),
                          include.lowest = TRUE))

category_ranges <- data_processed %>%
  group_by(UGCI_level) %>%
  summarize(Min_Value = min(UGCI), Max_Value = max(UGCI))

print(category_ranges)

# Visualization
library(ggplot2)
ggplot(data_processed, aes(x = x, y = y, fill = UGCI_level)) +
  geom_tile() +
  scale_fill_manual(values = c("blue", "green", "yellow", "red")) +
  labs(fill = "UGCI level") +
  theme_minimal() +
  coord_equal()


forest_fin <- mask(UGCI_masked_weighted, forest_filtered, maskvalues = NA)
agri_fin <- mask(UGCI_masked_weighted, agri_mask, maskvalues = NA)
park_fin <- mask(UGCI_masked_weighted, park_mask, maskvalues = NA)
roadside_fin <- mask(UGCI_masked_weighted, roadside_mask, maskvalues = NA)

plot(roadside_fin)


writeRaster(forest_fin, filename = "C:/Users/user/Desktop/UGCI_raw/OUTPUT/forest_fin_UGCI.tif", format = "GTiff", overwrite = TRUE)
writeRaster(agri_fin, filename = "C:/Users/user/Desktop/UGCI_raw/OUTPUT/agri_fin_UGCI.tif", format = "GTiff", overwrite = TRUE)
writeRaster(park_fin, filename = "C:/Users/user/Desktop/UGCI_raw/OUTPUT/park_fin_UGCI.tif", format = "GTiff", overwrite = TRUE)
writeRaster(roadside_fin, filename = "C:/Users/user/Desktop/UGCI_raw/OUTPUT/roadside_fin_UGCI.tif", format = "GTiff", overwrite = TRUE)


forest_data <- as.data.frame(forest_fin, xy = TRUE) %>%
  rename(UGCI = "layer") %>%
  mutate(LandCover = "Forest")

agri_data <- as.data.frame(agri_fin, xy = TRUE) %>%
  rename(UGCI = "layer") %>%
  mutate(LandCover = "Agricultural area")

park_data <- as.data.frame(park_fin, xy = TRUE) %>%
  rename(UGCI = "layer") %>%
  mutate(LandCover = "Urban park")

roadside_data <- as.data.frame(roadside_fin, xy = TRUE) %>%
  rename(UGCI = "layer") %>%
  mutate(LandCover = "Roadside vegetation")


combined_data <- bind_rows(forest_data, agri_data, park_data, roadside_data)


clean_data <- combined_data %>%
  filter(!is.na(UGCI) & is.finite(UGCI))


processed_data <- clean_data %>%
  mutate(UGCI_level = cut(UGCI, breaks = quantile(UGCI, probs = c(0, 0.25, 0.5, 0.75, 1)),
                          labels = c("Extremely Low", "Low", "Moderate", "High"), include.lowest = TRUE))


land_cover_summary <- processed_data %>%
  group_by(LandCover, UGCI_level) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count)) #%>%
 # mutate(Proportion = format(Proportion, digits = 8, nsmall = 6))


print(land_cover_summary)

library(scales)

land_cover_summary_plot <- land_cover_summary %>%
  mutate(
    LandCover = factor(LandCover, levels = c("Forest", "Agricultural area", "Urban park", "Roadside vegetation")),
    UGCI_level = factor(UGCI_level, levels = c("Extremely Low", "Low", "Moderate","High"))
  )

# Stacked bar chart (proportion graph)
ggplot(land_cover_summary_plot, aes(x = LandCover, y = Proportion, fill = UGCI_level)) +
  geom_bar(stat = "identity") +
  labs(x = "Land cover type", y = "Proportion", fill = "UGCI level") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Extremely Low" = "#CB181D", "Low" = "#FB6A4A", "Moderate" = "#FCAE91", "High" = "#FEE5D9")) +
  theme_minimal()





land_cover_mean <- processed_data %>%
  group_by(LandCover) %>%
  summarise(Average_UGCI = mean(UGCI, na.rm = TRUE))

print(land_cover_mean)


palette_list <- list(
  "Forest" = c("#e6ffd6", "#caffa8", "#7bdb3d", "#2a6b00"),
  "Agricultural area" = c("#ffefd6", "#ffaa73", "#fc7e2b", "#d45400"),
  "Urban park" = c("#d6ecff", "#8cb3ff", "#2e5cb8", "#0d1494"),
  "Roadside vegetation" = c("#e2ccff","#c294ff","#8334eb","#420a8c")
)

level_index <- c("Extremely low" = 1, "Low" = 2, "Moderate" = 3, "High" = 4)


processed_data_color <- processed_data %>%
  mutate(Color = mapply(function(landcover, level) {
    palette_list[[landcover]][level_index[[level]]]
  }, LandCover, UGCI_level))


unique_data <- processed_data_color %>%
  distinct(Color, LandCover, UGCI_level)

color_breaks <- unique_data$Color
color_labels <- paste(unique_data$LandCover, unique_data$UGCI_level, sep = " - ")

ggplot(processed_data_color, aes(x = x, y = y, fill = Color)) +
  geom_tile() +
  scale_fill_identity(
    name = "UGCI Level",
    breaks = color_breaks,
    labels = color_labels,
    guide = "legend"
  ) +
  labs(title = "Integrated UGCI Levels by Land Cover Type") +
  theme_bw() +
  coord_equal()


values <- values(forest_fin)
values <- values[!is.na(values)]
hist(values, 
     main = "Distribution of CMPI Values:roadside", 
     xlab = "CMPI", 
     ylab = "Frequency", 
     col = "lightblue", 
     border = "white", 
     breaks = 50,
     xlim = c(0,1))  


##### Urban park processing #####################
# Calculate patch sizes for Urban Parks
park_clumps <- clump(park_fin, directions = 8)  # Group pixels into patches
park_freq <- freq(park_clumps)  # Frequency of clump IDs
park_patch_sizes <- as.data.frame(park_freq) %>%
  filter(!is.na(value)) %>%
  mutate(Area_ha = count * (res(park_fin)[1] * res(park_fin)[2] / 10000))  # Area in hectares

# Average UGCI per patch
park_ugci_zonal <- zonal(park_fin, park_clumps, fun = 'mean', na.rm = TRUE)
park_ugci_data <- as.data.frame(park_ugci_zonal) %>%
  filter(!is.na(zone)) %>%
  rename(Patch_ID = zone, Avg_UGCI = mean)

# Merge CMPI data with patch sizes
park_analysis <- merge(park_patch_sizes, park_ugci_data, by.x = "value", by.y = "Patch_ID")

# Correlation
correlation <- cor(park_analysis$Area_ha, park_analysis$Avg_UGCI, use = "complete.obs")
print(correlation)

ggplot(park_analysis, aes(x = log(Area_ha), y = Avg_UGCI)) +
  geom_point(color = "gray", size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(x = "Log(Area, ha)", y = "UGCI") +
  theme_minimal()



# TERNARY PLOT
new_for <- as.data.frame(stack(forest_fin, VCS_masked_norm, VCU_masked_norm, SCS_masked_norm), xy = TRUE)
names(new_for) <- c("x", "y", "UGCI", "VCS","CUP", "SCS")

new_for_cleaned <- new_for %>%
  filter(!is.na(UGCI) & !is.na(VCS) & !is.na(CUP) & !is.na(SCS)) %>%
  filter_all(all_vars(is.finite(.))) %>%
  mutate(
    UGCI_level = cut(
      UGCI, 
      breaks = c(0, 0.14, 0.26, 0.49, 1),
      labels = c("Extremely Low", "Low", "Moderate", "High"),
      include.lowest = TRUE
    )
  )

shifted_data <- new_for_cleaned %>%
  filter(UGCI_level == c("Extremely Low", "Low"))

library(ggtern)
ggtern(shifted_data, aes(x = VCS, y = CUP, z = SCS, color = UGCI_level)) +
  theme_rgbw() +
  geom_point(alpha = 1) +
  labs(title = "Forest", x = "Vegetation C storage", y = "Net C uptake", z = "Soil C storage") +
  theme(
    axis.title = element_text(size = 0),  
    axis.text = element_text(size = 18),   
    #legend.position = "none"               
  )



new_agr <- as.data.frame(stack(agri_fin, VCS_masked_norm, VCU_masked_norm, SCS_masked_norm), xy = TRUE)
names(new_agr) <- c("x", "y", "UGCI", "VCS","CUP", "SCS")

new_agr_cleaned <- new_agr %>%
  filter(!is.na(UGCI) & !is.na(VCS) & !is.na(CUP) & !is.na(SCS)) %>%
  filter_all(all_vars(is.finite(.))) %>%
  mutate(
    UGCI_level = cut(
      UGCI, 
      breaks = c(0, 0.14, 0.26, 0.49, 1),
      labels = c("Extremely Low", "Low", "Moderate", "High"),
      include.lowest = TRUE
    )
  )

shifted_data <- new_agr_cleaned %>%
  filter(UGCI_level == c("Extremely Low", "Low"))

library(ggtern)
ggtern(shifted_data, aes(x = VCS, y = CUP, z = SCS, color = UGCI_level)) +
  theme_rgbw() +
  geom_point(alpha = 1) +
  labs(title = "Agri", x = "Vegetation C storage", y = "Net C uptake", z = "Soil C storage") +
  scale_color_manual(
    values = c(
      "Extremely Low" = "#00BFC4",
      "Low" = "#F8766D"
    )
  ) +
  theme(
    axis.title = element_text(size = 0), 
    axis.text = element_text(size = 18),  
    #legend.position = "none"             
  )



new_prk <- as.data.frame(stack(park_fin, VCS_masked_norm, VCU_masked_norm, SCS_masked_norm), xy = TRUE)
names(new_prk) <- c("x", "y", "UGCI", "VCS","CUP", "SCS")

new_prk_cleaned <- new_prk %>%
  filter(!is.na(UGCI) & !is.na(VCS) & !is.na(CUP) & !is.na(SCS)) %>%
  filter_all(all_vars(is.finite(.))) %>%
  mutate(
    UGCI_level = cut(
      UGCI, 
      breaks = c(0, 0.14, 0.26, 0.49, 1),
      labels = c("Extremely Low", "Low", "Moderate", "High"),
      include.lowest = TRUE
    )
  )

shifted_data <- new_prk_cleaned %>%
   filter(UGCI_level == c("Extremely Low", "Low"))

library(ggtern)
ggtern(shifted_data, aes(x = VCS, y = CUP, z = SCS, color = UGCI_level)) +
  theme_rgbw() +
  geom_point(alpha = 1) +
  labs(title = "Park", x = "Vegetation C storage", y = "Net C uptake", z = "Soil C storage") +
  scale_color_manual(
    values = c(
      "Extremely Low" = "#00BFC4",
      "Low" = "#F8766D"
    )
  ) +
  theme(
  axis.title = element_text(size = 0),  
  axis.text = element_text(size = 18),   
  #legend.position = "none"              
)






new_road <- as.data.frame(stack(roadside_fin, VCS_masked_norm, VCU_masked_norm, SCS_masked_norm), xy = TRUE)
names(new_road) <- c("x", "y", "UGCI", "VCS","CUP", "SCS")

new_road_cleaned <- new_road %>%
  filter(!is.na(UGCI) & !is.na(VCS) & !is.na(CUP) & !is.na(SCS)) %>%
  filter_all(all_vars(is.finite(.))) %>%
  mutate(
    UGCI_level = cut(
      UGCI, 
      breaks = c(0, 0.14, 0.26, 0.49, 1),
      labels = c("Extremely Low", "Low", "Moderate", "High"),
      include.lowest = TRUE
    )
  )

shifted_data <- new_road_cleaned %>%
  filter(UGCI_level == c("Extremely Low", "Low"))

library(ggtern)
ggtern(shifted_data, aes(x = VCS, y = CUP, z = SCS, color = UGCI_level)) +
  theme_rgbw() +
  geom_point(alpha = 1) +
  labs(title = "Road", x = "Vegetation C storage", y = "Net C uptake", z = "Soil C storage") +
  scale_color_manual(
    values = c(
      "Extremely Low" = "#00BFC4",
      "Low" = "#F8766D"
    )
  ) +
  theme(
    axis.title = element_text(size = 0), 
    axis.text = element_text(size = 18),
  )



# --- 1. Identify Forest Patches and Their Boundaries ---
forest_clumps <- clump(forest_filtered, directions = 8)

# Extract the boundaries (edge cells) of these forest patches.
forest_boundary <- boundaries(forest_clumps, type = "inner", directions = 8)

# Create a mask for forest edge cells: retain only cells with a boundary value of 1.
forest_edge_mask <- forest_boundary
forest_edge_mask[forest_edge_mask != 1] <- NA
plot(forest_edge_mask)

# --- 2. Define Forest Edge and Interior Rasters ---
forest_edge <- mask(forest_filtered, forest_edge_mask)

forest_interior <- forest_filtered
forest_interior[!is.na(forest_edge_mask)] <- NA
plot(forest_interior)

forest_edge_ugci <- mask(UGCI_masked_weighted, forest_edge_mask)

# Extract CMPI values corresponding to forest interior cells.
forest_interior_ugci <- mask(UGCI_masked_weighted, forest_interior)

# Combine forest CMPI data (edge and interior)
forest_edge_df <- as.data.frame(forest_edge_ugci, xy = TRUE) %>%
  rename(UGCI = layer) %>%
  mutate(LandCover = "Forest", Location = "Edge")

forest_int_df <- as.data.frame(forest_interior_ugci, xy = TRUE) %>%
  rename(UGCI = layer) %>%
  mutate(LandCover = "Forest", Location = "Non-edge")


# Urban park edge vs. interior extraction using park_fin and park_clumps
park_boundary <- boundaries(park_clumps, type = "inner", classes = TRUE, directions = 8)
park_edge_raster <- park_boundary
park_edge_raster[park_edge_raster != 1] <- NA  # keep only boundary cells

park_interior_raster <- park_boundary
park_interior_raster[park_interior_raster != 0] <- NA  # keep only interior cells

park_edge <- mask(park_fin, park_edge_raster)
park_int  <- mask(park_fin, park_interior_raster)

park_edge_df <- as.data.frame(park_edge, xy = TRUE) %>%
  rename(UGCI = layer) %>%
  mutate(LandCover = "Urban park", Location = "Edge")

park_int_df <- as.data.frame(park_int, xy = TRUE) %>%
  rename(UGCI = layer) %>%
  mutate(LandCover = "Urban park", Location = "Non-edge")


# Combine data frames from Forest and Urban park
combined_box <- bind_rows(
  forest_edge_df, forest_int_df,
  park_edge_df, park_int_df
) %>%
  filter(!is.na(UGCI) & is.finite(UGCI))

# Ensure the factors are ordered correctly for plotting
combined_box <- combined_box %>%
  mutate(
    LandCover = factor(LandCover, levels = c("Forest", "Urban park")),
    Location  = factor(Location, levels = c("Edge", "Non-edge"))
  )

# Create the box plot with ggpubr
ggboxplot(combined_box, 
          x = "Location", y = "UGCI", 
          color = "Location", palette = "jco",
          add = "jitter", 
          add.params = list(alpha = 0.05, size = 0.5),  # Set jitter transparency and size
          facet.by = "LandCover", short.panel.labs = TRUE) +
  # Add mean value as a red diamond
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "white", color = "red", show.legend = FALSE) +
  # Add mean value as text labels (rounded to 2 decimals) in bold
  stat_summary(fun = mean, geom = "text", 
               aes(label = round(..y.., 2)), 
               vjust = -1, color = "black", fontface = "bold", show.legend = FALSE) +
  # Add t-test significance (stars) between Edge and Interior for each land cover type
  stat_compare_means(method = "t.test", label = "p.signif", hide.ns = TRUE) +
  labs(x="",y = "UGCI Value")


impervious_grid <- mask(LC_rp, remove_mask)

forest_edge_dist <- distance(forest_edge_mask)
park_edge_dist <- distance(park_edge)
impervious_dist <- distance(impervious_grid)

threshold <- 30


forest_impervious_edge <- forest_edge_mask
forest_impervious_edge[impervious_dist > threshold] <- NA


forest_park_edge <- forest_edge_mask
forest_park_edge[forest_impervious_edge] <- NA


park_impervious_edge <- park_edge
park_impervious_edge[impervious_dist > threshold] <- NA


park_forest_edge <- park_edge
park_forest_edge[park_impervious_edge] <- NA


forest_park_edge_ugci <- mask(UGCI_masked_weighted, forest_park_edge)
forest_impervious_edge_ugci <- mask(UGCI_masked_weighted, forest_impervious_edge)
park_forest_edge_ugci <- mask(UGCI_masked_weighted, park_forest_edge)
park_impervious_edge_ugci <- mask(UGCI_masked_weighted, park_impervious_edge)


df_forest_park <- as.data.frame(forest_park_edge_ugci, xy = TRUE) %>%
  rename(UGCI = layer) %>%
  mutate(AdjacentLC = "Park", Category = "Forest edge") %>%
  filter(!is.na(UGCI) & is.finite(UGCI))

df_forest_impervious <- as.data.frame(forest_impervious_edge_ugci, xy = TRUE) %>%
  rename(UGCI = layer) %>%
  mutate(AdjacentLC = "Impervious", Category = "Forest edge") %>%
  filter(!is.na(UGCI) & is.finite(UGCI))


df_park_forest <- as.data.frame(park_forest_edge_ugci, xy = TRUE) %>%
  rename(UGCI = layer) %>%
  mutate(AdjacentLC = "Forest", Category = "Urban park edge") %>%
  filter(!is.na(UGCI) & is.finite(UGCI))

df_park_impervious <- as.data.frame(park_impervious_edge_ugci, xy = TRUE) %>%
  rename(UGCI = layer) %>%
  mutate(AdjacentLC = "Impervious", Category = "Urban park edge") %>%
  filter(!is.na(UGCI) & is.finite(UGCI))


combined_box <- bind_rows(df_forest_park, df_forest_impervious,
                          df_park_forest, df_park_impervious) %>%
  mutate(
    EdgeType = case_when(
      Category == "Forest edge" & AdjacentLC == "Park" ~ "Forest-Park",
      Category == "Forest edge" & AdjacentLC == "Impervious" ~ "Forest-Impervious",
      Category == "Urban park edge" & AdjacentLC == "Forest" ~ "Park-Forest",
      Category == "Urban park edge" & AdjacentLC == "Impervious" ~ "Park-Impervious"
    ),
    Category = factor(Category, levels = c("Forest edge", "Urban park edge")),
    EdgeType = factor(EdgeType, 
                      levels = c("Forest-Park", "Forest-Impervious",
                                 "Park-Forest", "Park-Impervious"))
  )


ggboxplot(combined_box, 
          x = "EdgeType", y = "UGCI", 
          color = "AdjacentLC", palette = c("green4", "orange", "green4", "orange"),
          add = "jitter", 
          add.params = list(alpha = 0.1, size = 1)) +
  facet_wrap(~Category, scales = "free_x") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "white", color = "red", show.legend = FALSE) +
  stat_summary(fun = mean, geom = "text", 
               aes(label = round(..y.., 2)), 
               vjust = -1, color = "black", fontface = "bold", show.legend = FALSE) +
 
  stat_compare_means(aes(x = EdgeType, y = UGCI), 
                     method = "t.test", label = "p.signif", hide.ns = TRUE) +
  labs(x = "", y = "UGCI Value")

park_edge_data <- combined_box %>% filter(Category == "Urban park edge")

t_test_park <- t.test(UGCI ~ EdgeType, data = park_edge_data)
print(t_test_park)

