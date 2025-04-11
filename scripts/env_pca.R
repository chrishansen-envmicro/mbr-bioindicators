### Environmental PCA - MH vs. GNP ###

#Source of initial biplot script -> https://github.com/vqv/ggbiplot

#load libraries
library(tidyverse)
library(GGally)
library(ggbiplot)
library(wesanderson)

setwd("/Users/chrishansen/Desktop/hydro/hydro_git/data")
dat <- "Master_Data_Sheet_11112021_consolidated.csv"

# Import metadata
data <- read.csv(dat)

# Remove abbreviations and replace with site names
data$location <- recode(data$location, MH = "Mount Hood") %>%
  recode(GNP = "Glacier National Park")

# Filter to extract key variables and remove duplicate env values from biological replicates.
# Visibility into selected variables: https://docs.google.com/spreadsheets/d/1-Q---tBrG7Pugg2pc2jJyQNlYCLFRIXH/edit?gid=1346821870#gid=1346821870

data_filt <- select(data, location, sample_type, alt_sample_name, bio_rep, elevation_m, ysi_temp, ysi_do, ysi_spc,
                   ysi_ph, ysi_orp, alkalinity_caco3, bicarbonate_hco3,
                   calcium, chloride, magnesium, potassium, sio2, sodium,
                   strontium, sulfate, tds, fi_best) %>%
  
  # if bio_rep = 0; no biological replicates
  # if bio_rep = 1; 1st of 2+ biological replicates
  # if bio+rep = NA; No biological data
  
  filter(is.na(bio_rep) | bio_rep == 0 | bio_rep == 1) %>%
  filter(sample_type == "Spring" | sample_type == "Stream") %>%
  
  # Sulfate and chloride are the primary restrictions on n(total)
  drop_na(location, sample_type, elevation_m, ysi_temp, ysi_do, ysi_spc,
          ysi_ph, ysi_orp, alkalinity_caco3, bicarbonate_hco3,
          calcium, chloride, magnesium, potassium, sio2, sodium,
          strontium, sulfate, tds, fi_best)

# Separate site data
data_loc <- select(data_filt, location)
data_loc$location <- as.factor(data_loc$location)

# Separate environmental data
data_env <- select(data_filt, elevation_m, ysi_temp, ysi_do, ysi_spc,
                  ysi_ph, ysi_orp, alkalinity_caco3, 
                  bicarbonate_hco3, calcium, chloride,
                  magnesium, potassium, sio2, sodium,
                  strontium, sulfate, tds, fi_best) %>%
  
  # Rename variables for readability
  dplyr::rename(Elevation = elevation_m, Temp = ysi_temp, DO = ysi_do, SPC = ysi_spc,
                pH = ysi_ph, ORP = ysi_orp, CaCO3 = alkalinity_caco3, 
                Ca = calcium, HCO3 = bicarbonate_hco3, Cl = chloride,
                Mg = magnesium, K = potassium, SiO2 = sio2, 
                Na = sodium, Sr = strontium, SO4 = sulfate, TDS = tds, Fi = fi_best) %>%
  
  # If left-skew = x^2; If right-skew = log10(x). Informed by ggpairs(data_env)
  mutate(Temp = log10(Temp), DO = log10(DO), SPC = log10(SPC), CaCO3 = log10(CaCO3), 
         HCO3 = log10(HCO3), Ca = log10(Ca), Cl = log10(Cl), Mg = log10(Mg),
         K = log10(K), SiO2 = log10(SiO2), Na = log10(Na), Sr = log10(Sr), 
         SO4 = log10(SO4), TDS = log10(TDS))

# Assess normality of pairwise environmental variables
ggpairs(data_env)

# Remove multicollinearity across samples (R^2 >= 0.9; expections noted in Methods)
# Grp 1: CaCO3, HCO3, Mg, Ca, SPC
# Grp 2: SPC, TDS
# Grp 3: Sr, TDS
# Grp 4: Na, K, Cl, SiO2
vars_to_remove <- c("Cl", "HCO3", "CaCO3", "TDS", "Mg", "K", "SPC")
data_env_filt <- data_env %>% select(-all_of(vars_to_remove))

# Perform PCA
data_env_pca <- prcomp(data_env_filt, scale. = TRUE)

# Create biplot
p <- ggbiplot::ggbiplot(data_env_pca, obs.scale = 1, var.scale = 1,
                        group = data_loc$location, ellipse = TRUE, shape = data_loc$location,
                        ellipse.fill = TRUE, ellipse.prob = 0.95, ellipse.linewidth = 0,
                        point.size	= 3, varname.adjust = 1.5, ,
                        varname.size = 1, varname.color = "slategray", alpha = 0.7) +
  scale_color_manual(name = '', values=wes_palette("Moonrise2")) +
  scale_fill_manual(name = '', values=wes_palette("Moonrise2")) +
  scale_shape_manual(name = '', values=c(16,17)) +
  theme_light() +
  theme(axis.title = element_text(size = 14),
        legend.position="top",
        legend.text=element_text(size=12))

# Extract vector and text components for customization
seg <- which(sapply(p$layers, function(x) class(x$geom)[1] == 'GeomSegment'))
txt <- which(sapply(p$layers, function(x) class(x$geom)[1] == 'GeomText'))

# Customize vectors (thickness, color, arrows)
p$layers[[seg]] <- geom_segment(
  data = p$layers[[seg]]$data,
  aes(x = 0, y = 0, xend = xvar, yend = yvar),
  color = "slategray", # Vector color
  linewidth = 0.3, # Thickness (use 'size' in older ggplot2)
  arrow = arrow(length = unit(0.2, "cm") # Rounded line ends
  ))

# Modify text
p$layers[[txt]] <- geom_label(aes(x = xvar, y = yvar, label = varname,
                                  angle = angle, hjust = hjust), 
                              label.size = NA,
                              data = p$layers[[txt]]$data, 
                              fill = NA,
                              size=4)
p # print figure

### Run statistical analysis ###
library(vegan)
locs <- data_loc$location # Create locs variable
pca_scores <- as.data.frame(data_env_pca$x)  # Extract PC scores
pca_scores$locs <- data_loc$location  # Add group labels (e.g., MH vs. GNP)

# Euclidean distance on PCA scores (use first n PCs explaining most variance)
dist_matrix <- dist(pca_scores[,1:2])  # PC1 and PC2 only (or include more axes)

set.seed(503)
#Perform PERMANOVA to evaluate centroid differences
permanova_result <- adonis2(
  dist_matrix ~ locs,  # Group = MG vs. GNP
  data = pca_scores, 
  permutations = 999,  # Iterations
  method = "euclidean"  # Ideal distance metric for environmental variables; Bray-Curtis for 16S data
)

# View results
print(permanova_result) # Significant clustering differences

# Evalaute homogeneity of dispersion assumptions for MH and GNP around centroids
permdisp <- betadisper(dist_matrix, group = data_loc$location)
permutest(permdisp, pairwise = FALSE, permutations = 999)


# Separate environmental data, without correcting via log10 transformation
data_compare <- select(data_filt, elevation_m, ysi_temp, ysi_do, ysi_spc,
                  ysi_ph, ysi_orp, alkalinity_caco3, 
                  bicarbonate_hco3, calcium, chloride,
                  magnesium, potassium, sio2, sodium,
                  strontium, sulfate, tds, fi_best) %>%
  
  # Rename variables for readability
  dplyr::rename(Elevation = elevation_m, Temp = ysi_temp, DO = ysi_do, SPC = ysi_spc,
                pH = ysi_ph, ORP = ysi_orp, CaCO3 = alkalinity_caco3, 
                Ca = calcium, HCO3 = bicarbonate_hco3, Cl = chloride,
                Mg = magnesium, K = potassium, SiO2 = sio2, 
                Na = sodium, Sr = strontium, SO4 = sulfate, TDS = tds, Fi = fi_best)


data_compare$location <- as.factor(data_loc$location) # Add location column to data_env for summaries


#Summarize mean/median differences. Median differences applied for right-skewed data.
data_compare %>%
  group_by(location) %>%
  summarise(ele = mean(Elevation),
            temp = median(Temp),
            do = median(DO),
            spc = median(SPC),
            ph = mean(pH),
            orp = mean(ORP),
            ca = median(Ca),
            hco3 = median(HCO3),
            cl = median(Cl),
            mg = median(Mg),
            k = median(K),
            sio = median(SiO2),
            na = median(Na),
            sr = median(Sr),
            so4 = median(SO4),
            tds = median(TDS),
            caco3 = median(CaCO3))

### END ###