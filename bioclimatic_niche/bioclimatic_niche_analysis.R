# Hillieae bioclimatic niche description using niche package and worldclim bioclimatic variables + elevation (2.5 min resolution)
# niche package: https://rdrr.io/rforge/ade4/man/niche.html

setwd("/Users/laymonball/Desktop/niche_analysis")

library(dismo)
library(raster)
library(ade4)
library(ggplot2)
library(factoextra)
library(adegraphics)
library(ggrepel)
library(ENMTools)
library(virtualspecies)
library(maps)
library(sp)
library(devtools)
library(ggforce)
library(dplyr)
library(sf)
library(mixtools)
library(gplots)
library(scales)
library(ggtreeExtra)


########################################################################
## First, need to test for correlations between 20 Worldclim variables ##
########################################################################
setwd("~/Desktop")

# Load the environmental variables (resolution: 2.5 arc-minute)
files <- list.files(path = paste("niche_analysis/wc2.1_2.5m_bio/", sep = ""), pattern = "tif$", full.names = TRUE)
files

predictors <- stack(files)
# predictors
names(predictors) <- c("bio01","bio10","bio11","bio12","bio13","bio14","bio15","bio16",
                       "bio17","bio18","bio19","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09",
                       "elev")

# Assess the rasters
# Calculate correlation coefficient between rasters, and then plot; can decide whether to remove highly correlated variables
pears_cor_plots <- raster.cor.plot(predictors, method = "pearson")
# pears_cor_plots
# pears_cor_plots$cor.heatmap$data
# write.table(pears_cor_plots$cor.heatmap$data, "niche_analysis/wc2.1_2.5m_bio/clipped_asc/pears_cor_plot_heatmapVals.txt", quote = FALSE, sep = "\t")

vars_collin_tree <- removeCollinearity(
  predictors,
  multicollinearity.cutoff = 0.6,
  select.variables = FALSE,
  sample.points = FALSE,
  plot = TRUE,
  method = "pearson"
)

# Based on above outputs, I'm going to keep bio1, bio4, bio12, bio15, and elevation.
vars_filt <- subset(predictors, c(1, 4, 7, 14, 20))
# vars_filt

vars_filt_dframe <- as.data.frame(subset(predictors, c(1, 4, 7, 14, 20)))



#########################################################
## Extract bioclim values for species occurrence points ##
#########################################################
# Load rasters
ras_files <- list.files(path = paste("~/Desktop/niche_analysis/wc2.1_2.5m_bio/filtered_clipped/", sep = ""), pattern = "sdat$", full.names = TRUE)
predictors <- stack(ras_files)
names(predictors) <- c("bio01","bio12","bio15", "bio04","elev")

# Load species coordinates and extract bioclim values for each coord
points <- unique(read.csv("~/Desktop/niche_analysis/all_Hillieae_coords.csv", h = T))[1:3]
values <- raster::extract(predictors, SpatialPoints(points[c(2:3)]), df = T)
niche_dat <- na.omit(cbind(values[2:6], points[1]))

# Take log10 of bioclim values, except for elevation
niche_dat[1:4] <- log10(niche_dat[1:4])

# Calculate median value for climatic variables for each species
median_niche_values <- aggregate(niche_dat[1:5], niche_dat[6], median)
# write.table(median_niche_values, "niche_analysis/median_climate_data.csv", quote = FALSE, sep = "\t")

# mean_niche <- aggregate(niche_dat[1:5], niche_dat[6], mean)
# SD <- aggregate(niche_dat[1:5], niche_dat[6], sd)



#########################################################
## Describe niche following methods from Alexandre 2017 ##
#########################################################
# Script adapted from Alexandre_2017_supplementary_data/bioclimatic_niches_script.html

# Load the environmental variables (resolution: 2.5 arc-minute)
ras_files <- list.files(path = paste("~/Desktop/niche_analysis/wc2.1_2.5m_bio/filtered_clipped/", sep = ""), pattern = "sdat$", full.names = TRUE)
predictors <- stack(ras_files)
# predictors
names(predictors) <- c("bio01", "bio12", "bio15", "bio04", "elev")
pred <- as.data.frame(predictors)
PRED <- cbind(pred, rownames(pred))

# Load region polygon shapefiles
Amazon <- st_read("~/Desktop/niche_analysis/region_polygon_shapefiles/Amazon_dissolved.shp")
# plot(Amazon)
Andes <- st_read("~/Desktop/niche_analysis/region_polygon_shapefiles/Andes_dissolved.shp")
NCA <- st_read("~/Desktop/niche_analysis/region_polygon_shapefiles/NCA_dissolved.shp")
SCA <- st_read("~/Desktop/niche_analysis/region_polygon_shapefiles/SCA_dissolved.shp")

# Get lat/long coordinates from raster cells
PRED$coords <- xyFromCell(predictors, as.numeric(PRED$`rownames(pred)`))
PRED$region <- NA
# Iterate over each set of coordinates and check if it's within any of the polygons
for (i in 1:nrow(PRED$coords)) {
  point <- st_point(PRED$coords[i,])
  Am_overlap <- as.numeric(st_within(point, Amazon))
  Am_overlap[is.na(Am_overlap)] <- 0
  And_overlap <- as.numeric(st_within(point, Andes))
  And_overlap[is.na(And_overlap)] <- 0
  NCA_overlap <- as.numeric(st_within(point, NCA))
  NCA_overlap[is.na(NCA_overlap)] <- 0
  SCA_overlap <- as.numeric(st_within(point, SCA))
  SCA_overlap[is.na(SCA_overlap)] <- 0
  if (Am_overlap == 1) {
    PRED$region[i] <- "Amazon"
  } else if (And_overlap == 1) {
    PRED$region[i] <- "Andes"
  } else if (NCA_overlap == 1) {
    PRED$region[i] <- "NCA"
  } else if (SCA_overlap == 1) {
    PRED$region[i] <- "SCA"
  }
}

# Copy "regions" column from PRED to pred variable
pred$region <- PRED$region

# Define sampling units (su). Su are all the pixels present over Hillieae range (CA + SA south to Uruguay)
su <- rownames(na.omit(pred))

# Extract environmental variables for each su
env <- na.omit(pred[rownames(pred) %in% su, ])
colnames(env) <- c("bio1", "bio12", "bio15", "bio04", "elev", "region")
# write.table(env, file = "~/Desktop/niche_analysis/env_2.5min.txt", sep = "\t")

# Perform a PCA for environmental variables and each su
pca_su <- dudi.pca(env[, 1:5], scannf = F, nf = 2)
# write.table(pca_su$li, file = "~/Desktop/niche_analysis/niche_PCA_2.5min.txt", sep = "\t")
# manually add "region" column from "env" variable to table, with excel

# Get variables contribution to axes
inertia <- inertia.dudi(pca_su, row.inertia = TRUE,
                        col.inertia = TRUE)
var.contrib <- inertia$col.abs/100

# Get variance percentage explained by each PC
100 * pca_su$eig / sum(pca_su$eig)
# 55.264014 26.657415 11.054585  6.533337  0.490648

# Extract the occurrence data
occtot <- unique(read.csv("~/Desktop/niche_analysis/all_Hillieae_coords.csv", h = T))
names(occtot) <- c("species", "long", "lat", "syndrome", "region") # A=Andes, B=SCA, C=NCA, D=Amazon

# Extract unique pixels from each region for each species
presence_unique <- unique(data.frame(occtot$species, cellFromXY(predictors, occtot[, c(2, 3)])))

# Make a list of species with a least 5 pixels
sp_good <- occtot[occtot$species %in% as.list(names(table(presence_unique[, 1])[table(presence_unique[, 1]) >= 5])), ]
sp_good

# Make a matrix with presence/absence of each species in each pixel
# Get pixel for each occurrence of sp_good
pres_matrix <- matrix(nrow = length(rownames(env)),
                      ncol = length(unique(sp_good$species)), data = 0)
species_unique <- unique(sp_good$species)
rownames(pres_matrix) <- rownames(env)
colnames(pres_matrix) <- species_unique

for (i in 1:length(species_unique)) {
  temp_cells <- cellFromXY(predictors, sp_good[sp_good$species == species_unique[i], c(2, 3)])
  pres_matrix[rownames(pres_matrix) %in% temp_cells, i] <- 1
}

# Presence matrix with no missing data (at least one presence for each species)
p <- as.data.frame(pres_matrix[, colSums(pres_matrix) != 0])
names(p) <- colnames(p)
# write.table(p, "niche_analysis/species_presence_matrix.txt", row.names = T, col.names = T, sep = "\t")

# Get a table with, for each pixel the value of environmental PC and the presence (1) or
# Absence (0) of each species
ind_pix <- cbind(pca_su$li[, 1:2], pres_matrix[, colnames(pres_matrix) %in% as.list(colnames(pres_matrix)[colSums(pres_matrix) >= 1])])
ind_pix$region <- env$region
# write.table(ind_pix, "niche_PCA_2.5min.txt", row.names = TRUE, col.names = TRUE, sep = "\t")

# Reformat presence absence matrix; save species names in first column and positions of their pixels on Axis1 and Axis2 of the PCA in the next 2 columns
# Initialize an empty data frame to store the results
species_pca_coords <- data.frame(Species = character(),
                                 Axis1 = numeric(),
                                 Axis2 = numeric(),
                                 stringsAsFactors = FALSE)
species_pca_coords[1, ] <- colnames(species_pca_coords)

# Loop over columns of ind_pix
for (col_index in 3:27) {
  # Get the column name
  col_name <- colnames(ind_pix)[col_index]
  
  # Find rows where the value is greater than 0
  nonzero_rows <- which(ind_pix[, col_index] > 0)
  
  # Extract values from the first two columns of ind_pix for non-zero cells
  for (row_index in nonzero_rows) {
    first_value <- ind_pix[row_index, 1]
    second_value <- ind_pix[row_index, 2]
    
    # Append the results to the result_df
    species_pca_coords <- rbind(species_pca_coords, list(col_name, first_value, second_value))
  }
}

species_pca_coords <- species_pca_coords[-1, ]

axes <- c("Axis1", "Axis2")
metrics <- c("mean", "min", "max", "sd")

# Initialize an empty list to store results
niche_metrics <- list()

# Loop through each axis and compute the metrics
for (axis in axes) {
  for (metric in metrics) {
    # Create the appropriate function (mean, min, max, sd)
    func <- match.fun(metric)
    
    # Aggregate the data by species and apply the function
    result <- aggregate(x = as.numeric(species_pca_coords[[axis]]),
                        by = list(species_pca_coords$Species),
                        FUN = func)
    
    # Create a column name based on the axis and metric
    colname <- paste0(axis, "_", metric)
    colnames(result) <- c("Species", colname)
    
    # Store the result in the list
    niche_metrics[[colname]] <- result
  }
}

# Combine all results into one data frame
niche_data <- Reduce(function(x, y) merge(x, y, by = "Species"), niche_metrics)
syndrome <- read.table("species_syndrome.txt", sep = "\t", h = T)
niche_data <- merge(x = niche_data, y = syndrome, by.x = "Species", by.y = "species")
niche_data$region <- c("C", "ABD", "ABD", "BC", "B", "B", "AD", "B", "AD", "A", "B", "BC", "AD", "B", "A", "ABD", "B", "BC", "AD", "A", "BC", "B", "B", "ABD", "A")
# write.table(niche_data, "niche_identity_breadth.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
# niche_data <- read.table("niche_identity_breadth.txt", sep = "\t", h = T)

ggplot(niche_data, aes(x = Axis1_mean, y = Axis2_mean)) +
  geom_point() +
  geom_text(aes(label = Species))

# Test for significant differences between bat, HB, HM individuals' position on each axis
# https://stackoverflow.com/questions/51132259/phylogenetic-model-using-multiple-entries-for-each-species
library(mulTree)

phy <- readNexus(file = "adv.ase_with_HPD.tre")
phy$tip.label <- c("Hillia_killipii", "Hillia_pumila", "Hillia_parasitica", "Hillia_macromeris", 
                   "Hillia_macrophylla_SA", "Hillia_wurdackii", "Hillia_macrophylla_CA", "Hillia_bonoi", 
                   "Hillia_macbridei", "Hillia_allenii", "Hillia_triflora_var._pittieri", "Hillia_longifilamentosa", 
                   "Hillia_triflora_var._triflora", "Hillia_grayumii", "Hillia_foldatsii", "Hillia_illustris", 
                   "Hillia_loranthoides", "Hillia_maxonii", "Hillia_palmana", "Hillia_tetrandra", "Hillia_ulei", 
                   "Balmea_stormiae", "Cosmibuena_grandiflora", "Cosmibuena_macrocarpa", "Cosmibuena_matudae", 
                   "Cosmibuena_valerii")

species_pca_coords2 <- merge(x = species_pca_coords, y = syndrome, by.x = "Species", by.y = "species")
species_pca_coords2$Axis1 <- as.numeric(species_pca_coords2$Axis1)
species_pca_coords2$Axis2 <- as.numeric(species_pca_coords2$Axis2)
# Generate unique names by appending a sequence number
# species_pca_coords2$Taxa <- paste0(species_pca_coords2$Species, ave(species_pca_coords2$Species, species_pca_coords2$Species, FUN = seq_along))

comp_data <- as.mulTree(data = species_pca_coords2, tree = phy, taxa = "Species", clean.data = FALSE)
comp_data$data <- comp_data$data[comp_data$data$sp.col %in% phy$tip.label, ]

my_formula1 <- Axis1 ~ syndrome
my_formula2 <- Axis2 ~ syndrome

## Extracting the comparative data
mcmc_data <- comp_data$data

## MCMCglmmm
mod_mcmc1 <- MCMCglmm(fixed = my_formula1, 
                      random = ~ animal, 
                      family = "gaussian",
                      pedigree = phy, 
                      data = mcmc_data #,
                      # nitt = nitt,
                      # burnin = burnin,
                      # thin = thin,
                      # prior = prior
)
summary(mod_mcmc2)

mod_mcmc2 <- MCMCglmm(fixed = my_formula2, 
                      random = ~ animal, 
                      family = "gaussian",
                      pedigree = phy, 
                      data = mcmc_data #,
                      # nitt = nitt,
                      # burnin = burnin,
                      # thin = thin,
                      # prior = prior
)
summary(mod_mcmc2)



## ANOVA
anova_PC1_result <- aov(Axis1 ~ syndrome + Error(Species), data = species_pca_coords2)
summary(anova_PC1_result)
# TukeyHSD(anova_PC1_result)

anova_PC2_result <- aov(Axis2 ~ syndrome + Error(Species), data = species_pca_coords2)
summary(anova_PC2_result)
# TukeyHSD(anova_PC2_result)



########################################################
## Niche overlap analysis and test of niche equivalency ##
########################################################
# Script adapted from Alexandre_2017_supplementary_data/bioclimatic_niches_script.html

# Measure niche overlap along the first two niche PCA axes
niches <- read.table("~/Desktop/niche_analysis/niche_PCA_2.5min.txt", h = T, sep = "\t", check.names = F, row.names = 1)

# Only keep species present in at least 5 pixels
good <- names(niches[, 5:29])[colSums(niches[, 5:29]) >= 5]

res <- as.data.frame(matrix(ncol = 12, nrow = 0))  
names(res) <- c("species1", "species2",
                "D", "p_eqD", "CI_eqD_lower", "CI_eqD_upper",
                "p_sim1D", "CI_sim1D_lower", "CI_sim1D_upper",
                "p_sim2D", "CI_sim2D_lower", "CI_sim2D_upper")

for (spa in 1:(length(good) - 1)) {
  for (spb in (spa + 1):length(good)) {
    sp1 <- good[spa]
    sp2 <- good[spb]
    
    # Predict the scores on the axes
    region_sp1 <- unique(niches$region[niches[, sp1] == 1]) # regions where species 1 occurs
    scores.clim1 <- niches[niches$region %in% region_sp1, 3:4] # climatic data for all regions associated with the species
    
    region_sp2 <- unique(niches$region[niches[, sp2] == 1]) 
    scores.clim2 <- niches[niches$region %in% region_sp2, 3:4]
    
    scores.clim12 <- unique(rbind(scores.clim1, scores.clim2))
    
    scores.sp1 <- niches[niches[, sp1] == 1, 3:4] # climatic data for exact sites where the species occurs
    scores.sp2 <- niches[niches[, sp2] == 1, 3:4] 
    
    # Calculate occurrence density and test niche equivalency and similarity 
    R <- 100       
    z1 <- ecospat.grid.clim.dyn(scores.clim12, scores.clim1, scores.sp1, R)
    z2 <- ecospat.grid.clim.dyn(scores.clim12, scores.clim2, scores.sp2, R)
    a <- ecospat.niche.equivalency.test(z1, z2, rep = 100)
    b <- ecospat.niche.similarity.test(z1, z2, rep = 100)
    b2 <- ecospat.niche.similarity.test(z2, z1, rep = 100)
    
    # Append results
    res <- rbind(res, data.frame(
      species1 = sp1,
      species2 = sp2,
      D = a$obs$D,
      p_eqD = a$p.D,
      CI_eqD_lower = quantile(a$sim$D, probs = 0.05),
      CI_eqD_upper = quantile(a$sim$D, probs = 0.95),
      p_sim1D = b$p.D,
      CI_sim1D_lower = quantile(b$sim$D, probs = 0.05),
      CI_sim1D_upper = quantile(b$sim$D, probs = 0.95),
      p_sim2D = b2$p.D,
      CI_sim2D_lower = quantile(b2$sim$D, probs = 0.05),
      CI_sim2D_upper = quantile(b2$sim$D, probs = 0.95)
    ))
  }
}

# Save results
# write.csv(res, "results_comparison.csv", row.names = FALSE)



#################################################
## D comparisons using generalized linear models ##
#################################################
# Script adapted from Alexandre 2017
library(DHARMa)
library(glmmTMB)
library(effects)
library(phytools)
library(ape)

# All species, all regions
comp <- read.csv("results_comparison.csv", h = T)
type <- read.table("species_syndrome.txt", sep = "\t", h = T)
COMP <- merge(merge(x = comp, y = type, by.x = "species1", by.y = "species"), type, by.x = "species2", by.y = "species")
names(COMP)[13:14] <- c("syndrome_sp1", "syndrome_sp2")

pres <- read.table("species_region.txt", sep = "\t", h = T)
comp1 <- merge(merge(COMP, pres, by.x = "species1", by.y = "X"), pres, by.x = "species2", by.y = "X")
names(comp1)[15:22] <- c("Andes1", "Amazon1", "NCA1", "SCA1", 
                         "Andes2", "Amazon2", "NCA2", "SCA2")
for (i in 1:nrow(comp1)) {
  if ((comp1$Andes1[i] == 1 && comp1$Andes2[i] == 1) ||
      (comp1$Amazon1[i] == 1 && comp1$Amazon2[i] == 1) ||
      (comp1$NCA1[i] == 1 && comp1$NCA2[i] == 1) ||
      (comp1$SCA1[i] == 1 && comp1$SCA2[i] == 1)) {
    comp1$same_region[i] <- "yes"
  } else {
    comp1$same_region[i] <- "no"
  }
}

comp1$cat <- paste(comp1$syndrome_sp1, comp1$syndrome_sp2)
comp1$cat[comp1$cat == "hawkmoth bat"] <- "bat hawkmoth"
comp1$cat[comp1$cat == "hummingbird bat"] <- "bat hummingbird"
comp1$cat[comp1$cat == "hummingbird hawkmoth"] <- "hawkmoth humingbird"

comp1$same_pol_mode <- "no"
comp1$same_pol_mode[comp1$cat %in% c("bat bat", "hawkmoth hawkmoth", "hummingbrid hummingbird")] <- "yes"
comp1$cat[comp1$cat %in% list("bat hawkmoth", "bat hummingbird", "hawkmoth hummingbird")] <- "DIFF"

# Generalized linear models
# Visualize ditribution of niche overlap values 
ggplot(comp1, aes(x = D)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", alpha = 0.5) + # histogram
  geom_density(color = "red", size = 1) + # density curve
  labs(title = "Distibution of pairwise niche overlap",
       x = "D",
       y = "Density") +
  theme_minimal()

# Since a few D values are 0, and non-positive values are not allowed for the 'Beta' family,
# (and niche overlap is likely never absolutely 0) adjust the 0 values so that they are very low, but non-0
comp1$D[comp1$D == 0] <- 1e-6

# Load tree
tree <- readNexus(file = "~/Desktop/niche_analysis/adv.ase_with_HPD.tre")
tree$tip.label <- c("Hillia_macrophylla_SA", "Hillia_wurdackii", "Hillia_macromeris", "Hillia_parasitica",
                    "Hillia_pumila", "Hillia_killipii", "Hillia_bonoi", "Hillia_macbridei",
                    "Hillia_macrophylla_CA", "Hillia_allenii", "Hillia_triflora_var._pittieri", "Hillia_longifilamentosa",
                    "Hillia_triflora_var._triflora", "Hillia_grayumii", "Hillia_foldatsii", "Hillia_illustris",
                    "Hillia_loranthoides", "Hillia_maxonii", "Hillia_palmana", "Hillia_tetrandra",
                    "Hillia_ulei", "Balmea_stormiae", "Cosmibuena_grandiflora", "Cosmibuena_macrocarpa",
                    "Cosmibuena_matudae", "Cosmibuena_valerii")
tree <- drop.tip(tree, c("Hillia_bonoi", "Hillia_macbridei", "Hillia_foldatsii"))

# Calculate pairwise phylogenetic distances
phy_dist <- cophenetic.phylo(tree)
phy_cov <- vcv.phylo(tree, corr = TRUE)

# Remove non-phylo species from comp1
comp1_phy <- comp1[!(comp1$species1 %in% c("Hillia_costanensis", "Hillia_panamensis") | comp1$species2 %in% c("Hillia_costanensis", "Hillia_panamensis")), ]
comp1_phy$phylo_distance <- NA
comp1_phy$phylo_varcov <- NA

# Loop through each row in comp1 and add phylogenetic distance for species pair
for (i in 1:nrow(comp1_phy)) {
  comp1_phy$phylo_distance[i] <- phy_dist[comp1_phy$species1[i], comp1_phy$species2[i]]
}

# Loop through each row in comp1 and add the variance-covariance for species pair
for (i in 1:nrow(comp1_phy)) {
  comp1_phy$phylo_varcov[i] <- phy_cov[comp1_phy$species1[i], comp1_phy$species2[i]]
}

# GLM including phylogenetic correlation 
m1 <- (glmmTMB(D ~ 1, family = beta_family(link = "logit"), data = comp1_phy))
m2 <- (glmmTMB(D ~ same_pol_mode, family = beta_family(link = "logit"), data = comp1_phy))
m3 <- (glmmTMB(D ~ same_region, family = beta_family(link = "logit"), data = comp1_phy))
m4 <- (glmmTMB(D ~ phylo_varcov, family = beta_family(link = "logit"), data = comp1_phy))
m5 <- (glmmTMB(D ~ same_pol_mode + same_region, family = beta_family(link = "logit"), data = comp1_phy))
m6 <- (glmmTMB(D ~ same_pol_mode + phylo_varcov, family = beta_family(link = "logit"), data = comp1_phy))
m7 <- (glmmTMB(D ~ same_region + phylo_varcov, family = beta_family(link = "logit"), data = comp1_phy))
m8 <- (glmmTMB(D ~ same_pol_mode + same_region + phylo_varcov, family = beta_family(link = "logit"), data = comp1_phy))
m9 <- (glmmTMB(D ~ (same_pol_mode * phylo_varcov) + (same_region * phylo_varcov), family = beta_family(link = "logit"), data = comp1_phy))
m10 <- (glmmTMB(D ~ (same_region * same_pol_mode) + (phylo_varcov * same_pol_mode), family = beta_family(link = "logit"), data = comp1_phy)) 
m11 <- (glmmTMB(D ~ (phylo_varcov * same_region) + (same_pol_mode * same_region), family = beta_family(link = "logit"), data = comp1_phy))
m12 <- (glmmTMB(D ~ same_pol_mode * same_region + phylo_varcov, family = beta_family(link = "logit"), data = comp1_phy))
m13 <- (glmmTMB(D ~ same_pol_mode * phylo_varcov + same_region, family = beta_family(link = "logit"), data = comp1_phy))
m14 <- (glmmTMB(D ~ same_region * phylo_varcov + same_pol_mode, family = beta_family(link = "logit"), data = comp1_phy))
m15 <- (glmmTMB(D ~ same_pol_mode * same_region, family = beta_family(link = "logit"), data = comp1_phy))
m16 <- (glmmTMB(D ~ same_pol_mode * phylo_varcov, family = beta_family(link = "logit"), data = comp1_phy))
m17 <- (glmmTMB(D ~ same_region * phylo_varcov, family = beta_family(link = "logit"), data = comp1_phy))
m18 <- (glmmTMB(D ~ same_pol_mode * same_region * phylo_varcov, family = beta_family(link = "logit"), data = comp1_phy))

summary(m10)
plot(simulateResiduals(m10))
plot(allEffects(m10))

library(performance)
r2(m18) # proportion of variance in the data explained by the model
