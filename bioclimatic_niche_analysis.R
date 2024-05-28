#Hillieae bioclimatic niche description
#modified from Alexandre 2017 Supplemental Material bioclimatic_niches_script.html
setwd("~/Desktop")

library(dismo)
library(raster)
library(ade4)


####First need to test for correlations between 20 variables####
#Load the environmental variables (resolution: 2.5 arc-minute)
files <- list.files(path = paste("niche_analysis/wc2.1_2.5m_bio/Neotropics_clipped_asc/", sep = ""), pattern = "asc$", full.names = TRUE)
files

predictors <- stack(files)
#predictors

names(predictors) <- c("bio01","bio10","bio11","bio12","bio13","bio14","bio15","bio16",
                       "bio17","bio18","bio19","bio02","bio03","bio04","bio05","bio06","bio07","bio08","bio09",
                       "elev")



####assess rasters####
#calculate correlation coefficient between rasters, and then plot; can decide whether to remove highly correlated variables
pears_cor_plots <- raster.cor.plot(predictors, method = "pearson")
#pears_cor_plots
#pears_cor_plots$cor.heatmap$data

#write.table(pears_cor_plots$cor.heatmap$data, "niche_analysis/wc2.1_2.5m_bio/clipped_asc/pears_cor_plot_heatmapVals.txt", quote = FALSE, sep = "\t")

vars_collin_tree <- removeCollinearity(
  predictors,
  multicollinearity.cutoff = 0.6,
  select.variables = FALSE,
  sample.points = FALSE,
  plot = TRUE,
  method = "pearson"
)

#based on above outputs, I'm going to keep bio1, bio4, bio12, bio15, and elevation.
vars_filt <- subset(predictors, c(1, 4, 7, 14, 20))
#vars_filt

vars_filt_dframe <- as.data.frame(subset(predictors, c(1, 4, 7, 14, 20)))



####Following methods from Alexandre 2017####
#Extract the environmental variables (resolution: 2.5 arc-minute)
ras_files <- list.files(path=paste("~/Desktop/niche_analysis/wc2.1_2.5m_bio/filtered_clipped/", sep=""),pattern="sdat$",full.names=TRUE)
predictors <- stack(ras_files)
predictors
names(predictors) <- c("bio01","bio12","bio15", "bio04","elev")
pred <- as.data.frame(predictors)
PRED <- cbind(pred,rownames(pred))

#Load region polygon shapefiles
Amazon <- st_read("~/Desktop/niche_analysis/region_polygon_shapefiles/Amazon_dissolved.shp")
Andes <- st_read("~/Desktop/niche_analysis/region_polygon_shapefiles/Andes_dissolved.shp")
NCA <- st_read("~/Desktop/niche_analysis/region_polygon_shapefiles/NCA_dissolved.shp")
SCA <- st_read("~/Desktop/niche_analysis/region_polygon_shapefiles/SCA_dissolved.shp")

#Get lat/long coordinates from raster cells
PRED$coords <- xyFromCell(predictors, as.numeric(PRED$`rownames(pred)`))

PRED$region <- NA

#Iterate over each set of coordinates and check if it's within any of the polygons; if it is, assign the corresponding label (this takes a while)
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
    if (Am_overlap ==  1 ) {
      PRED$region[i] <- "Amazon"
    } else if (And_overlap == 1) {
      PRED$region[i] <- "Andes"
    } else if (NCA_overlap == 1) {
      PRED$region[i] <- "NCA"
    } else if (SCA_overlap == 1) {
      PRED$region[i] <- "SCA"
    }
  }

#Copy "regions" column from PRED to pred variable
pred$region <- PRED$region

#Define sampling units (su) SU are all the pixel present over the carribean
su <- rownames(na.omit(pred))

# extract environmental variables for each su
env <- na.omit(pred[rownames(pred) %in% su,])
colnames(env) <- c("bio1", "bio12", "bio15", "bio04", "elev", "region")



#Perform a PCA for environmental variables and each su
pca_su <- dudi.pca(env[, 1:5],scannf = F, nf = 2)

#Get variables contribution to axes
inertia <- inertia.dudi(pca_su, row.inertia = TRUE,
                        col.inertia = TRUE)
var.contrib <- inertia$col.abs/100

#Get variance percentage explained by each PC
100*pca_su$eig/sum(pca_su$eig)


####Project species onto the environmental PCA####
#Extract the occurrence data
occtot <- unique(read.csv("~/Desktop/niche_analysis/phylo_Hillieae_coords.csv",h=T))
names(occtot) <- c("species","long","lat","syndrome","region") #A=Andes, B=SCA, C=NCA, D=Amazon

#Extract unique pixels from each region for each species
presence_unique <- unique(data.frame(occtot$species, cellFromXY(predictors, occtot[,c(2,3)])))

#Make a list of species with a least 1 pixels (or, can change the minimum number of points per species)
sp_good <- occtot[occtot$species%in%as.list(names(table(presence_unique[,1])[table(presence_unique[,1])>=1])),]
sp_good

#Make a matrix with presence/absence of each species in each pixel
#Get pixel for each occurrence of sp_good
pres_matrix <- matrix(nrow = length(rownames(env)),
                    ncol = length(unique(sp_good$species)), data = 0)
species_unique <- unique(sp_good$species)
rownames(pres_matrix) <- rownames(env)
colnames(pres_matrix) <- species_unique

for (i in 1:length(species_unique)) {
  temp_cells <- cellFromXY(predictors,sp_good[sp_good$species==species_unique[i],c(2,3)])
  pres_matrix[rownames(pres_matrix) %in% temp_cells,i] <- 1
}

#Get presence matrix with no missing data
p <- as.data.frame(pres_matrix[,colSums(pres_matrix)!=0])
names(p) <- colnames(p)

#Get a table with, for each pixel, the value of environmental PC and the presence (1) or
# absence (0) of each species
ind_pix <- cbind(pca_su$li[,1:2],pres_matrix[,colnames(pres_matrix)%in%as.list(colnames(pres_matrix)[colSums(pres_matrix)>=1])])

#Reformat presence absence matrix; save species names in first column and positions of their pixels on Axis1 and Axis2 of the PCA in the next 2 columns
# Initialize an empty data frame to store the results
species_pca_coords <- data.frame(Species = character(),
                        Axis1 = numeric(),
                        Axis2 = numeric(),
                        stringsAsFactors = FALSE)
species_pca_coords[1, ] <- colnames(species_pca_coords)

#Loop over columns of ind_pix
for (col_index in 3:28) {
  #Get the column name
  col_name <- colnames(ind_pix)[col_index]
  
  #Find rows where the value is greater than 0
  nonzero_rows <- which(ind_pix[, col_index] > 0)
  
  #Extract values from the first two columns of ind_pix for non-zero cells
  for (row_index in nonzero_rows) {
    first_value <- ind_pix[row_index, 1]
    second_value <- ind_pix[row_index, 2]
    
    #Append the results to the result_df
    species_pca_coords <- rbind(species_pca_coords, list(col_name, first_value, second_value))
  }
}

species_pca_coords <- species_pca_coords[-1,]

#For each species, calculate mean Axis1 and Axis2 positions and standard deviation; calculate mean and std by group
niche_identity_1 <- aggregate(x= as.numeric(species_pca_coords$Axis1),
                       # Specify group indicator
                       by = list(species_pca_coords$Species),      
                       # Specify function (i.e. mean)
                       FUN = mean)
colnames(niche_identity_1) <- c("Species", "Axis1_id")

niche_identity_2 <- aggregate(x= as.numeric(species_pca_coords$Axis2),
                              # Specify group indicator
                              by = list(species_pca_coords$Species),      
                              # Specify function (i.e. mean)
                              FUN = mean)
colnames(niche_identity_2) <- c("Species", "Axis2_id")

niche_breadth_1 <- aggregate(x= as.numeric(species_pca_coords$Axis1),
                             # Specify group indicator
                             by = list(species_pca_coords$Species),      
                             # Specify function (i.e. mean)
                             FUN = sd)
colnames(niche_breadth_1) <- c("Species", "Axis1_breadth")

niche_breadth_2 <- aggregate(x= as.numeric(species_pca_coords$Axis1),
                             # Specify group indicator
                             by = list(species_pca_coords$Species),      
                             # Specify function (i.e. mean)
                             FUN = sd)
colnames(niche_breadth_2) <- c("Species", "Axis2_breadth")

niche_id_br <- merge(niche_identity_1, niche_identity_2, by = intersect(names(niche_identity_1), names(niche_identity_2)))
niche_id_br <- merge(niche_id_br, niche_breadth_1, by = intersect(names(niche_id_br), names(niche_breadth_1)))
niche_id_br <- merge(niche_id_br, niche_breadth_2, by = intersect(names(niche_id_br), names(niche_breadth_2)))
#write.table(niche_id_br, "niche_identity_breadth.txt", row.names = TRUE, col.names = TRUE, sep = "\t")


