#Use this script to calculate the proportion of gene trees in which pollination syndromes are monophyletic

#Load the required library
library(ape)

#Load trees
gene_trees <- read.tree("../all_rooted_gene_trees.tre") #gene trees using the pxrr function from phyx
class(gene_trees) <- "multiPhylo"

#Drop outgroup tips from trees
drop_tips <- c("Pal_attenuata_Fuentes4647", "Psy_panamensis_Stevens32285", "COU_hexandra_Callejas4587", "HOF_phoenicopoda_Wendt3392", "HAM_patens_LagomarsinoSN") 

pruned.trees <- lapply(gene_trees, drop.tip, tip = c(drop_tips))
class(pruned.trees) <- "multiPhylo"

#Load the data
syndrome_df <- read.csv("syndromes_all_samps.csv")

#Specify the flower colors (syndromes) you want to check for monophyly
colors_to_check <- c("B")

########################
# CALCULATE PROPORTION #
########################

#Initialize a variable to count the number of gene trees where specified colors are monophyletic
count_monophyletic <- 0

#Loop through each gene tree
for (i in 1:length(pruned.trees)) {
  gene_tree <- pruned.trees[[i]]
  
  #Loop through the specified colors
  for (color in colors_to_check) {
    #Get species with the specified color
    species_with_color <- syndrome_df$taxon[syndrome_df$syndrome == color]
    
    #Check if they are monophyletic in the current gene tree
    if (is.monophyletic(gene_tree, species_with_color, reroot = F)) {
      count_monophyletic <- count_monophyletic + 1
      break
    }
  }
}

#Calculate the proportion of gene trees where specified colors are monophyletic
proportion_monophyletic <- count_monophyletic / length(pruned.trees)

#Print the proportion
cat("Proportion of gene trees where specified colors are monophyletic:", proportion_monophyletic, "\n")


##############################
# EXTRACT MONOPHYLETIC TREES #
##############################

#Initialize a list to store gene trees where specified colors are monophyletic
monophyletic_trees <- list()

#Loop through each gene tree
for (i in 1:length(pruned.trees)) {
  gene_tree <- pruned.trees[[i]]
  
  #Loop through the specified colors
  for (color in colors_to_check) {
    #Get species with the specified color
    species_with_color <- syndrome_df$taxon[syndrome_df$syndrome == color]
    
    #Check if they are monophyletic in the current gene tree
    if (is.monophyletic(gene_tree, species_with_color, reroot = F)) {
      monophyletic_trees[[paste("Tree", i)]] <- gene_tree
      break
    }
  }
}

#Print the names of extracted gene trees
cat("Names of gene trees where specified colors are monophyletic:\n")
print(names(monophyletic_trees))
#plot(monophyletic_trees$`Tree 1320`)
