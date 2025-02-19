setwd("/Users/laymonball/Library/Mobile Documents/com~apple~CloudDocs/Manuscripts/Hillieae_phylogenetics_Ch1_paper")

library(l1ou)
library(ape)
library(treeio)
library(phangorn)
library(phytools)

tree <- read.nexus("Hillieae_DTE_fixed_topology.tree")
tree$tip.label <- c("Hillia_macrophylla_SA", "Hillia_wurdackii", "Hillia_macromeris", "Hillia_parasitica",
                    "Hillia_pumila", "Hillia_killipii", "Hillia_bonoi", "Hillia_macbridei",
                    "Hillia_macrophylla_CA", "Hillia_allenii", "Hillia_triflora_var._pittieri", "Hillia_longifilamentosa",
                    "Hillia_triflora_var._triflora", "Hillia_grayumii", "Hillia_foldatsii", "Hillia_illustris",
                    "Hillia_loranthoides", "Hillia_maxonii", "Hillia_palmana", "Hillia_tetrandra",
                    "Hillia_ulei", "Balmea_stormiae", "Cosmibuena_grandiflora", "Cosmibuena_macrocarpa",
                    "Cosmibuena_matudae", "Cosmibuena_valerii")
plot(tree)
is.ultrametric(tree)

##########################################################
## fix ultrametric tree:                                ##
## http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html ##
##########################################################
force.ultrametric <- function(tree, method = c("nnls", "extend")) {
  method <- method[1]
  if (method == "nnls") 
    tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE, trace = 0)
  else if (method == "extend") {
    h <- diag(vcv(tree))
    d <- max(h) - h
    ii <- sapply(1:Ntip(tree), function(x, y) which(y == x),
                 y = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

tree_ult <- ladderize(force.ultrametric(tree, method = "extend"), F)
plot(tree_ult)
is.ultrametric(tree_ult)

##########################################################
## Run l1ou using PC axes                               ##
##########################################################
niche_dat <- read.table("niche_identity_breadth.txt", header = T, sep = "\t")
rownames(niche_dat) <- niche_dat[, 1]
niche_dat[1]

hil <- adjust_data(tree_ult, as.matrix(niche_dat[, 2:3]))
colnames(hil$Y) <- c("PC1", "PC2")
eModel <- estimate_shift_configuration(hil$tree, hil$Y, alpha.upper = 5, criterion = "AICc") # raised alpha.upper to 5 to help the models converge
plot(eModel)

# png(filename = "l1ou.png", width = 1000, height = 1000)
plot(eModel, palette = c("#DDCC77", "#88CCEE", "gray90"), asterisk = F, edge.shift.ann = F, 
     cex = 1.5, edge.width = 4, label.offset = 0.05, bar.axis = T)
# dev.off()

result <- l1ou_bootstrap_support(eModel, nItrs = 1000, multicore = TRUE, nCores = 4)
result$detection.rate

##########################################################
## Run l1ou using separate variables                    ##
##########################################################
niche_dat <- read.csv("~/Desktop/niche_analysis/median_climate_data.csv", header = T, row.names = 1)

hil <- adjust_data(tree_ult, as.matrix(niche_dat))
colnames(hil$Y) <- c("bio1", "bio12", "bio15", "bio4", "elev")
eModel <- estimate_shift_configuration(hil$tree, hil$Y, criterion = "pBIC", alpha.upper = 8)
plot(eModel, edge.shift.ann = F)
edgelabels(seq_along(tree$edge[, 1]), adj = c(0.5, -0.5), frame = "n", cex = 0.8, col = "blue")

png(filename = "l1ou.png", width = 1800, height = 1000)
plot(eModel, palette = c("#88CCEE", alpha("#88CCEE", 0.35), alpha("#88CCEE", 0.35), "#DDCC77", "gray90"), 
     asterisk = T, edge.shift.ann = F, cex = 0.75, edge.width = 5.5, label.offset = 0.2, bar.axis = F, 
     x.lim = c(0, 1))
dev.off()

result_pBIC <- l1ou_bootstrap_support(eModel, nItrs = 1000, multicore = TRUE, nCores = 4)
result_pBIC$detection.rate
