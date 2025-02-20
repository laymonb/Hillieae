setwd("/Users/laymonball/Library/Mobile Documents/com~apple~CloudDocs/Manuscripts/Hillieae_phylogenetics_Ch1_paper/Morphology/Revbayes")

library(RevGadgets)
library(ggplot2)
library(gridExtra)
library(ggtree)
library(dplyr)
library(phytools)
library(IRanges)
library(plotrix)
library(treeio)



#tree_file = "output/syndromes_ase_ERM.tree"

# process the ancestral states
#ase <- processAncStates(tree_file,
                        # Specify state labels.
                        # These numbers correspond to
                        # your input data file.
                        #state_labels = c("0" = "Hawkmoth", "1" = "Hummingbird", "2" = "Bat"))

#ase_data <- as.data.frame(ase@data)
#write.table(ase_data, "ase_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#to add HPD bars to plot, need to combine t_beast@data$height_0.95_HPD with ase@data
#time_tree <- read.beast("data/Hillieae_DTE_fixed_topology.tree")
#time_tree_data <- time_tree@data[c("age_0.95_HPD", "node")] #select 0.95 HPD and node columns from t_beast@data3time_tree_data <- time_tree_data %>% arrange(as.numeric(time_tree_data$node)) #reorder rows according to node column
#time_tree_data <- time_tree_data %>% arrange(as.numeric(time_tree_data$node)) #reorder rows according to node column

#add index numbers to the tibble; these are given in the Revbayes output ase tree; 
#there is one per node, but these are not in the same order as the node numbers that you get when you plot node labels with ggtree (here, the node numbers in the ase and t_beast tree are not the same either).
#to figure out which index number corresponds to each node number, I plotted both with ggtree, compared them and ase_data$index, and then added the index numbers to t_beast_data manually:
#ggtree(ase@phylo) + geom_text(aes(label = node))
#ggtree(time_tree) + geom_text(aes(label = index))
#time_tree_data$index <- c(17, 18, 16, 21, 20, 19, 15, 14, 13, 10, 9, 11, 12, 8, 6, 7, 25, 26, 23, 24, 22, 5, 4, 3, 2, 1, 51, 47, 46, 41, 34, 33, 31, 28, 27, 30, 29, 32, 40, 38, 37, 35, 36, 39, 45, 44, 42, 43, 50, 49, 48)
#then use combine_annotations.R script to combine the ase and t_beast tree annotations. I manually reformatted the new tree (adv.ase_with_HPD.tre) to NEXUS format so that it could be read by RevGadgets.

#now update tree_fn variable so inclue tree with both sets of annotations
tree_file = "output/syndromes_ase_ERM_with_HPD.tree"

# process the ancestral states
ase <- processAncStates(tree_file,
                        # Specify state labels.
                        # These numbers correspond to
                        # your input data file.
                        state_labels = c("0" = "Hawkmoth", "1" = "Hummingbird", "2" = "Bat"))
ase@data <- ase@data %>% arrange(as.numeric(ase@data$node)) #reorder rows according to node column
#df <- apply(data.frame(ase@data), 2, as.character)
#write.table(df, "ase_data_w_HPD.txt", sep = "\t", row.names = F)

#To plot posterior probabilities from stochastic character mapping  
#Load and summarize simmap trees
sim = describe.simmap(read.simmap(file = "simmaps.tree", format = "phylip"))
plot(sim)

#Save data table with posterior probabilities
#new_pp <- data.frame(t(apply(t(sim$ace), 2, function(col) sort(col, decreasing = TRUE))))[c(1:25), ]
#colnames(new_pp) <- c("anc_state_1_pp", "anc_state_2_pp", "anc_state_3_pp")
#write.table(new_pp, "stoch_mapping_pp.txt", sep = "\t", row.names = T)
new_pp <- read.table("./stoch_mapping_pp.txt")


#Update ase@data 
ase@data$anc_state_1_pp[27:51] <- new_pp$anc_state_1_pp[1:25]
ase@data$anc_state_2_pp[27:51] <- new_pp$anc_state_2_pp[1:25]
ase@data$anc_state_3_pp[27:51] <- new_pp$anc_state_3_pp[1:25]



ase@phylo$tip.label <- c("H. macrophylla Col", "H. wurdackii", "H. macromeris", "H. parasitica",
                         "H. pumila", "H. killipii", "H. bonoi", "H. macbridei",
                         "H. macrophylla CR", "H. allenii", "H. triflora var. pittieri", "H. longifilamentosa",
                         "H. triflora var. triflora", "H. grayumii", "H. foldatsii", "H. illustris",
                         "H. loranthoides", "H. maxonii", "H. palmana", "H. tetrandra",
                         "H. ulei", "B. stormiae", "C. grandiflora", "C. macrocarpa",
                         "C. matudae", "C. valerii")

colors <- c("#999933", "#888888", "#882255")

pp_map2=plotAncStatesMAP2(t = ase, 
                          node_color = colors, node_size = c(3, 9), #branch.color = "white",
                          tip_labels_size = 6, #I edited the original plotAncStatesMAP function so that it could take an age_max parameter, allowing the 95% HPD intervals to be visible on the plot. This is saved as a separate R script that just needs to be ran prior to using plotAncStatesMAP2.
                          age_max = 27,
                          timeline = T, geo = F, time_bars = F,
                          state_transparency = 1,
                          tip_labels_italics = T,
                          #time_bars = F,
                          node_labels_as = NULL, node_labels_centered = T,
                          node_labels_offset = 0,
#                          tip_labels_states_size = 6, tip_states = T,
#                          tip_states_shape = c(19, 19, 19, 19, 19, 19, 19, 17, 19, 19, 19, 19, 17, 19, 15, 15, 15, 15, 17, 19, 19, 19, 19, 19, 15, 17),
                          tip_states_size = 4,
                          cladogenetic = F, tip_labels_offset = 0.5, size = 1.1) + #cladogenetic = T to show probabilites/states at branch corners
  ggplot2::theme(legend.position = c(0.08, 0.67),
                 legend.key = element_blank(), legend.title = element_blank(),
                 legend.text = element_text(size = 20),
                 legend.background = element_rect(fill="transparent"),
                 legend.key.size = unit(0.75, "cm"))

pp_map2

#set the color for each 95% HPD interval
#colors2 <- c("#888888", "#888888", "#888888", "#888888", "#888888", "#888888", "#888888", "#888888", "#888888", "#888888", "#888888", "#888888", 
 #            "#999933", "#999933", "#882255", "#882255", "#882255", "#999933", "#888888", "#888888", "#888888", "#888888", "#888888", "#888888",
  #           "#888888") #colors for 95% HPD bars

#plot2 <- pp_map2 + geom_range(range = 'height_0.95_HPD', color = colors2,  alpha = .7, size = 2.3)# + geom_text(aes(label = node))
#plot2


#save
#ggsave(file=plot_fn, plot=pp_map, device="pdf", height=15, width=12, useDingbats=F)
ggsave(file="Fig2b.pdf", plot=pp_map2, device="pdf", height=14, width=16, useDingbats=F)



#Plot simmap 
sim2 <- read.simmap("output/character.tree", format = "phylip")

colors = vector()
for (i in 1:length( sim2$maps ) ) { 
  colors = c(colors, names(sim2$maps[[i]]) )
}
colors = sort(as.numeric(unique(colors)))
cols = setNames( c("#888888", "#882255", "#999933"), colors)

pdf(file = "./syndrome_simmap.pdf", width = 16, height = 14, bg = "transparent")
#par(mar = c(8, 8, 1.5, 1.5), mgp = c(4.5, 1.5, 0))
#par(cex.axis = 2.5, cex.lab = 3) 
plotSimmap(ladderize.simmap(sim2, F), cols, fsize = 0, lwd = 4, split.vertical = TRUE, ftype = "i", xlim = c(0, 28))

dev.off()



#Plot density plots for # of shifts
df <- data.frame(sim$count)[, c(2:7)]
colnames(df) <- c("HM2HB", "HM2Bat", "HB2HM", "HB2Bat", "Bat2HM", "Bat2HB")

plot_hist <- function(column, name) {
  ggplot(data = data.frame(Value = column), aes(x = Value)) +
    geom_histogram(color = "black", fill = "gray90", bins = 10) +
    labs(title = name, x = "n Shifs", y = "") +
    theme_minimal()
}

plot_list <- list(
  plot_hist(df$HM2HB, "HM to HB"),
  plot_hist(df$HM2Bat, "HM to Bat"),
  plot_hist(df$HB2HM, "HB to HM"),
  plot_hist(df$HB2Bat, "HB to Bat"),
  plot_hist(df$Bat2HM, "Bat to HM"),
  plot_hist(df$Bat2HB, "Bat to HB")
)

grid.arrange(
  grobs = plot_list,
  nrow = 2,
  ncol = 3
)

#Make table with mean and 95% CI
results <- apply(df, 2, function(x) {
  c(mean = median(x), lower_CI = quantile(x, 0.025), upper_CI = quantile(x, 0.975))
})

results <- t(results)  # Transpose for clarity
colnames(results) <- c("Median","Lower CI", "Upper CI")
print(results)



# Test for diifferent evolutionary rates
tree <- readNexus(file = "~/Desktop/niche_analysis/adv.ase_with_HPD.tre")

fit.single <- fitMk(tree, sim$tree[[1]])




