library(phytools)
library(ape)
setwd("/Users/laymonball/Library/Mobile Documents/com~apple~CloudDocs/Manuscripts/Hillieae_phylogenetics_Ch1_paper")

tree <- read.nexus(file = "Hillieae_DTE_fixed_topology.tree")

x <- read.csv("Morphology/syndromes_all_samps.csv", h = TRUE, row.names = 1)
synds <- setNames(x$syndrome, rownames(x))
#dotTree(tree, synds)

fitER <- fitMk(tree, synds, model = "ER")
fitSYM <- fitMk(tree, synds, model = "SYM")
fitARD <- fitMk(tree, synds, model = "ARD")

hrmER1 <- fitHRM(tree, synds, ncat = 2, model = "ER", umbral = TRUE, niter = 10, pi = "fitzjohn", opt.method = "nlimnb")
hrmER2 <- fitHRM(tree, synds, ncat = 3, model = "ER", umbral = TRUE, niter = 10, pi = "fitzjohn", opt.method = "nlimnb")
hrmSYM1 <- fitHRM(tree, synds, ncat = 2, model = "SYM", umbral = TRUE, niter = 10, pi = "fitzjohn", opt.method = "nlimnb")
hrmSYM2 <- fitHRM(tree, synds, ncat = 3, model = "SYM", umbral = TRUE, niter = 10, pi = "fitzjohn", opt.method = "nlimnb")
hrmARD1 <- fitHRM(tree, synds, ncat = 2, model = "ARD", umbral = TRUE, niter = 10, pi = "fitzjohn", opt.method = "nlimnb")
hrmARD2 <- fitHRM(tree, synds, ncat = 3, model = "ARD", umbral = TRUE, niter = 10, pi = "fitzjohn", opt.method = "nlimnb")

gamER <- fitgammaMk(tree, synds, model = "ER")
gamSYM <- fitgammaMk(tree, synds, model = "SYM")
gamARD <- fitgammaMk(tree, synds, model = "ARD")


AIC <- setNames(sapply(list(fitER, fitSYM, fitARD, hrmER1, hrmER2, hrmSYM1, hrmSYM2, hrmARD1, hrmARD2, gamER, gamSYM, gamARD), AIC), 
                c("ER", "SYM", "ARD", "hrmER1", "hrmER2", "hrmSYM1", "hrmSYM2", "hrmARD1", "hrmARD2", "gamER", "gamSYM", "gamARD"))
AIC
aic.w(AIC)
