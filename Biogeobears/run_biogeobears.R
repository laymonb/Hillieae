library(ape)
library(optimx) 
library(GenSA)
library(rexpokit)
library(cladoRcpp)
library(snow)    
library(parallel)
library(BioGeoBEARS)

wd <- "/Users/laymonball/Library/Mobile Documents/com~apple~CloudDocs/Manuscripts/Hillieae_phylogenetics_Ch1_paper/Biogeography/biogeobears"
setwd(wd)


############################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE #
############################################
trfn = "Hillieae_time_tree.tre"
tr = read.tree(trfn)
plot(tr, cex=0.5)
axisPhylo()


geogfn = "Hillieae_ranges.phy"
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
max_range_size <- max(rowSums(dfnums_to_numeric(tipranges@df)))

# Get your states list
areas = getareas_from_tipranges_object(tipranges)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

ranges_list = NULL
for (i in 1:length(states_list_0based))
{    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

ranges_list

# Manually modify list of allowed ranged (exclude disjunct ranges)
nonadjacent=c("AC", "BD", "CD", "ACD", "BCD")
keepTF = ranges_list %in% nonadjacent == FALSE

ranges_list_NEW = ranges_list[keepTF]
length(ranges_list_NEW) 

states_list_0based_NEW = states_list_0based[keepTF]
length(states_list_0based_NEW)     



###########
# Run DEC #
###########

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE  
BioGeoBEARS_run_object$distsfn = "distances.txt"
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0.0

# Update list of allowed ranges
BioGeoBEARS_run_object$states_list = states_list_0based_NEW

# Run this to check inputs. Read the error messages if you get them!
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)


runslow = TRUE
resfn = "Hillieae_DEC.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}



# Plot
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
plot_BioGeoBEARS_results(results_object)

pdf("DEC.pdf", height=6, width=6)
plot_BioGeoBEARS_results(results_object)
dev.off()



###############
# Run DIVALIKE#
###############
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$distsfn = "distances.txt"
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0.0

# Update list of allowed ranges
BioGeoBEARS_run_object$states_list = states_list_0based_NEW

# Run this to check inputs. Read the error messages if you get them!
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Hillieae_DIVALIKE.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}



# Plot
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
plot_BioGeoBEARS_results(results_object)

pdf("DIVALIKE.pdf", height=6, width=6)
plot_BioGeoBEARS_results(results_object)
dev.off()



###################
# Run BAYAREALIKE #
###################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE
BioGeoBEARS_run_object$distsfn = "distances.txt"
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = TRUE
BioGeoBEARS_run_object$num_cores_to_use = 1
BioGeoBEARS_run_object$force_sparse = FALSE
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0.0

# Update list of allowed ranges
BioGeoBEARS_run_object$states_list = states_list_0based_NEW

# Check the inputs
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

runslow = TRUE
resfn = "Hillieae_BAYAREALIKE.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resBAYAREALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resBAYAREALIKE = res
}



# Plot 
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
plot_BioGeoBEARS_results(results_object)

pdf("BAYAREALIKE.pdf", height=6, width=6)
plot_BioGeoBEARS_results(results_object)
dev.off()



######################################################
# Statistics -- Comparing DEC, DIVA, and BAYAREALIKE #
######################################################
aic_DEC <- getAIC(get_LnL_from_BioGeoBEARS_results_object(resDEC), 2)
aic_DIVA <- getAIC(get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE), 2)
aic_BAYAREA <- getAIC(get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE), 2)



###############################################
# Probabilities of states/ranges at each node #
###############################################
# To get the probabilities of each state/range at each node:
# What you want, if "res" is your results object, is:
resDIVALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node

# In this table:
# - columns are states/ranges
# - rows are nodes, in APE order (tips, then root, then internal)

#  You can see the node numbers in the same APE order with:
trtable = prt(tr, printflag=FALSE)
head(trtable)
tail(trtable)

# You can plot APE node labels with:
plot(tr)
axisPhylo()
nodelabels()
tiplabels(1:length(tr$tip.label))

# Make the node numbers the row names
# Make the range_list the column names
range_probabilities = as.data.frame(resDIVALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node)
row.names(range_probabilities) = trtable$node
names(range_probabilities) = ranges_list_NEW

# Look at the table (first six rows)
head(range_probabilities)

# Write the table to a tab-delimited text file (for Excel etc.)
write.table(range_probabilities, file="range_probabilities.txt", quote=FALSE, sep="\t")

















