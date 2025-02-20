# Script to estimate relaxed clock tree 
# using a partitioned GTR+gamma substition model
# and  fixed Topology
# by Isaac lichter-Marck and Mike May (https://github.com/ilichtermarck/DTE_fixed_topology/blob/main/code/phylo-relaxed-clocked-with-fixed-topology.Rev)
# modified from the original by Will Freyman
# modified again by Laymon

#converted .fasta files to .nex with:
#for file in ./*.fasta ; do
#    seqmagick convert --output-format nexus --alphabet dna $file ${file%%.fasta}.nex
#done

#######################
# Reading in the Data #
#######################
analysis = "Hillieae_DTE_fixed_topology"
# read in the alignments
filenames = v("Assembly_000000016823_Contig_1_488_CINxNeolamarckia_samps_rem.nex",
"Assembly_000000024465_Contig_1_722_CINxHamelia_samps_rem.nex",
"Assembly_000000029849_Contig_2_980_CINxNeolamarckia_paralogs_no_chimeras.NT.fs_samps_rem.nex",
"Assembly_000000063341_Contig_1_592_CINxHamelia_samps_rem.nex",
"Assembly_000000085047_Contig_1_519_CINxNeolamarckia_paralogs_no_chimeras.NT.fs_samps_rem.nex",
"Assembly_000000094771_Contig_1_679_CINxNeolamarckia_samps_rem.nex",
"Assembly_000000094993_Contig_1_1358_CINxNeolamarckia_samps_rem.nex",
"Assembly_000000098777_Contig_1_640_CINxHamelia_samps_rem.nex",
"Assembly_000000100891_Contig_1_1030_CINxNeolamarckia_paralogs_no_chimeras.NT.fs_samps_rem.nex",
"Assembly_000000102135_Contig_2_478_CINxNeolamarckia_samps_rem.nex",
"Assembly_000000107621_Contig_2_634_CINxNeolamarckia_samps_rem.nex",
"Assembly_000000107637_Contig_3_949_CINxNeolamarckia_samps_rem.nex",
"Assembly_000000114851_Contig_1_698_CINxCinchona_paralogs_no_chimeras.NT.fs_samps_rem.nex",
"Assembly_000000238591_Contig_1_642_CINxCinchona_samps_rem.nex",
"Assembly_000000249273_Contig_1_582_CINxCinchona_samps_rem.nex")

# Get the tree
observed_phylogeny = readTrees("chronos_starting_tree2.tre")[1]

# read the data
n_data_subsets = filenames.size()
for (i in 1:n_data_subsets) {
    data[i] = readDiscreteCharacterData(filenames[i])
}

# get the taxa in the tree
taxa = observed_phylogeny.taxa()

# get the names of taxa on the tree
for(i in 1:taxa.size()) {
    taxa_names[i] = taxa[i].getName()
}

# first, exclude all taxa
for(i in 1:n_data_subsets) {
    data[i].excludeTaxa(data[i].taxa())
}

# now, include taxa in the tree
for(i in 1:n_data_subsets) {
    data[i].includeTaxa(taxa_names)
}

# make sure we include missing data (add taxa that are missing in the alignment but appear in the tree)
for(i in 1:n_data_subsets) {
    data[i].addMissingTaxa(taxa_names)
}

#some useful data
n_tips <- data[1].ntaxa()
taxa <- data[1].taxa()
n_branches <- 2 * n_tips - 2
mi = 0

#set up moves and monitors vectors

moves    = VectorMoves()
monitors = VectorMonitors()

##################################
#      Substitution Model        #
#   Loop over each data subset   #
##################################

for (i in 1:n_data_subsets) {

    # exchangeability rates for partition i
    er_prior[i] <- v(1,1,1,1,1,1)
    er[i] ~ dnDirichlet(er_prior[i])
    moves.append( mvSimplexElementScale(er[i], weight=1) )

    # stationary frequencies for partition i
    pi_prior[i] <- v(1,1,1,1)
    pi[i] ~ dnDirichlet(pi_prior[i])
    moves.append( mvSimplexElementScale(pi[i], weight=1) )

    # rate matrix for partition i
    Q[i] := fnGTR(er[i],pi[i]) 

    # +Gamma for partition i
    alpha_prior <- 0.05
	alpha[i] ~ dnExponential( alpha_prior )
    gamma_rates[i] := fnDiscretizeGamma( alpha[i], alpha[i], 4, false )

    # add moves for the alpha parameter
    moves.append( mvScale(alpha[i],weight=1) )

    # the probability of a site being invariable
    pinvar[i] ~ dnBeta(1,1)
    moves.append( mvBetaProbability(pinvar[i], delta=10, tune=true, weight=2.0) )

}

##############################
# Partition rate multipliers #
##############################

# specify a rate multiplier for each partition
part_rate_mult ~ dnDirichlet( rep(1.0, n_data_subsets) )
moves.append( mvBetaSimplex(part_rate_mult, alpha=1.0, tune=true, weight=n_data_subsets) )
moves.append( mvDirichletSimplex(part_rate_mult, alpha=1.0, tune=true, weight=2.0) )

# note that we use here a vector multiplication, 
# i.e., multiplying each element of part_rate_mult by n_data_subsets
part_rate := part_rate_mult * n_data_subsets


##############
# Tree model #
##############

## if sampling is incomplete, enter the total species diversity of the ingroup clade into NUM_LINEAGES; if complete, enter n_tips
NUM_LINEAGES <- 30
# in this analysis the BDP is conditioned on the root time, which is a normal distribution 
root_time ~ dnNormal(20.9, 4.0)
root_time.setValue(observed_phylogeny.rootAge())
moves.append( mvScale(root_time, weight=2) )# the birth rate is a stochastic random variable drawn from a lognormal prior
# MCMC samples this variable using a scale proposal
speciation_mean <- ln( abs(ln(abs(NUM_LINEAGES)/2.0)) / 15.6 )
speciation_sd <- 0.587405
speciation ~ dnLognormal(mean=speciation_mean, sd=speciation_sd) 
moves.append( mvScale(speciation, lambda=1.0, tune=true, weight=3.0) )
extinction_mean <- ln( abs(ln(abs(NUM_LINEAGES)/2.0)) / 15.6 )
extinction_sd <- 0.587405
extinction ~ dnLognormal(mean=extinction_mean, sd=extinction_sd) 
moves.append( mvScale(extinction, lambda=1.0, tune=true, weight=3.0) )
diversification := speciation - extinction
# rho is the probability of sampling species at the present
rho := Probability(n_tips/NUM_LINEAGES)
# the time tree is a stochastic node modeled by the constant rate birth-death process (dnBDP)

bdp = dnBDP(lambda=speciation, mu=extinction, rho=rho, rootAge=abs(root_time), taxa=taxa)
phylogeny ~ bdp
phylogeny.setValue(observed_phylogeny)
moves.append( mvSubtreeScale(phylogeny, weight=5.0) )
moves.append( mvNodeTimeSlideUniform(phylogeny, weight=15.0) )
moves.append( mvNodeTimeScale(phylogeny, weight=15.0) )


###################
# Molecular clock #
###################
clock_mean ~ dnLoguniform(1e-6,1)
clock_mean.setValue(0.05)
moves.append( mvScale(clock_mean, weight=5.0) )
clock_sd ~ dnExponential(abs(1 / 0.587405))
moves.append( mvScale(clock_sd, weight=5.0) )
# use a discretized lognormal

for(i in 1:n_branches) {

    # draw the branch rate from a lognormal
    branch_rates[i] ~ dnLognormal( ln(clock_mean) - clock_sd * clock_sd * 0.5, clock_sd)
    moves.append( mvScale(branch_rates[i], weight=1.0) )

}
mean_rt := mean(branch_rates)
# some joint moves
speciation.setValue(0.5)
extinction.setValue(0.15)
up_down_scale_div = mvUpDownScale(lambda=1.0, weight=20)
up_down_scale_div.addVariable( speciation, TRUE )
up_down_scale_div.addVariable( extinction, TRUE )
moves.append( up_down_scale_div )
up_down_rate_scale = mvUpDownScale(lambda=1.0, weight=10)
up_down_rate_scale.addVariable( clock_mean,   TRUE )
up_down_rate_scale.addVariable( branch_rates, TRUE )
moves.append( up_down_rate_scale )
up_down_scale_tree = mvUpDownScale(lambda=1.0, weight=20)
up_down_scale_tree.addVariable( phylogeny,    TRUE )
up_down_scale_tree.addVariable( clock_mean,   FALSE )
up_down_scale_tree.addVariable( branch_rates, FALSE )
moves.append( up_down_scale_tree )
rate_age_proposal = mvRateAgeProposal(phylogeny, weigh=20)
rate_age_proposal.addRates( branch_rates )
moves.append( rate_age_proposal )



###################
###################
# PhyloCTMC Model #
###################

for (i in 1:n_data_subsets) {
    phyloSeq[i] ~ dnPhyloCTMC(tree=phylogeny, Q=Q[i], branchRates=part_rate[i] * branch_rates, siteRates=gamma_rates[i], pInv=pinvar[i], type="DNA")
    phyloSeq[i].clamp(data[i])
}

#############
# THE Model #
#############
mymodel = model(Q)
monitors.append( mnModel(filename="output/" + analysis + ".log",printgen=10, separator = TAB) )
monitors.append( mnFile(filename="output/" + analysis + ".trees",printgen=10, separator = TAB, phylogeny) )
monitors.append( mnScreen(printgen=10, root_time) )
mymcmc = mcmc(mymodel, monitors, moves)
mymcmc.run(generations=300000)



# Now, we will analyze the tree output.
# Let us start by reading in the tree trace
treetrace = readTreeTrace("output/" + analysis + ".trees", treetype="clock", burnin=0.25)
# and get the summary of the tree trace
#treetrace.summarize()
map_tree = mapTree(treetrace,"output/" + analysis + ".tree")
# you may want to quit RevBayes now
q()