#################################################################################################
# AIM 5: uniform distribution of null p-values
#################################################################################################
## Deps
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("./") # <=
@everywhere using CpelNano
using StatsBase
using StatsPlots
using Distributions
using Combinatorics

## Source code
root_dir = dirname(dirname(pathof(CpelNano)))
include("$(root_dir)/calls/PaperSimulations/Aim-5/Aim-5-Source.jl")

## Constants

# IO
const aim_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-5/" # <=
const blind_friend_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]

# Set up
@everywhere const M = 5
@everywhere const POBS = 1.0
@everywhere const LEVEL_STD_M = 2.0
@everywhere const LEVEL_STD_U = 2.0
@everywhere const LEVEL_MEAN_M = 84.0
@everywhere const LEVEL_MEAN_U = 80.0
@everywhere const EF_PERC_RNG = 0.1:0.1:0.1
@everywhere const LMAX_TWO_SAMP_AIM_5 = 20

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#################################################################################################
# Calls
#################################################################################################

# Args
n_sim_reps = parse(Int64, ARGS[1])

# Generate null p-values
pvals_tmml, pvals_tnme, pvals_tcmd = run_null_sim(n_sim_reps)

## Plot histogram

# Plot
plt = plt_null_pvals_hist(pvals_tmml, pvals_tnme, pvals_tcmd)

# Store 
savefig("$(aim_dir)/Null-Pvals/Histogram-Null-Pvals.pdf")

## Plot ECDF

# Plot
plt = plt_null_pvals_ecdf(pvals_tmml, pvals_tnme, pvals_tcmd)

# Store 
savefig("$(aim_dir)/Null-Pvals/ECDF-Null-Pvals.pdf")
