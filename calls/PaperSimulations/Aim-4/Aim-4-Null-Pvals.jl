#################################################################################################
# AIM 4: uniform distribution of null p-values
#################################################################################################
## Deps
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("./") # <=
@everywhere using CpelNano
using StatsBase
using StatsPlots
using Distributions

## Source code
root_dir = dirname(dirname(pathof(CpelNano)))
include("$(root_dir)/calls/PaperSimulations/Aim-4/Aim-4-Source.jl")

## Constants

# IO
const aim_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-4/" # <=
const blind_friend_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]

# Set up
@everywhere const M = 5
@everywhere const POBS = 1.0
@everywhere const LEVEL_STD_M = 2.0
@everywhere const LEVEL_STD_U = 2.0
@everywhere const LEVEL_MEAN_M = 84.0
@everywhere const LEVEL_MEAN_U = 80.0
@everywhere const EF_PERC_RNG = 0.1:0.1:0.1

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#################################################################################################
# Calls
#################################################################################################

# Args
matched = parse(Bool, ARGS[1])
n_sim_reps = parse(Int64, ARGS[2])
n_grp_reps = parse(Int64, ARGS[3])

# Generate null p-values
pvals_tmml, pvals_tnme, pvals_tcmd = run_null_sim(n_sim_reps, n_grp_reps, matched)

## Plot histogram

# Plot
plt = plt_null_pvals_hist(pvals_tmml, pvals_tnme, pvals_tcmd, matched)

# Store 
if matched
    savefig("$(aim_dir)/Null-Pvals/Histogram-Null-Pvals-Matched.pdf")
else
    savefig("$(aim_dir)/Null-Pvals/Histogram-Null-Pvals-Unmatched.pdf")
end

## Plot ECDF

# Plot
plt = plt_null_pvals_ecdf(pvals_tmml, pvals_tnme, pvals_tcmd, matched)

# Store 
if matched
    savefig("$(aim_dir)/Null-Pvals/ECDF-Null-Pvals-Matched.pdf")
else
    savefig("$(aim_dir)/Null-Pvals/ECDF-Null-Pvals-Unmatched.pdf")
end
