#################################################################################################
# AIM 5: two-sample comparison performance
#################################################################################################
## Deps
using StatsPlots
using Distributed
using Combinatorics
using Distributions
using DelimitedFiles
@everywhere using Pkg
@everywhere Pkg.activate("./") # <=
@everywhere using CpelNano

## Source code
root_dir = dirname(dirname(pathof(CpelNano)))
include("$(root_dir)/calls/PaperSimulations/Aim-5/Aim-5-Source.jl")

## Constants

# IO
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-5/" # <=
const blind_friend_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]

# Set up
@everywhere const M = 5
@everywhere const POBS = 1.0
@everywhere const Δ_LEVEL_MEAN = 5.0
@everywhere const LEVEL_STD_M = 3.0 / sqrt(4000 / 450)
@everywhere const LEVEL_STD_U = 3.0 / sqrt(4000 / 450)
@everywhere const LEVEL_MEAN_M = Δ_LEVEL_MEAN / 2.0
@everywhere const LEVEL_MEAN_U = -Δ_LEVEL_MEAN / 2.0
@everywhere const EF_PERC_RNG = 0.05:0.05:0.25

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#################################################################################################
# Calls
#################################################################################################

# Arguments
eff_type = ARGS[1]
n_sim_reps = parse(Int64, ARGS[2])
CpelNano.print_log("Running ROC simulations for two-sample test $(eff_type) w/ $(n_sim_reps) reps")

# Run
eff_size_arr, tpr, fpr = run_sim(n_sim_reps, eff_type)

# Plot
plt = plt_roc(eff_size_arr, tpr, fpr, "ROC - $(eff_type) - two-sample test ")
savefig(plt, "$(data_dir)/ROC-Two-Sample-Test-$(eff_type).pdf")
