#################################################################################################
# AIM 4: group comparison performance
#################################################################################################
## Deps
using StatsPlots
using Distributed
using Distributions
using DelimitedFiles
@everywhere using Pkg
@everywhere Pkg.activate("./") # <=
@everywhere using CpelNano

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
@everywhere const Δ_LEVEL_MEAN = 0.13
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
matched = parse(Bool, ARGS[2])
n_sim_reps = parse(Int64, ARGS[3])
CpelNano.print_log("Running ROC simulations for $(matched) $(eff_type) w/ $(n_sim_reps) reps")

# Run
if matched

    ### Matched
    
    # Vars
    n_grp_reps = 6

    # Simulations
    eff_size_arr, tpr, fpr = run_roc_sim(n_sim_reps, n_grp_reps, eff_type, matched)

    # Plot
    plt = plt_roc(eff_size_arr, tpr, fpr, "ROC - $(eff_type) - matched ")
    savefig(plt, "$(aim_dir)/ROC/ROC-Matched-$(eff_type).pdf")

else

    ### Unmatched

    # Vars
    n_grp_reps = 5

    # Simulations
    eff_size_arr, tpr, fpr = run_roc_sim(n_sim_reps, n_grp_reps, eff_type, matched)

    # Plot
    plt = plt_roc(eff_size_arr, tpr, fpr, "ROC - $(eff_type) - unmatched ")
    savefig(plt, "$(aim_dir)/ROC/ROC-Unmatched-$(eff_type).pdf")

end
