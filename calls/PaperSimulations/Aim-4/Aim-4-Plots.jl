#################################################################################################
# AIM 4: group comparison performance
#################################################################################################
## Deps
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("./") <=
@everywhere using CpelNano
using StatsPlots
using Distributions
using DelimitedFiles

## Constants

# IO
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-4/" <=
const blind_friend_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]

# Set up
@everywhere const M = 5
@everywhere const pobs = 1.0
@everywhere const level_std_m = 2.0
@everywhere const level_std_u = 2.0
@everywhere const level_mean_m = 84.0
@everywhere const level_mean_u = 80.0
@everywhere const eff_perc_rng = 0.05:0.05:0.25

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#################################################################################################
# Functions
#################################################################################################
function get_ϕs_eff_rng(eff_type, config)
    
    # Pick max values
    if eff_type == "mml"
        ϕ1_max = [-1.0,-10.0,2.0]
        ϕ2_max = [+1.0,+10.0,2.0]
    elseif eff_type == "nme"
        ϕ1_max = [0.0,0.0,0.0]
        ϕ2_max = [0.0,0.0,5.0]
    elseif eff_type == "cmd" 
        ϕ1_max = [-1.0,-10.0,2.0]
        ϕ2_max = [+1.0,+10.0,2.0]
    else
        CpelNano.print_log("Effect type $(eff_type) type not recognized.)")
    end
    
    # Loop over percentages
    i = 1
    eff_perc = 0.0
    eff_size_arr = Vector{Float64}()
    ϕ1_arr = Vector{Vector{Float64}}()
    ϕ2_arr = Vector{Vector{Float64}}()
    next_eff_size = minimum(eff_perc_rng)
    eff_size_check = fill(false, length(eff_perc_rng))
    while ! all(eff_size_check)
        
        # Scale vector 
        ϕ1 = eff_perc * ϕ1_max
        ϕ2 = eff_perc * ϕ2_max
        
        # Baseline struct
        rs_g1 = create_baseline_struct(config)
        rs_g2 = create_baseline_struct(config)

        # Set vector
        rs_g1.ϕhat = ϕ1
        rs_g2.ϕhat = ϕ2
        
        # Compute stats
        proc_rep_rs!(rs_g1, config)
        proc_rep_rs!(rs_g2, config)

        # Effect size
        if eff_type == "mml"
            eff_size = mean(abs.(rs_g1.mml - rs_g2.mml))
        elseif eff_type == "nme"
            eff_size = mean(abs.(rs_g1.nme - rs_g2.nme))
        elseif eff_type == "cmd"
            eff_size = mean(CpelNano.comp_cmd(rs_g1, rs_g2))
        end
        
        # println("eff_size=$(eff_size)")

        # Record vector
        if eff_size >= next_eff_size
            
            # Record vector
            push!(ϕ1_arr, ϕ1)
            push!(ϕ2_arr, ϕ2)
            push!(eff_size_arr, eff_size)

            # Break
            i == length(eff_perc_rng) && break
            
            # Next effect size
            next_eff_size = eff_perc_rng[i + 1]
    
            # Increase counter
            i += 1

        end

        # Increase effect percentage
        eff_perc += 0.0005
        
    end

    # Return vector
    return eff_size_arr, ϕ1_arr, ϕ2_arr

end
@everywhere function create_baseline_struct(config::CpelNano.CpelNanoConfig)

    # Define
    rs = CpelNano.RegStruct()
    rs.chrst = 1
    rs.chrend = 3000
    rs.N = 30
    rs.L = 30
    rs.Nl = fill(1.0, rs.L)
    rs.ρl = fill(0.1, rs.L)
    rs.ρn = fill(0.1, rs.N)
    rs.dl = fill(100.0, rs.L - 1)
    rs.dn = fill(100.0, rs.N - 1)
    rs.cpg_pos = [rs.chrst + 100 * n - 50 for n = 1:rs.N]
    CpelNano.get_nls_reg_info!(rs, config)

    # Return struct
    return rs
    
end
@everywhere function cp_rs_flds!(rs::CpelNano.RegStruct, rs_aux::CpelNano.RegStruct)
    
    # Fill fields
    rs_aux.L = rs.L
    rs_aux.Nl = rs.Nl
    rs_aux.ρl = rs.ρl
    rs_aux.dl = rs.dl
    rs_aux.N = rs.N
    rs_aux.ρn = rs.ρn
    rs_aux.dn = rs.dn
    rs_aux.cpg_pos = rs.cpg_pos
    rs_aux.nls_rgs = rs.nls_rgs

    # Return nothing
    return nothing
    
end
@everywhere function gen_grp_data(n_grp_reps::Int64, ϕ1::Vector{Float64}, ϕ2::Vector{Float64}, config::CpelNano.CpelNanoConfig)

    ## Set baseline struct
    rs = create_baseline_struct(config)
    αs1, βs1 = CpelNano.get_αβ_from_ϕ(ϕ1, rs)
    αs2, βs2 = CpelNano.get_αβ_from_ϕ(ϕ2, rs)
    
    ## Generate group data
    ms_g1 = Vector{CpelNano.RegStruct}()
    ms_g2 = Vector{CpelNano.RegStruct}()
    for rep = 1:n_grp_reps

        ## Group 1
        
        # Generate calls
        rs_aux = CpelNano.cpel_samp_ont(M, αs1, βs1, pobs, level_mean_m, level_std_m, level_mean_u, level_std_u)
        
        # Fill fields
        cp_rs_flds!(rs, rs_aux)

        # Push struct
        push!(ms_g1, rs_aux)

        ## Group 2
        
        # Generate data
        rs_aux = CpelNano.cpel_samp_ont(M, αs2, βs2, pobs, level_mean_m, level_std_m, level_mean_u, level_std_u)
        
        # Fill fields
        cp_rs_flds!(rs, rs_aux)
        
        # Push struct
        push!(ms_g2, rs_aux)

    end

    # Returns group data
    return ms_g1, ms_g2

end
@everywhere function proc_rep_rs!(rs::CpelNano.RegStruct, config::CpelNano.CpelNanoConfig)

    # Estimate parameter vector ϕ
    if length(rs.ϕhat) == 0
        CpelNano.get_ϕhat!(rs, config)
        length(rs.ϕhat) > 0 || return nothing
    end
    
    # Re-scale to per cpg site resolution
    CpelNano.rscle_grp_mod!(rs)
    
    # Get matrices
    CpelNano.get_rs_lgtrck_mats!(rs)
    
    # Partition function
    CpelNano.get_rs_logZ!(rs)
    
    # Get them all
    CpelNano.get_rs_log_gs!(rs)
    
    # Compute E[log g(X)]
    CpelNano.get_rs_exp_log_g1!(rs)
    CpelNano.get_rs_exp_log_g2!(rs)
    
    # Compute E[X] and E[XX]
    CpelNano.get_rs_exps!(rs)
    
    # Compute μ(X)
    CpelNano.comp_mml!(rs)

    # Compute h(X)
    CpelNano.comp_nme!(rs)

    # Set processed
    rs.proc = true

    # Return nothing
    return nothing
    
end
function get_pvals(test_out::CpelNano.RegStatTestStruct, eff_type::String)

    # Get
    if eff_type == "mml"
        pvals = [x[2] for x in test_out.tests.tmml_test]
    elseif eff_type == "nme"
        pvals = [x[2] for x in test_out.tests.tnme_test]
    elseif eff_type == "cmd"
        pvals = [x[2] for x in test_out.tests.tcmd_test]
    else
        CpelNano.print_log("Wrong effect type $(eff_type)")
    end

    # Return p-values
    return pvals
    
end
@everywhere function pmap_run_sim(n_grp_reps::Int64, ϕ1::Vector{Float64}, ϕ2::Vector{Float64}, config::CpelNano.CpelNanoConfig)

    # Choose whether effect or not (50% chance)
    eff_bool_j = rand() > 0.5
        
    # Generate data
    if eff_bool_j
        # Introduce effect
        ms_g1, ms_g2 = gen_grp_data(n_grp_reps, ϕ1, ϕ2, config)
    else
        # Do not introduce effect
        ms_g1, ms_g2 = gen_grp_data(n_grp_reps, ϕ1, ϕ1, config)
    end

    # Estimate parameters
    for rep = 1:n_grp_reps
        proc_rep_rs!(ms_g1[rep], config)
        proc_rep_rs!(ms_g2[rep], config)
    end

    # Find valid reps
    n_val = 0
    for rep = 1:n_grp_reps
        rs_1 = ms_g1[rep]
        rs_2 = ms_g2[rep]
        n_val += (rs_1.proc && rs_2.proc) ? 1 : 0
    end
    n_val == n_grp_reps || return nothing

    # Perform test
    test_out = config.matched ? CpelNano.mat_est_reg_test(ms_g1, ms_g2) : CpelNano.unmat_est_reg_test(ms_g1, ms_g2)
    
    # Return test out and effect bool
    return eff_bool_j, test_out

end
function run_sim(n_sim_reps::Int64, n_grp_reps::Int64, eff_type::String, matched::Bool)

    # Init
    pval_thrs_rng = 0.0:0.001:1.0
    tpr = Array{Float64}(undef, length(eff_perc_rng), length(pval_thrs_rng))
    fpr = Array{Float64}(undef, length(eff_perc_rng), length(pval_thrs_rng))

    # CpelNano config
    config = CpelNano.CpelNanoConfig()
    config.max_size_subreg = 350
    config.matched = matched
    config.max_em_iters = 25
    config.max_em_init = 5
    config.verbose = false

    # Get vectors
    eff_size_arr, ϕ1_arr, ϕ2_arr = get_ϕs_eff_rng(eff_type, config)

    # Loop over percent of effect size
    for i = 1:length(ϕ1_arr)
        
        # Init
        ϕ1 = ϕ1_arr[i]
        ϕ2 = ϕ2_arr[i]
        CpelNano.print_log("Effect percentage $(eff_size_arr[i])")
        
        # Generate data for multiple reps
        pmap_out = pmap(i -> pmap_run_sim(n_grp_reps, ϕ1, ϕ2, config), 1:n_sim_reps)

        # Get pvals
        pvals_neg = []
        pvals_pos = []
        for x in pmap_out
            isnothing(x) && continue
            if x[1]
                push!(pvals_pos, get_pvals(x[2], eff_type))
            else
                push!(pvals_neg, get_pvals(x[2], eff_type))
            end
        end
        
        # Flatten vectors
        pvals_pos = vcat(pvals_pos...)
        pvals_neg = vcat(pvals_neg...)

        # Remove NaN
        is_nan = isnan.(pvals_pos)
        pvals_pos = pvals_pos[.!is_nan]
        is_nan = isnan.(pvals_neg)
        pvals_neg = pvals_neg[.!is_nan]
        
        # Store
        for (j, p_thrs) in enumerate(pval_thrs_rng) 
           
            # Count true positives and true negatives
            if length(pvals_pos) > 0
                tp = sum(pvals_pos .<= p_thrs)
                fn = sum(pvals_pos .> p_thrs)
            else
                tp = NaN
                fn = NaN
            end
            if length(pvals_neg) > 0
                tn = sum(pvals_neg .> p_thrs)
                fp = sum(pvals_neg .<= p_thrs)
            else
                tn = NaN
                fp = NaN
            end

            # Store TPR & FPR
            tpr[i,j] = tp / (tp + fn + 1)
            fpr[i,j] = fp / (fp + tn + 1)

        end
        
    end
   
    return eff_size_arr, tpr, fpr

end
function plt_roc(eff_size_arr, tpr, fpr, ttl)
    
    # Add each effect size
    plt = plot(xlabel="FPR", ylabel="TPR", title=ttl, xlim=(0, 1), ylim=(0, 1), size=(700, 700), legend=:bottomright)
    for (i, eff_size) in enumerate(eff_size_arr)
        plot!(plt, fpr[i,:], tpr[i,:],  seriestype=:line, label="$(round(eff_size;digits=2))", color=blind_friend_col[i])
    end
    
    # Return plot
    return plt

end
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
    eff_size_arr, tpr, fpr = run_sim(n_sim_reps, n_grp_reps, eff_type, matched)

    # Plot
    plt = plt_roc(eff_size_arr, tpr, fpr, "ROC - $(eff_type) - matched ")
    savefig(plt, "$(data_dir)/ROC-Matched-$(eff_type).pdf")


else

    ### Unmatched

    # Vars
    n_grp_reps = 5

    # Simulations
    eff_size_arr, tpr, fpr = run_sim(n_sim_reps, n_grp_reps, eff_type, matched)

    # Plot
    plt = plt_roc(eff_size_arr, tpr, fpr, "ROC - $(eff_type) - unmatched ")
    savefig(plt, "$(data_dir)/ROC-Unmatched-$(eff_type).pdf")

end
