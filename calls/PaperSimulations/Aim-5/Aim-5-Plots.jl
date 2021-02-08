#################################################################################################
# AIM 5: two-sample comparison performance
#################################################################################################
## Deps
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("./") # <=
@everywhere using CpelNano
using StatsPlots
using Distributions
using Combinatorics

## Constants

# IO
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-5/" # <=
const blind_friend_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]

# Set up
@everywhere const M = 5
@everywhere const POBS = 1.0
@everywhere const LEVEL_STD_M = 2.0
@everywhere const LEVEL_STD_u = 2.0
@everywhere const LEVEL_MEAN_M = 84.0
@everywhere const LEVEL_MEAN_u = 80.0
@everywhere const LMAX_TWO_SAMP_AIM_5 = 100
@everywhere const EFF_PERC_RNG = 0.05:0.05:0.25

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
    next_eff_size = minimum(EFF_PERC_RNG)
    eff_size_check = fill(false, length(EFF_PERC_RNG))
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
            i == length(EFF_PERC_RNG) && break
            
            # Next effect size
            next_eff_size = EFF_PERC_RNG[i + 1]
    
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
@everywhere function gen_grp_data(ϕ1::Vector{Float64}, ϕ2::Vector{Float64}, config::CpelNano.CpelNanoConfig)

    ## Set baseline struct
    rs = create_baseline_struct(config)
    αs1, βs1 = CpelNano.get_αβ_from_ϕ(ϕ1, rs)
    αs2, βs2 = CpelNano.get_αβ_from_ϕ(ϕ2, rs)
    
    ## Group 1
    
    # Generate calls
    ms_g1 = CpelNano.cpel_samp_ont(M, αs1, βs1, POBS, LEVEL_MEAN_M, LEVEL_STD_M, LEVEL_MEAN_u, LEVEL_STD_u)
    
    # Fill fields
    cp_rs_flds!(rs, ms_g1)

    ## Group 2
    
    # Generate data
    ms_g2 = CpelNano.cpel_samp_ont(M, αs2, βs2, POBS, LEVEL_MEAN_M, LEVEL_STD_M, LEVEL_MEAN_u, LEVEL_STD_u)
    
    # Fill fields
    cp_rs_flds!(rs, ms_g2)
    
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
@everywhere function pmap_two_samp_null_stats(rs_ref::CpelNano.RegStruct, calls_s1::Vector{Vector{CpelNano.MethCallCpgGrp}}, calls_s2::Vector{Vector{CpelNano.MethCallCpgGrp}}, perm_ids::Vector{Int64}, config::CpelNano.CpelNanoConfig)::NTuple{3,Vector{Float64}}

    # Init output
    tmml = fill(NaN, rs_ref.nls_rgs.num)
    tnme = fill(NaN, rs_ref.nls_rgs.num)
    tcmd = fill(NaN, rs_ref.nls_rgs.num)

    # Create two aux structs
    aux_rs_s1 = CpelNano.RegStruct()
    aux_rs_s2 = CpelNano.RegStruct()
    cp_rs_flds!(rs_ref, aux_rs_s1)
    cp_rs_flds!(rs_ref, aux_rs_s2)
    aux_rs_s1.m = length(calls_s1)
    aux_rs_s2.m = length(calls_s2)

    # Assign calls
    aux_calls = vcat(calls_s1, calls_s2)
    aux_rs_s1.calls = aux_calls[perm_ids]
    deleteat!(aux_calls, perm_ids)
    aux_rs_s2.calls = aux_calls

    # Do estimation
    proc_rep_rs!(aux_rs_s1, config)
    proc_rep_rs!(aux_rs_s2, config)
    (aux_rs_s1.proc && aux_rs_s2.proc) || return tmml, tnme, tcmd
        
    # Compute stats
    tmml = abs.(aux_rs_s1.mml - aux_rs_s2.mml)
    tnme = abs.(aux_rs_s1.nme - aux_rs_s2.nme)
    tcmd = CpelNano.comp_cmd(aux_rs_s1, aux_rs_s2)

    # Return tuple
    return tmml, tnme, tcmd

end
function pmap_diff_two_samp_comp(mod_s1::CpelNano.RegStruct, mod_s2::CpelNano.RegStruct, config::CpelNano.CpelNanoConfig)::CpelNano.RegStatTestStruct

    CpelNano.print_log("Processing estimation region")
    
    # Init output struct
    test_struct = CpelNano.RegStatTestStruct(mod_s1)
    
    # Observed statistics
    tmml_obs = abs.(mod_s1.mml - mod_s2.mml)
    tnme_obs = abs.(mod_s1.nme - mod_s2.nme)
    tcmd_obs = CpelNano.comp_cmd(mod_s1, mod_s2)

    # Init p-values
    tmml_pvals = fill(NaN, length(tmml_obs))
    tnme_pvals = fill(NaN, length(tnme_obs))
    tcmd_pvals = fill(NaN, length(tcmd_obs))

    # Compute number of possible randomizations
    L = binomial(length(mod_s1.calls) + length(mod_s2.calls), length(mod_s1.calls))

    # If enough permutations
    if L > 20
        
        #  Check if exact test
        exact = L < LMAX_TWO_SAMP_AIM_5
    
        # Store calls in aux vars
        calls_s1 = deepcopy(mod_s1.calls)
        calls_s2 = deepcopy(mod_s2.calls)

        # Create iteratable object with all combinations
        comb_iter = combinations(1:(length(calls_s1) + length(calls_s2)), length(calls_s1))

        # Get sample label combinations to use
        comb_iter_used = []
        if exact
            # Use all sample assignments
            comb_iter_used = comb_iter
        else
            # Use Lmax group assignments
            ind_subset = rand(1:L, LMAX_TWO_SAMP_AIM_5)
            @inbounds for (ind, comb) in enumerate(comb_iter)
                (ind in ind_subset) && push!(comb_iter_used, comb)
            end
        end

        ## Null statistics

        # Compute null statistics
        pmap_out = pmap(perm -> pmap_two_samp_null_stats(mod_s1, calls_s1, calls_s2, perm, config), comb_iter_used)

        # Distribute statistics
        tmml_perms = [x[1] for x in pmap_out]
        tnme_perms = [x[2] for x in pmap_out]
        tcmd_perms = [x[3] for x in pmap_out]
        CpelNano.print_log("tmml_perms=$(tmml_perms)")
        CpelNano.print_log("tnme_perms=$(tnme_perms)")
        CpelNano.print_log("tcmd_perms=$(tcmd_perms)")

        ## P-value computation

        # Loop over analysis regions
        @inbounds for k in 1:length(tmml_obs)

            # Check if data
            isnan(tmml_obs[k]) && continue

            # Get permutation stats from k-th analysis region
            tmml_perms_k = [perm[k] for perm in tmml_perms]
            tnme_perms_k = [perm[k] for perm in tnme_perms]
            tcmd_perms_k = [perm[k] for perm in tcmd_perms]

            # Clean NaNs in null stats
            tmml_perms_k = tmml_perms_k[.!isnan.(tmml_perms_k)]
            tnme_perms_k = tnme_perms_k[.!isnan.(tnme_perms_k)]
            tcmd_perms_k = tcmd_perms_k[.!isnan.(tcmd_perms_k)]

            # Check enough null stats after filtering NaNs
            # length(tmml_perms_k) > 20 || continue

            # Get number of permutation stats equal or above observed
            tmml_pval_k = sum(abs.(tmml_perms_k) .>= abs(tmml_obs[k]))
            tnme_pval_k = sum(abs.(tnme_perms_k) .>= abs(tnme_obs[k]))
            tcmd_pval_k = sum(tcmd_perms_k .>= tcmd_obs[k])

            # Compute p-values
            tmml_pvals[k] = exact ? tmml_pval_k / length(tmml_perms_k) : (1.0 + tmml_pval_k) / (1.0 + length(tmml_perms_k))
            tnme_pvals[k] = exact ? tnme_pval_k / length(tnme_perms_k) : (1.0 + tnme_pval_k) / (1.0 + length(tnme_perms_k))
            tcmd_pvals[k] = exact ? tcmd_pval_k / length(tcmd_perms_k) : (1.0 + tcmd_pval_k) / (1.0 + length(tcmd_perms_k))

        end

    end

    # Fill return object
    @inbounds for k in 1:length(tmml_obs)
        test_struct.tests.tmml_test[k] = (tmml_obs[k], tmml_pvals[k])
        test_struct.tests.tnme_test[k] = (tnme_obs[k], tnme_pvals[k])
        test_struct.tests.tcmd_test[k] = (tcmd_obs[k], tcmd_pvals[k])
    end
    CpelNano.print_log("test_struct.tests=$(test_struct.tests)")

    # Return test struct
    return test_struct

end
@everywhere function map_run_sim(ϕ1::Vector{Float64}, ϕ2::Vector{Float64}, config::CpelNano.CpelNanoConfig)

    # Choose whether effect or not (50% chance)
    eff_bool_j = rand() > 0.5
        
    # Generate data
    if eff_bool_j
        # Introduce effect
        m_g1, m_g2 = gen_grp_data(ϕ1, ϕ2, config)
    else
        # Do not introduce effect
        m_g1, m_g2 = gen_grp_data(ϕ1, ϕ1, config)
    end

    # Estimate parameters
    proc_rep_rs!(m_g1, config)
    proc_rep_rs!(m_g2, config)

    # Find valid reps
    (m_g1.proc && m_g2.proc) || return nothing

    # Perform test
    test_out = pmap_diff_two_samp_comp(m_g1, m_g2, config)
    
    # Return test out and effect bool
    return eff_bool_j, test_out

end
function run_sim(n_sim_reps::Int64, eff_type::String)

    # Init
    pval_thrs_rng = 0.0:0.001:1.0
    tpr = Array{Float64}(undef, length(EFF_PERC_RNG), length(pval_thrs_rng))
    fpr = Array{Float64}(undef, length(EFF_PERC_RNG), length(pval_thrs_rng))

    # CpelNano config
    config = CpelNano.CpelNanoConfig()
    config.max_size_subreg = 350
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
        map_out = map(i -> map_run_sim(ϕ1, ϕ2, config), 1:n_sim_reps)

        # Get pvals
        pvals_neg = []
        pvals_pos = []
        for x in map_out
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
        CpelNano.print_log("pvals_pos=$(pvals_pos)")
        CpelNano.print_log("pvals_neg=$(pvals_neg)")

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
n_sim_reps = parse(Int64, ARGS[2])
CpelNano.print_log("Running ROC simulations for two-sample test $(eff_type) w/ $(n_sim_reps) reps")

# Run
eff_size_arr, tpr, fpr = run_sim(n_sim_reps, eff_type)

# Plot
plt = plt_roc(eff_size_arr, tpr, fpr, "ROC - $(eff_type) - two-sample test ")
savefig(plt, "$(data_dir)/ROC-Two-Sample-Test-$(eff_type).pdf")
