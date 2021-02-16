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
        ϕ1_max = [+1.0,+10.0,2.0]
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
    next_eff_size = minimum(EF_PERC_RNG)
    eff_size_check = fill(false, length(EF_PERC_RNG))
    while ! all(eff_size_check)
        
        # Scale vector 
        if eff_type == "mml"
            ϕ1 = [eff_perc * ϕ1_max[1],eff_perc * ϕ1_max[2],ϕ1_max[3]]
            ϕ2 = [eff_perc * ϕ2_max[1],eff_perc * ϕ2_max[2],ϕ2_max[3]]
        elseif eff_type == "nme"
            ϕ1 = eff_perc * ϕ1_max
            ϕ2 = eff_perc * ϕ2_max
        else
            ϕ1 = ϕ1_max * ( 1.0 .+ eff_perc)
            ϕ2 = ϕ2_max / ( 1.0 .+ eff_perc)
        end
        
        # Baseline struct
        rs_g1 = create_baseline_struct(config)
        rs_g2 = create_baseline_struct(config)

        # Set vector
        rs_g1.ϕhat = ϕ1
        rs_g2.ϕhat = ϕ2
        
        # Compute stats
        proc_rep_rs!(rs_g1, config)
        proc_rep_rs!(rs_g2, config)

        # Non-empty analysis regions
        kp_in = [int != 0:0 for int in rs_g1.nls_rgs.cpg_ind]
        
        # Effect size
        if eff_type == "mml"
            eff_size = mean(abs.(rs_g1.mml[kp_in] - rs_g2.mml[kp_in]))
        elseif eff_type == "nme"
            eff_size = mean(abs.(rs_g1.nme[kp_in] - rs_g2.nme[kp_in]))
        elseif eff_type == "cmd"
            eff_size = mean(CpelNano.comp_cmd(rs_g1, rs_g2)[kp_in])
        end
        
        println("eff_size=$(eff_size)")

        # Record vector
        if eff_size >= next_eff_size
            
            # Record vector
            push!(ϕ1_arr, ϕ1)
            push!(ϕ2_arr, ϕ2)
            push!(eff_size_arr, eff_size)

            # Break
            i == length(EF_PERC_RNG) && break
            
            # Next effect size
            next_eff_size = EF_PERC_RNG[i + 1]
    
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

    # Initialize
    rs = CpelNano.RegStruct()

    # Read in file
    rs_data = readdlm("$(aim_dir)/estimation_region_chr_chr1_42900254_42903250.txt")

    # Get values
    rs.chrst = 42900254
    rs.chrend = 42903250
    rs.cpg_pos = Int64.(rs_data[:,1])
    rs.ρn = Float64.(rs_data[:,2])
    rs.dn = Float64.(rs_data[:,3])

    # Get rest of values
    rs.N = length(rs.cpg_pos)
    rs.cpg_grps = CpelNano.get_grps_from_cgs(rs.cpg_pos, config.min_grp_dist)
    rs.Nl = map(grp -> length(grp.cpg_ind), rs.cpg_grps)
    rs.ρl = CpelNano.get_ρl(rs.ρn, rs.cpg_grps)
    rs.dl = CpelNano.get_dl(rs.cpg_pos, rs.cpg_grps)
    rs.L = length(rs.cpg_grps)
    
    # Get analysis region info 
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
        rs_aux = CpelNano.cpel_samp_ont(M, αs1, βs1, POBS, LEVEL_MEAN_M, LEVEL_STD_M, LEVEL_MEAN_U, LEVEL_STD_U)
        
        # Fill fields
        cp_rs_flds!(rs, rs_aux)

        # Push struct
        push!(ms_g1, rs_aux)

        ## Group 2
        
        # Generate data
        rs_aux = CpelNano.cpel_samp_ont(M, αs2, βs2, POBS, LEVEL_MEAN_M, LEVEL_STD_M, LEVEL_MEAN_U, LEVEL_STD_U)
        
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
@everywhere function pmap_run_roc_sim(n_grp_reps::Int64, ϕ1::Vector{Float64}, ϕ2::Vector{Float64}, config::CpelNano.CpelNanoConfig)

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
@everywhere function pmap_run_null_sim(n_grp_reps::Int64, ϕ::Vector{Float64}, config::CpelNano.CpelNanoConfig)

    CpelNano.print_log("Producing set of p-values...")

    # Generate data
    ms_g1, ms_g2 = gen_grp_data(n_grp_reps, ϕ, ϕ, config)
        
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
    
    # Return test out
    return test_out

end
function run_roc_sim(n_sim_reps::Int64, n_grp_reps::Int64, eff_type::String, matched::Bool)

    # Init
    pval_thrs_rng = 0.0:0.001:1.0
    tpr = Array{Float64}(undef, length(EF_PERC_RNG), length(pval_thrs_rng))
    fpr = Array{Float64}(undef, length(EF_PERC_RNG), length(pval_thrs_rng))

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
        pmap_out = pmap(i -> pmap_run_roc_sim(n_grp_reps, ϕ1, ϕ2, config), 1:n_sim_reps)

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
function run_null_sim(n_sim_reps::Int64, n_grp_reps::Int64, matched::Bool)
    
    # Init
    ϕ = fill(0.0, 3)

    # CpelNano config
    config = CpelNano.CpelNanoConfig()
    config.max_size_subreg = 350
    config.matched = matched
    config.max_em_iters = 25
    config.max_em_init = 5
    config.verbose = false
    
    # Generate data for multiple reps
    pmap_out = pmap(i -> pmap_run_null_sim(n_grp_reps, ϕ, config), 1:n_sim_reps)

    # Get pvals
    pvals_tmml = []
    pvals_tnme = []
    pvals_tcmd = []
    for x in pmap_out
        isnothing(x) && continue
        push!(pvals_tmml, get_pvals(x, "mml"))
        push!(pvals_tnme, get_pvals(x, "nme"))
        push!(pvals_tcmd, get_pvals(x, "cmd"))
    end
    
    # Flatten vectors
    pvals_tmml = vcat(pvals_tmml...)
    pvals_tnme = vcat(pvals_tnme...)
    pvals_tcmd = vcat(pvals_tcmd...)

    # Remove NaN
    is_nan = isnan.(pvals_tmml)
    pvals_tmml = pvals_tmml[.!is_nan]
    is_nan = isnan.(pvals_tnme)
    pvals_tnme = pvals_tnme[.!is_nan]
    is_nan = isnan.(pvals_tcmd)
    pvals_tcmd = pvals_tcmd[.!is_nan]
   
    # Return pvals
    return pvals_tmml, pvals_tnme, pvals_tcmd

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
function plt_null_pvals_hist(pvals_tmml, pvals_tnme, pvals_tcmd, matched)
    
    # Vars
    y_lim = 0.5
    n_bins = 20
    y_label = "percentage (%)"

    # Tmml
    h = fit(Histogram, pvals_tmml, 0:1 / n_bins:1.0)
    midpts = collect(h.edges[1]) .- 2 / n_bins / n_bins
    percs = h.weights ./ sum(h.weights)
    plt_tmml = plot(midpts, percs,  seriestype=:bar, label="", xlabel="", 
        ylabel=y_label, title="Tmml", xlim=(0, 1), ylim=y_lim, size=(600, 500))

    # Tnme
    h = fit(Histogram, pvals_tnme, 0:1 / n_bins:1.0)
    midpts = collect(h.edges[1]) .- 2 / n_bins
    percs = h.weights ./ sum(h.weights)
    if matched
        plt_tnme = plot(midpts, percs,  seriestype=:bar, label="", xlabel="p-value", 
            ylabel=y_label, title="Tnme", xlim=(0, 1), ylim=y_lim, size=(600, 500))
    else
        plt_tnme = plot(midpts, percs,  seriestype=:bar, label="", xlabel="", 
            ylabel=y_label, title="Tnme", xlim=(0, 1), ylim=y_lim, size=(600, 500))
    end

    # Tcmd
    h = fit(Histogram, pvals_tcmd, 0:1 / n_bins:1.0)
    midpts = collect(h.edges[1]) .- 2 / n_bins
    percs = h.weights ./ sum(h.weights)
    plt_tcmd = plot(midpts, percs,  seriestype=:bar, label="", xlabel="p-value", 
        ylabel=y_label, title="Tcmd", xlim=(0, 1), ylim=y_lim, size=(600, 500))
    
    # Collage
    if matched
        plt = plot(plt_tmml, plt_tnme, layout=(2, 1), size=(600, 800))
    else
        plt = plot(plt_tmml, plt_tnme, plt_tcmd, layout=(3, 1), size=(600, 1200))
    end

    # Return plot
    return plt

end
function plt_null_pvals_ecdf(pvals_tmml, pvals_tnme, pvals_tcmd, matched)
    
    # Vars
    y_lim = (0, 1.0)
    x_lim = (0, 1.0)
    y_label = "ecdf"
    pval_rng = 0.01:0.01:1.0

    # Tmml
    F = ecdf(pvals_tmml)
    plt_tmml = plot(xlabel=" ", ylabel=y_label, title="Tmml", xlim=x_lim, ylim=y_lim, size=(600, 500))
    plot!(plt_tmml, pval_rng, pval_rng, seriestype=:line, label="Ideal")
    plot!(pval_rng, F.(pval_rng), seriestype=:line, label="Empirical")

    # Tnme
    F = ecdf(pvals_tnme)
    if matched
        plt_tnme = plot(xlabel="p-value", ylabel=y_label, title="Tnme", xlim=x_lim, ylim=y_lim, size=(600, 500))
        plot!(plt_tnme, pval_rng, pval_rng, seriestype=:line, label="Ideal")
        plot!(pval_rng, F.(pval_rng), seriestype=:line, label="Empirical")
    else
        plt_tnme = plot(xlabel=" ", ylabel=y_label, title="Tnme", xlim=x_lim, ylim=y_lim, size=(600, 500))
        plot!(plt_tnme, pval_rng, pval_rng, seriestype=:line, label="Ideal")
        plot!(pval_rng, F.(pval_rng), seriestype=:line, label="Empirical")
    end

    # Tcmd
    F = ecdf(pvals_tcmd)
    plt_tcmd = plot(xlabel="p-value", ylabel=y_label, title="Tcmd", xlim=x_lim, ylim=y_lim, size=(600, 500))
    plot!(plt_tcmd, pval_rng, pval_rng, seriestype=:line, label="Ideal")
    plot!(pval_rng, F.(pval_rng), seriestype=:line, label="Empirical")
    
    # Collage
    if matched
        plt = plot(plt_tmml, plt_tnme, layout=(2, 1), size=(600, 800))
    else
        plt = plot(plt_tmml, plt_tnme, plt_tcmd, layout=(3, 1), size=(600, 1200))
    end

    # Return plot
    return plt

end
