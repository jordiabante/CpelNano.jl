"""
    `creat_cp_rs_struct(RS_REF)`

    Function that returns a RegStruct object RS with fields taken from RS_REF.

    # Examples
    ```julia-repl
    julia> CpelNano.creat_cp_rs_struct(rs_ref)
    ```
"""
function cp_rs_flds(rs_ref::RegStruct)::RegStruct

    # Create structure
    rs = deepcopy(rs_ref)

    # Reset fields
    rs.m = 0
    rs.nme = []
    rs.mml = []
    rs.ϕhat = []
    rs.calls = []
    rs.logZ = NaN
    rs.tm = TransferMat()
    rs.exps = Expectations()
    rs.nls_rgs = AnalysisRegions()

    # Return
    return rs

end
"""
    `pmap_two_samp_null_stats(MOD_REF,CALLS_S1,CALLS_S2,PERM,CONFIG)`

    Function that computes set of null statistics for a given data permutation.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_two_samp_null_stats(mod_ref,calls_s1,calls_s2,perm,config)
    ```
"""
function pmap_two_samp_null_stats(rs_ref::RegStruct, calls_s1::Vector{Vector{MethCallCpgGrp}}, calls_s2::Vector{Vector{MethCallCpgGrp}}, perm_ids::Vector{Int64}, config::CpelNanoConfig)::NTuple{3,Vector{Float64}}

    # Create two aux structs
    aux_rs_s1 = cp_rs_flds(rs_ref)
    aux_rs_s2 = cp_rs_flds(rs_ref)
    aux_rs_s1.m = length(calls_s1)
    aux_rs_s2.m = length(calls_s2)

    # Assign calls
    aux_calls = vcat(calls_s1, calls_s2)
    aux_rs_s1.calls = aux_calls[perm_ids]
    deleteat!(aux_calls, perm_ids)
    aux_rs_s2.calls = aux_calls

    # Init output
    tmml = fill(NaN, rs_ref.nls_rgs.num)
    tnme = fill(NaN, rs_ref.nls_rgs.num)
    tcmd = fill(NaN, rs_ref.nls_rgs.num)

    # Do estimation
    get_ϕhat!(aux_rs_s1, config)
    get_ϕhat!(aux_rs_s2, config)

    # Leave if we couldn't estimate ϕ
    (length(aux_rs_s1.ϕhat) > 0 && length(aux_rs_s2.ϕhat) > 0) || return tmml, tnme, tcmd
    print_log("ϕhat1=$(aux_rs_s1.ϕhat); ϕhat2=$(aux_rs_s2.ϕhat)")

    # Re-scale to per cpg site resolution
    rscle_grp_mod!(aux_rs_s1)
    rscle_grp_mod!(aux_rs_s2)
    
    # Statistical summaries
    get_stat_sums!(aux_rs_s1, config)
    get_stat_sums!(aux_rs_s2, config)

    # Compute stats
    tmml = abs.(aux_rs_s1.mml - aux_rs_s2.mml)
    tnme = abs.(aux_rs_s1.nme - aux_rs_s2.nme)
    tcmd = comp_cmd(aux_rs_s1, aux_rs_s2)
    print_log("tmml=$(tmml[1]); tnme=$(tnme[1]); tcmd=$(tcmd[1])")

    # Return tuple
    return tmml, tnme, tcmd

end
"""
    `pmap_diff_two_samp_comp(MODS_CHR_S1,MODS_CHR_S2,CONFIG)`

    Function that performs differential methylation analysis between two samples given two  
    CpelNano models.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_diff_two_samp_comp(mods_chr_s1,mods_chr_s2,config)
    ```
"""
function pmap_diff_two_samp_comp(mod_s1::RegStruct, mod_s2::RegStruct, nano_s1::String, nano_s2::String, fa_rec::FASTA.Record, config::CpelNanoConfig)::RegStatTestStruct

    print_log("Working on: $(mod_s1.chr)-$(mod_s1.chrst)-$(mod_s1.chrend)")
    
    # Init output struct
    test_struct = RegStatTestStruct(mod_s1)

    # Get group calls
    get_calls!(mod_s1, nano_s1, config)
    get_calls!(mod_s2, nano_s2, config)
    
    # Observed statistics
    tmml_obs = abs.(mod_s1.mml - mod_s2.mml)
    tnme_obs = abs.(mod_s1.nme - mod_s2.nme)
    tcmd_obs = comp_cmd(mod_s1, mod_s2)

    # Init p-values
    tmml_pvals = fill(NaN, length(tmml_obs))
    tnme_pvals = fill(NaN, length(tnme_obs))
    tcmd_pvals = fill(NaN, length(tcmd_obs))

    # Compute number of possible randomizations
    L = binomial(length(mod_s1.calls) + length(mod_s2.calls), length(mod_s1.calls))

    # If enough permutations
    if L > 20
        
        #  Check if exact test
        exact = L < LMAX_TWO_SAMP

        ## Retrieve region calls
        
        # Go back to group scale
        get_grp_info!(mod_s1, fa_rec, config.min_grp_dist)
        get_grp_info!(mod_s2, fa_rec, config.min_grp_dist)
    
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
            ind_subset = rand(1:L, LMAX_TWO_SAMP)
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
        # config.verbose && print_log("tmml_perms=$(tmml_perms)")
        # config.verbose && print_log("tnme_perms=$(tnme_perms)")
        # config.verbose && print_log("tcmd_perms=$(tcmd_perms)")

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
    # config.verbose && print_log("$(test_struct.tests)")
    # config.verbose && readline()

    # Return test struct
    return test_struct

end
"""
    `diff_two_samp_comp(MODEL_FILE_G1,MODEL_FILE_G2,FASTA,CONFIG)`

    Function that performs differential methylation analysis between two samples given two 
    CpelNano model files.

    # Examples
    ```julia-repl
    julia> CpelNano.diff_two_samp_comp(mods_path_s1,mods_path_s2,fasta,config)
    ```
"""
function diff_two_samp_comp(mods_path_s1::String, mods_path_s2::String, nano_s1::String, nano_s2::String, fasta::String, config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names, chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

    # Get output file names
    config.out_diff_files = OutputDiffFiles(config.out_dir, config.out_prefix)

    # Check if output file exists
    check_diff_output_exists(config) && return nothing

    # Loop over chromosomes
    for (i, chr) in enumerate(chr_names)
        
        print_log("Processing chromosome $(chr)")

        ## Get fasta record

        fa_rec = get_chr_fa_rec(chr, fasta)

        ## Read in models

        mods_chr_s1 = read_mod_file_chr(mods_path_s1, chr, fa_rec, config)
        mods_chr_s2 = read_mod_file_chr(mods_path_s2, chr, fa_rec, config)

        ## Find regions with at least one sample per group
        
        comm_regs = find_comm_regs([mods_chr_s1], [mods_chr_s2])
        
        # If no data in chromosome continue
        length(comm_regs) > 0 || continue
        config.verbose && print_log("Comm regions: $(comm_regs)")

        ## Process each region
        pmap_out = map(x -> pmap_diff_two_samp_comp(mods_chr_s1[x], mods_chr_s2[x], nano_s1, nano_s2, fa_rec, config), comm_regs)

        ## Write chromosome output
        write_diff_output(pmap_out, config)

    end

    # Return nothing
    return nothing

end