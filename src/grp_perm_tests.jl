##################################################################################################
## Unmatched test
##################################################################################################
"""
    `comp_unmat_stat_cmd(mods_g1,mods_g2)`
    Function that computes a test statistic Tcmd for each analysis region in unmatched case 
    between groups.
    # Examples
    ```julia-repl
    julia> CpelNano.comp_unmat_stat_cmd(mods_g1,mods_g2)
    ```
"""
function comp_unmat_stat_cmd(mods_g1::Vector{RegStruct}, mods_g2::Vector{RegStruct})::Vector{Float64}
    
    # Compute statistic
    cmds = []
    @inbounds for s1 = 1:length(mods_g1)
        @inbounds for s2 = 1:length(mods_g2)
            push!(cmds, comp_cmd(mods_g1[s1], mods_g2[s2]))
        end
    end

    # Return vector of average CMD between groups
    return sum(cmds) / (length(mods_g1) * length(mods_g2))

end
"""
    `get_unmat_cmd_tbl(mods_g1,mods_g2)`

    Function that computes CMD between all possible pairs.

    # Examples
    ```julia-repl
    julia> CpelNano.get_unmat_cmd_tbl(mods_g1,mods_g2)
    ```
"""
function get_unmat_cmd_tbl(mods_g1::Vector{RegStruct}, mods_g2::Vector{RegStruct})::Dict{NTuple{2,Int64},Vector{Float64}}
    
    # Number of analysis regions
    rs = mods_g1[1]
    num_nls_reg = rs.nls_rgs.num

    # Number of replicates
    n_g1 = length(mods_g1)
    n_g2 = length(mods_g2)
    n_tot = n_g1 + n_g2

    # Init matrix 
    cmd_tbl = Dict{NTuple{2,Int64},Vector{Float64}}()

    # Compute statistic within group 1
    @inbounds for s1 = 1:n_g1
        @inbounds for s2 = (s1+1):n_g1
            cmd_tbl[(s1,s2)] = comp_cmd(mods_g1[s1], mods_g1[s2])
        end
    end

    # Compute statistic within group 2
    @inbounds for s1 = 1:n_g2
        @inbounds for s2 = (s1+1):n_g2
            cmd_tbl[(n_g1+s1,n_g1+s2)] = comp_cmd(mods_g2[s1], mods_g2[s2])
        end
    end

    # Compute statistic between groups
    @inbounds for s1 = 1:n_g1
        @inbounds for s2 = 1:n_g2
            cmd_tbl[(s1,n_g1+s2)] = comp_cmd(mods_g1[s1], mods_g2[s2])
        end
    end

    # Return matrix of CMD
    return cmd_tbl

end
"""
    `comp_unmat_perm_stat_cmd(N_G1,N_G2,NUM_NLS_RGS,CMD_TBL,PERM_IDS)`

    Function that produces a permutation statistic for Tcmd in unmatched case.

    # Examples
    ```julia-repl
    julia> CpelNano.comp_unmat_perm_stat_cmd(n_g1,n_g2,num_nls_rgs,cmd_tbl,ind_g1)
    ```
"""
function comp_unmat_perm_stat_cmd(n_g1::Int64,n_g2::Int64,num_nls_rgs::Int64,cmd_tbl::Dict{Tuple{Int64,Int64},Array{Float64,1}},ind_g1::Vector{Int64})::Vector{Float64}

    # Replicates
    ind_g2 = collect(1:(n_g1 + n_g2))
    ind_g2 = deleteat!(ind_g2,ind_g1)

    # Get statistic
    cmd = fill(0.0,num_nls_rgs)
    @inbounds for i in ind_g1
        @inbounds for j in ind_g2
            if haskey(cmd_tbl,(i,j))
                cmd .+= cmd_tbl[(i,j)]
            else
                cmd .+= cmd_tbl[(j,i)]
            end
        end
    end

    # Return stat for permutation
    return cmd / (n_g1 * n_g2)
    
end
"""
    `comp_unmat_stat_mml(μs_g1,μs_g2)`

    Function that computes MML difference between groups (unmatched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_unmat_stat_mml(μs_g1,μs_g2)
```
"""
function comp_unmat_stat_mml(μs_g1::Vector{Vector{Float64}}, μs_g2::Vector{Vector{Float64}})::Vector{Float64}

    # Compute means within groups
    mml1 = sum(μs_g1) / length(μs_g1)
    mml2 = sum(μs_g2) / length(μs_g2)

    # Return Tmml
    return mml1 - mml2

end
"""
    `comp_unmat_perm_stat_mml(μs_g1,μs_g2,PERM_IDS)`

    Function that produces a permutation statistic for Tmml (unmatched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_unmat_perm_stat_mml(μs_g1, μs_g2, perm_ids)
```
"""
function comp_unmat_perm_stat_mml(μs_g1::Vector{Vector{Float64}}, μs_g2::Vector{Vector{Float64}}, perm_ids::Vector{Int64})::Vector{Float64}

    # Get vectors for each group
    μs_g2p = vcat(μs_g1, μs_g2)
    μs_g1p = μs_g2p[perm_ids]
    deleteat!(μs_g2p, perm_ids)

    # Return
    return comp_unmat_stat_mml(μs_g1p, μs_g2p)

end
"""
    `comp_unmat_stat_nme(hs_g1,hs_g2)`

    Function that computes NME difference between groups (unmatched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_unmat_stat_nme(hs_g1,hs_g2)
```
"""
function comp_unmat_stat_nme(hs_g1::Vector{Vector{Float64}}, hs_g2::Vector{Vector{Float64}})::Vector{Float64}

    # Compute h within group
    nme1 = sum(hs_g1) / length(hs_g1)
    nme2 = sum(hs_g2) / length(hs_g2)

    # Return difference in nme
    return nme1 - nme2

end
"""
    `comp_unmat_perm_stat_nme(hs_g1,hs_g2,PERM_IDS)`

    Function that produces a permutation statistic for Tnme (unmatched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_unmat_perm_stat_nme(hs_g1,hs_g2,perm_ids)
```
"""
function comp_unmat_perm_stat_nme(hs_g1::Vector{Vector{Float64}}, hs_g2::Vector{Vector{Float64}}, perm_ids::Vector{Int64})::Vector{Float64}

    # Get vectors for each group
    hs_g2p = vcat(hs_g1, hs_g2)
    hs_g1p = hs_g2p[perm_ids]
    deleteat!(hs_g2p, perm_ids)

    # Return
    return comp_unmat_stat_nme(hs_g1p, hs_g2p)

end
"""
    `unmat_nls_reg_test(MODELS_G1,MODELS_G2)`

    Function that performs hypothesis testing in unmatched samples group comparison in an estimation region.

    # Examples
    ```julia-repl
    julia> CpelNano.unmat_nls_reg_test(ms_g1,ms_g2)
```
"""
function unmat_est_reg_test(ms_g1::Vector{RegStruct}, ms_g2::Vector{RegStruct})::RegStatTestStruct

    # Initialize output object
    test_struct = RegStatTestStruct(ms_g1[1])

    # Get stats for k-th analysis region
    μs_g1 = [x.mml for x in ms_g1]
    μs_g2 = [x.mml for x in ms_g2]
    hs_g1 = [x.nme for x in ms_g1]
    hs_g2 = [x.nme for x in ms_g2]

    # Compute observed stats
    tmml_obs = comp_unmat_stat_mml(μs_g1, μs_g2)
    tnme_obs = comp_unmat_stat_nme(hs_g1, hs_g2)
    tcmd_obs = comp_unmat_stat_cmd(ms_g1, ms_g2)

    # Init p-values
    tmml_pvals = fill(NaN, length(tmml_obs))
    tnme_pvals = fill(NaN, length(tnme_obs))
    tcmd_pvals = fill(NaN, length(tcmd_obs))
    
    # Compute number of possible randomizations
    n_g1 = length(μs_g1)
    n_g2 = length(μs_g2)
    L = binomial(n_g1 + n_g2, n_g1)

    # If enough data compute p-values (p<0.05)
    if L>20
        
        # Check if exact test
        exact = L < LMAX

        # Create iteratable object with all combinations
        comb_iter = combinations(1:(n_g1+n_g2), n_g1)

        # Get group label combinations to use
        comb_iter_used = []
        if exact
            # Use all group assignments
            comb_iter_used = comb_iter
        else
            # Use Lmax group assignments
            ind_subset = rand(1:L, LMAX)
            @inbounds for (ind, comb) in enumerate(comb_iter)
                (ind in ind_subset) && push!(comb_iter_used, comb)
            end
        end

        # Compute all CMD pairs
        cmd_tbl = get_unmat_cmd_tbl(ms_g1,ms_g2)

        # Use method for random permutation
        tmml_perms = map(x -> comp_unmat_perm_stat_mml(μs_g1, μs_g2, x), comb_iter_used)
        tnme_perms = map(x -> comp_unmat_perm_stat_nme(hs_g1, hs_g2, x), comb_iter_used)
        tcmd_perms = map(x -> comp_unmat_perm_stat_cmd(n_g1, n_g2, length(tmml_obs), cmd_tbl, x), comb_iter_used)

        # Compute p-values two-sided test
        @inbounds for k in 1:length(tmml_obs)

            # Check if data
            isnan(tmml_obs[k]) && continue
            
            # Get permutation stats from k-th analysis region
            tmml_perms_k = [perm[k] for perm in tmml_perms]
            tnme_perms_k = [perm[k] for perm in tnme_perms]
            tcmd_perms_k = [perm[k] for perm in tcmd_perms]

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
        test_struct.tests.tmml_test[k] = (tmml_obs[k],tmml_pvals[k])
        test_struct.tests.tnme_test[k] = (tnme_obs[k],tnme_pvals[k])
        test_struct.tests.tcmd_test[k] = (tcmd_obs[k],tcmd_pvals[k])
    end

    # Return output structure
    return test_struct
    
end
##################################################################################################
## Matched test
##################################################################################################
"""
    `comp_mat_stat_tcmd(MODS_G1,MODS_G2)`

    Function that computes CMD between matched pairs in two groups (matched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_mat_stat_tcmd(mods_g1,mods_g2)
    ```
"""
function comp_mat_stat_tcmd(mods_g1::Vector{RegStruct}, mods_g2::Vector{RegStruct})::Vector{Float64}
    
    # Compute stat
    cmds = []
    @inbounds for s = 1:length(mods_g1)
        push!(cmds, comp_cmd(mods_g1[s], mods_g2[s]))
    end
    
    # Return mean CMD
    return sum(cmds) / length(mods_g1)

end
"""
    `mat_reg_test_tcmd(mods_g1,mods_g2)`

    Function that performs hypothesis testing in matched samples group comparison for Tcmd.

    # Examples
    ```julia-repl
    julia> CpelNano.mat_reg_test_tcmd(mods_g1,mods_g2)
    ```
"""
function mat_reg_test_tcmd(mods_g1::Vector{RegStruct}, mods_g2::Vector{RegStruct})::NTuple{2,Vector{Float64}}

    # Return stat-pval pair
    return (comp_mat_stat_tcmd(mods_g1, mods_g2), fill(NaN, length(mods_g1)))
    
end
"""
    `comp_mat_diff_mml(MMLs_G1,MMLs_G2)`

    Function that computes MML differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_mat_diff_mml(μs_g1,μs_g2)
    ```
"""
function comp_mat_diff_mml(μs_g1::Vector{Vector{Float64}}, μs_g2::Vector{Vector{Float64}})::Vector{Vector{Float64}}
    
    # Compute mean differences
    diffs = []
    @inbounds for s = 1:length(μs_g1)
        push!(diffs, μs_g1[s] - μs_g2[s])
    end

    # Return vector of differences
    return diffs

end
"""
    `comp_mat_diff_nme(NMEs_G1,NMEs_G2)`

    Function that computes NME differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_mat_diff_mml(hs_g1,hs_g2)
    ```
"""
function comp_mat_diff_nme(hs_g1::Vector{Vector{Float64}}, hs_g2::Vector{Vector{Float64}})::Vector{Vector{Float64}}
    
    # Compute entropy differences
    diffs = []
    @inbounds for s = 1:length(hs_g1)
        push!(diffs, hs_g1[s] - hs_g2[s])
    end

    # Return vector of differences
    return diffs

end
"""
    `comp_mat_j_stat(diffs,j)`

    Function that computes permutation statistic for j-th sign assigment given vector of differences
    between matched pairs (matched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_mat_j_stat(diffs,j)
    ```
"""
function comp_mat_j_stat(diffs::Vector{Vector{Float64}}, j::Int64)::Vector{Float64}
    
    # Init diff in permutation
    diffs_perm = fill(0.0, length(diffs[1]))
    
    # Change sign if pertinent
    @inbounds for k = 1:length(diffs[1])
        j_k = j
        @inbounds for s = 1:length(diffs)
            diffs_perm[k] += Bool(j_k & 1) ? diffs[s][k] : -diffs[s][k]
            j_k >>= 1
        end
    end

    # Return permutation statistic
    return diffs_perm / length(diffs)

end
"""
    `comp_mat_perm_stats(DIFFs,Js)`

    Function that computes permutation statistics given a vector of differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_mat_perm_stats(diffs,js)
    ```
"""
function comp_mat_perm_stats(diffs::Vector{Vector{Float64}}, js::Vector{Int64})::Vector{Vector{Float64}}
    
    # Return all possible signed sums
    return [comp_mat_j_stat(diffs, j) for j in js]

end
"""
    `mat_est_reg_test(MODELS_G1,MODELS_G2)`

    Function that performs hypothesis testing in unmatched samples group comparison in an estimation region.

    # Examples
    ```julia-repl
    julia> CpelNano.mat_est_reg_test(ms_g1,ms_g2)
    ```
"""
function mat_est_reg_test(ms_g1::Vector{RegStruct}, ms_g2::Vector{RegStruct})::RegStatTestStruct

    # Initialize output object
    test_struct = CpelNano.RegStatTestStruct(ms_g1[1])

    # Get stats for k-th analysis region
    μs_g1 = [m.mml for m in ms_g1]
    μs_g2 = [m.mml for m in ms_g2]
    hs_g1 = [m.nme for m in ms_g1]
    hs_g2 = [m.nme for m in ms_g2]

    # Compute observed differences
    mml_diffs = CpelNano.comp_mat_diff_mml(μs_g1, μs_g2)
    nme_diffs = CpelNano.comp_mat_diff_nme(hs_g1, hs_g2)
    
    # Compute observed stats
    tmml_obs = sum(mml_diffs) / length(ms_g1)
    tnme_obs = sum(nme_diffs) / length(ms_g1)
    tcmd_obs = CpelNano.comp_mat_stat_tcmd(ms_g1, ms_g2)

    # Init p-values
    tmml_pvals = fill(NaN, length(tmml_obs))
    tnme_pvals = fill(NaN, length(tnme_obs))
    
    # If enough data compute p-values
    if (0.5^(length(μs_g1)-1))<0.05
    
        # Get group label combinations to use
        exact = 2^length(ms_g1) < LMAX
        js = exact ? collect(0:2^length(ms_g1) - 1) : rand(0:2^length(ms_g1) - 1, LMAX)

        # Compute permutation stats
        tmml_perms = CpelNano.comp_mat_perm_stats(mml_diffs, js)
        tnme_perms = CpelNano.comp_mat_perm_stats(nme_diffs, js)
        
        # Compute p-values two-sided test
        @inbounds for k in 1:length(tmml_obs)

            # Check if data
            isnan(tmml_obs[k]) && continue

            # Get permutation stats from k-th analysis region
            tmml_perms_k = [perm[k] for perm in tmml_perms]
            tnme_perms_k = [perm[k] for perm in tnme_perms]

            # Get number of permutation stats equal or above observed
            tmml_pval_k = sum(abs.(tmml_perms_k) .>= abs(tmml_obs[k]))
            tnme_pval_k = sum(abs.(tnme_perms_k) .>= abs(tnme_obs[k]))
            
            # Compute p-values
            tmml_pvals[k] = exact ? tmml_pval_k / length(tmml_perms_k) : (1.0 + tmml_pval_k) / (1.0 + length(tmml_perms_k))
            tnme_pvals[k] = exact ? tnme_pval_k / length(tnme_perms_k) : (1.0 + tnme_pval_k) / (1.0 + length(tnme_perms_k))

        end

    end

    # Fill return object
    @inbounds for k in 1:length(tmml_pvals)
        test_struct.tests.tmml_test[k] = (tmml_obs[k], tmml_pvals[k])
        test_struct.tests.tnme_test[k] = (tnme_obs[k], tnme_pvals[k])
        test_struct.tests.tcmd_test[k] = (tcmd_obs[k], NaN)
    end

    # Return output structure
    return test_struct
    
end
##################################################################################################
## User functions
##################################################################################################
"""
    `read_in_grp_mods(MOD_FILES,CHR,FA_REC,CONFIG)`

    Function that reads in model files of group in chromosome `CHR`.

    # Examples
    ```julia-repl
    julia> CpelNano.read_in_grp_mods(mod_fls,chr,fa_rec,config)
    ```
"""
function read_in_grp_mods(mod_fls::Vector{String}, chr::String, fa_rec::FASTA.Record, config::CpelNanoConfig)::Vector{Dict{String,RegStruct}}
    
    # Read in models
    ms = Vector{Dict{String,RegStruct}}()
    for f in mod_fls
        push!(ms, read_mod_file_chr(f, chr, fa_rec, config))
    end

    # Return models
    return ms

end
"""
    `find_comm_regs(MODS_G1,MODS_G2)`

    Function that finds common regions between the two groups.

    # Examples
    ```julia-repl
    julia> CpelNano.find_comm_regs(ms_g1,ms_g2)
    ```
"""
function find_comm_regs(ms_g1::Vector{Dict{String,RegStruct}}, ms_g2::Vector{Dict{String,RegStruct}})::Vector{String}
    
    # Deal with group 1
    g1_dic = Dict{String,Int64}()
    for rep in ms_g1
        for k in keys(rep)
            g1_dic[k] = k in keys(g1_dic) ? g1_dic[k] + 1 : 1
        end
    end
    
    # Deal with group 2
    g2_dic = Dict{String,Int64}()
    for rep in ms_g2
        for k in keys(rep)
            g2_dic[k] = k in keys(g2_dic) ? g2_dic[k] + 1 : 1
        end
    end
    
    # Common regions
    comm_regs = Vector{String}()
    for k in keys(g1_dic)
        k in keys(g2_dic) && push!(comm_regs, k)
    end

    # Return common regions
    return comm_regs

end
"""
    `pmap_diff_grp_comp(REG_ID,MODS_G1,MODS_G2,CONFIG)`

    Function that performs differential methylation analysis between groups given two groups of 
    CpelNano model files at a given estimation region `REG_ID`.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_diff_grp_comp(reg_id,mods_g1,mods_g2,config)
    ```
"""
function pmap_diff_grp_comp(reg_id::String, ms_g1::Vector{Dict{String,RegStruct}}, ms_g2::Vector{Dict{String,RegStruct}}, config::CpelNanoConfig)::RegStatTestStruct

    
    ## Perform testing
    
    print_log("Working on: $(reg_id)")
    
    # Matched vs Unmatched
    if config.matched

        ## Get available pairs for estimation region
        reg_ms_g1 = Vector{RegStruct}()
        reg_ms_g2 = Vector{RegStruct}()
        for rep_ind in 1:length(ms_g1)
            if haskey(ms_g1[rep_ind],reg_id) && haskey(ms_g2[rep_ind],reg_id)
                # Get sample dictionaries
                rep_dic_g1 = ms_g1[rep_ind]
                rep_dic_g2 = ms_g2[rep_ind]
                # Push pair
                push!(reg_ms_g1, rep_dic_g1[reg_id])
                push!(reg_ms_g2, rep_dic_g2[reg_id])
            end
        end
        
        # Matched test
        out_test = mat_est_reg_test(reg_ms_g1, reg_ms_g2)

    else

        ## Get models for specific region
    
        # Group 1
        reg_ms_g1 = Vector{RegStruct}()
        for dict in ms_g1
            haskey(dict, reg_id) && push!(reg_ms_g1, dict[reg_id])
        end
    
        # Group 2
        reg_ms_g2 = Vector{RegStruct}()
        for dict in ms_g2
            haskey(dict, reg_id) && push!(reg_ms_g2, dict[reg_id])
        end
        
        ## Unmatched test
        out_test = unmat_est_reg_test(reg_ms_g1, reg_ms_g2)

    end

    # Return test struct
    return out_test

end
"""
    `diff_grp_comp(MODEL_FILES_G1,MODEL_FILES_G2,FASTA,CONFIG)`

    Function that performs differential methylation analysis between groups given two groups of 
    CpelNano model files.

    # Examples
    ```julia-repl
    julia> CpelNano.diff_grp_comp(mod_fls_g1,mod_fls_g2,fasta,config)
    ```
"""
function diff_grp_comp(mod_fls_g1::Vector{String}, mod_fls_g2::Vector{String}, fasta::String, config::CpelNanoConfig)::Nothing

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

        ms_g1 = read_in_grp_mods(mod_fls_g1, chr, fa_rec, config)
        ms_g2 = read_in_grp_mods(mod_fls_g2, chr, fa_rec, config)

        ## Find regions with at least one sample per group
        
        comm_regs = find_comm_regs(ms_g1, ms_g2)
        
        # If no data in chromosome continue
        length(comm_regs) > 0 || continue
        config.verbose && print_log("Comm regions: $(comm_regs)")

        ## Process each region
        pmap_out = pmap(x -> pmap_diff_grp_comp(x, ms_g1, ms_g2, config), comm_regs)

        ## Write chromosome output
        write_diff_output(pmap_out, config)

    end

    # Perform multiple hypothesis testing
    print_log("Multiple hypothesis testing")
    mult_hyp_corr(config)

    # Return nothing
    return nothing

end