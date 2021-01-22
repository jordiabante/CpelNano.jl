##################################################################################################
## Unmatched test Tcmd analysis region
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
            # TODO: consider computing only min{s1,s2} pairs
            push!(cmds, comp_cmd(mods_g1[s1], mods_g2[s2]))
        end
    end

    # Return vector of average CMD between groups
    return sum(cmds) / (length(mods_g1) * length(mods_g2))

end
"""
    `comp_unmat_perm_stat_cmd(mods_g1,mods_g2,perm_ids)`

    Function that produces a permutation statistic for Tcmd in unmatched case.

    # Examples
    ```julia-repl
    julia> CpelNano.comp_unmat_perm_stat_cmd(mods_g1,mods_g2,perm_ids)
    ```
"""
function comp_unmat_perm_stat_cmd(mods_g1::Vector{RegStruct}, mods_g2::Vector{RegStruct}, perm_ids::Vector{Int64})::Vector{Float64}

    # Get vectors for each group
    mods_g2_p = vcat(mods_g1, mods_g2)
    mods_g1_p = mods_g2_p[perm_ids]
    deleteat!(mods_g2_p, perm_ids)

    # Return stat for permutation
    return comp_unmat_stat_cmd(mods_g1_p, mods_g2_p)
    
end
"""
    `unmat_reg_test_tcmd(mods_g1,mods_g2)`

    Function that performs hypothesis testing in unmatched group comparison for Tcmd.

    # Examples
    ```julia-repl
    julia> CpelNano.unmat_reg_test_tcmd(mods_g1,mods_g2)
    ```
"""
function unmat_reg_test_tcmd(mods_g1::Vector{RegStruct}, mods_g2::Vector{RegStruct})::NTuple{2,Vector{Float64}}

    # Compute observed stats
    tcmd_obs = CpelNano.comp_unmat_stat_cmd(mods_g1, mods_g2)
        # print_log("Tcmd=$(tcmd_obs)")

    # Compute number of possible randomizations
    L = binomial(length(mods_g1) + length(mods_g2), length(mods_g1))
    exact = L < LMAX

    # Create iteratable object with all combinations
    comb_iter = combinations(1:(length(mods_g1) + length(mods_g2)), length(mods_g1))

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

    # Use method for random permutation. TODO: consider parallelization here
    tcmd_perms = map(x -> CpelNano.comp_unmat_perm_stat_cmd(mods_g1, mods_g2, x), comb_iter_used)

    # Compute p-values two-sided test
    tcmd_pvals = fill(NaN, length(tcmd_obs))
    @inbounds for k in 1:length(tcmd_pvals)
        perms_k = [perm[k] for perm in tcmd_perms]
        pval_k = sum(perms_k .>= tcmd_obs[k])
        tcmd_pvals[k] = exact ? pval_k / length(perms_k) : (1.0 + pval_k) / (1.0 + length(perms_k))
    end
    
    # Return stat-pval pairs
    return (tcmd_obs, tcmd_pvals)
    
end
##################################################################################################
## Unmatched test Tmml & Tnme
##################################################################################################
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

    Function that produces a permutation statistic for Tmml in unmatched case.

    # Examples
    ```julia-repl
    julia> 
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
    julia> 
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

    Function that produces a permutation statistic for Tnme in unmatched case.

    # Examples
    ```julia-repl
    julia> 
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

    Function that performs hypothesis testing in unmatched samples group comparison in the k-th analysis region.

    # Examples
    ```julia-repl
    julia> 
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

    # Compute number of possible randomizations
    L = binomial(length(μs_g1) + length(μs_g2), length(μs_g1))
    exact = L < LMAX

    # Create iteratable object with all combinations
    comb_iter = combinations(1:(length(μs_g1) + length(μs_g2)), length(μs_g1))

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

    # Use method for random permutation
    tmml_perms = map(x -> comp_unmat_perm_stat_mml(μs_g1, μs_g2, x), comb_iter_used)
    tnme_perms = map(x -> comp_unmat_perm_stat_nme(hs_g1, hs_g2, x), comb_iter_used)
    tcmd_perms = map(x -> comp_unmat_perm_stat_cmd(ms_g1, ms_g2, x), comb_iter_used)

    # Compute p-values two-sided test
    tmml_pvals = fill(NaN, length(tmml_obs))
    tnme_pvals = fill(NaN, length(tnme_obs))
    tcmd_pvals = fill(NaN, length(tcmd_obs))
    @inbounds for k in 1:length(tmml_obs)

        # Check if data
        isnan(tmml_obs[k]) && continue
        
        # Set processed as true
        test_struct.tests[k].proc = true

        # Get permutation stats from k-th analysis region
        tmml_perms_k = [perm[k] for perm in tmml_perms]
        tnme_perms_k = [perm[k] for perm in tnme_perms]
        tcmd_perms_k = [perm[k] for perm in tcmd_perms]

        # Get number of permutation stats equal or above observed
        tmml_pval_k = sum(abs.(tmml_perms_k) .>= abs.(tmml_obs[k]))
        tnme_pval_k = sum(abs.(tnme_perms_k) .>= abs.(tnme_obs[k]))
        tcmd_pval_k = sum(tcmd_perms_k .>= tcmd_obs[k])
        
        # Compute p-values
        tmml_pvals[k] = exact ? tmml_pval_k / length(tmml_perms_k) : (1.0 + tmml_pval_k) / (1.0 + length(tmml_perms_k))
        tnme_pvals[k] = exact ? tnme_pval_k / length(tnme_perms_k) : (1.0 + tnme_pval_k) / (1.0 + length(tnme_perms_k))
        tcmd_pvals[k] = exact ? tcmd_pval_k / length(tcmd_perms_k) : (1.0 + tcmd_pval_k) / (1.0 + length(tcmd_perms_k))

    end
    
    # Fill return object
    @inbounds for k in 1:length(tmml_pvals)
        test_struct.tests[k].proc || continue
        test_struct.tests[k].tmml_test = (tmml_obs[k], tmml_pvals[k])
        test_struct.tests[k].tnme_test = (tnme_obs[k], tnme_pvals[k])
        test_struct.tests[k].tcmd_test = (tcmd_obs[k], tcmd_pvals[k])
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
##################################################################################################
## Matched test Tmml & Tnme analysis region
##################################################################################################
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
    `mat_est_reg_test(MODELS_G1,MODELS_G2,K)`

    Function that performs hypothesis testing in unmatched samples group comparison in the k-th analysis region.

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

    # Compute observed stats
    mml_diffs = CpelNano.comp_mat_diff_mml(μs_g1, μs_g2)
    nme_diffs = CpelNano.comp_mat_diff_nme(hs_g1, hs_g2)
    tcmd_obs = CpelNano.comp_mat_stat_tcmd(ms_g1, ms_g2)

    # Get group label combinations to use
    exact = 2^length(ms_g1) < LMAX
    js = exact ? collect(0:2^length(ms_g1) - 1) : rand(0:2^length(ms_g1) - 1,LMAX)
        # js = collect(0:2^length(ms_g1) - 1)

    # Compute permutation stats
    tmml_perms = CpelNano.comp_mat_perm_stats(mml_diffs, js)
    tnme_perms = CpelNano.comp_mat_perm_stats(nme_diffs, js)

    # Compute observed stats
    tmml_obs = sum(mml_diffs) / length(ms_g1)
    tnme_obs = sum(nme_diffs) / length(ms_g1)
    
    # Compute p-values two-sided test
    tmml_pvals = fill(NaN, length(tmml_obs))
    tnme_pvals = fill(NaN, length(tnme_obs))
    @inbounds for k in 1:length(tmml_obs)

        # Check if data
        isnan(tmml_obs[k]) && continue
        
        # Set processed as true
        test_struct.tests[k].proc = true

        # Get permutation stats from k-th analysis region
        tmml_perms_k = [perm[k] for perm in tmml_perms]
        tnme_perms_k = [perm[k] for perm in tnme_perms]

        # Get number of permutation stats equal or above observed
        tmml_pval_k = sum(abs.(tmml_perms_k) .>= abs.(tmml_obs[k]))
        tnme_pval_k = sum(abs.(tnme_perms_k) .>= abs.(tnme_obs[k]))
        
        # Compute p-values
        tmml_pvals[k] = exact ? tmml_pval_k / length(tmml_perms_k) : (1.0 + tmml_pval_k) / (1.0 + length(tmml_perms_k))
        tnme_pvals[k] = exact ? tnme_pval_k / length(tnme_perms_k) : (1.0 + tnme_pval_k) / (1.0 + length(tnme_perms_k))


    end

    # Fill return object
    @inbounds for k in 1:length(tmml_pvals)
        test_struct.tests[k].proc || continue
        test_struct.tests[k].tmml_test = (tmml_obs[k], tmml_pvals[k])
        test_struct.tests[k].tnme_test = (tnme_obs[k], tnme_pvals[k])
        test_struct.tests[k].tcmd_test = (tcmd_obs[k], NaN)
    end

    # Return output structure
    return test_struct
    
end