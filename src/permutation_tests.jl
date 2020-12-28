##################################################################################################
## Unmatched test Tpdm analysis region
##################################################################################################
"""
    `comp_unmat_stat_pdm(mods_g1,mods_g2)`

    Function that computes test statistic Tpdm in unmatched case between groups.

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_unmat_stat_pdm(mods_g1::Vector{RegStruct},mods_g2::Vector{RegStruct})::Float64
    
    # Compute stat
    gjsds = []
    @inbounds for s1=1:length(mods_g1)
        @inbounds for s2=1:length(mods_g2)
            push!(gjsds,comp_gjsd(mods_g1[s1],mods_g2[s2]))
        end
    end

    # Return average GJSD between groups
    return sum(gjsds)/(length(mods_g1)*length(mods_g2))

end
"""
    `comp_unmat_perm_stat_pdm(mods_g1,mods_g2,perm_ids)`

    Function that produces a permutation statistic for Tpdm in unmatched case.

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_unmat_perm_stat_pdm(mods_g1::Vector{RegStruct},mods_g2::Vector{RegStruct},perm_ids::Vector{Int64})::Float64

    # Get vectors for each group
    mods_g2_p = vcat(mods_g1,mods_g2)
    mods_g1_p = mods_g2_p[perm_ids]
    deleteat!(mods_g2_p,perm_ids)

    # Return stat for permutation
    return comp_unmat_stat_pdm(mods_g1_p,mods_g2_p)
    
end
"""
    `unmat_reg_test_tpdm(mods_g1,mods_g2)`

    Function that performs hypothesis testing in unmatched samples group comparison for Tpdm.

    # Examples
    ```julia-repl
    julia> using Random; Random.seed!(1234);
    ```
"""
function unmat_reg_test_tpdm(mods_g1::Vector{RegStruct},mods_g2::Vector{RegStruct})::NTuple{2,Float64}

    # Compute observed stats
    tpdm_obs = comp_unmat_stat_pdm(mods_g1,mods_g2)
        # print_log("Tpdm=$(tpdm_obs)")

    # Compute number of possible randomizations
    L = binomial(length(mods_g1)+length(mods_g2),length(mods_g1))
    exact = L<LMAX

    # Create iteratable object with all combinations
    comb_iter = combinations(1:(length(mods_g1)+length(mods_g2)),length(mods_g1))

    # Get group label combinations to use
    comb_iter_used = []
    if exact
        # Use all group assignments
        comb_iter_used = comb_iter
    else
        # Use Lmax group assignments
        ind_subset = rand(1:L,LMAX)
        @inbounds for (ind,comb) in enumerate(comb_iter)
            (ind in ind_subset) && push!(comb_iter_used,comb)
        end
    end

    # Use method for random permutation
    tpdm_perms = map(x->comp_unmat_perm_stat_pdm(mods_g1,mods_g2,x),comb_iter_used)

    # Compute p-values two-sided test
    tpdm_pval = sum(tpdm_perms.>=tpdm_obs)
    tpdm_pval = exact ? tpdm_pval/length(tpdm_perms) : (1.0+tpdm_pval)/(1.0+length(tpdm_perms))

    # Return stat-pval pair
    return (tpdm_obs,tpdm_pval)
    
end
##################################################################################################
## Unmatched test Tmml & Tnme α-subregion
##################################################################################################
"""
    `comp_unmat_stat_mml(μs_g1,μs_g2)`

    Function that computes MML difference between groups (unmatched case).

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_unmat_stat_mml(μs_g1::Vector{Float64},μs_g2::Vector{Float64})::Float64
    
    # Compute means within groups
    mml1 = sum(μs_g1)/length(μs_g1)
    mml2 = sum(μs_g2)/length(μs_g2)

    # Return Tmml
    return mml1-mml2

end
"""
    `comp_unmat_perm_stat_mml(μs_g1,μs_g2,PERM_IDS)`

    Function that produces a permutation statistic for Tmml in unmatched case.

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_unmat_perm_stat_mml(μs_g1::Vector{Float64},μs_g2::Vector{Float64},perm_ids::Vector{Int64})::Float64

    # Get vectors for each group
    μs_g2p = vcat(μs_g1,μs_g2)
    μs_g1p = μs_g2p[perm_ids]
    deleteat!(μs_g2p,perm_ids)

    # Return
    return comp_unmat_stat_mml(μs_g1p,μs_g2p)
    
end
"""
    `comp_unmat_stat_nme(hs_g1,hs_g2)`

    Function that computes NME difference between groups (unmatched case).

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_unmat_stat_nme(hs_g1::Vector{Float64},hs_g2::Vector{Float64})::Float64
    
    # Compute h within group
    nme1 = sum(hs_g1)/length(hs_g1)
    nme2 = sum(hs_g2)/length(hs_g2)

    # Return difference in nme
    return nme1-nme2

end
"""
    `comp_unmat_perm_stat_nme(hs_g1,hs_g2,PERM_IDS)`

    Function that produces a permutation statistic for Tnme in unmatched case.

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_unmat_perm_stat_nme(hs_g1::Vector{Float64},hs_g2::Vector{Float64},perm_ids::Vector{Int64})::Float64

    # Get vectors for each group
    hs_g2p = vcat(hs_g1,hs_g2)
    hs_g1p = hs_g2p[perm_ids]
    deleteat!(hs_g2p,perm_ids)

    # Return
    return comp_unmat_stat_nme(hs_g1p,hs_g2p)
    
end
# """
#     `comp_unmat_stat_pdm(n,ϕ1s,ϕ2s,∇logZ1s,∇logZ2s)`

#     Function that computes test statistic Tgsd in unmatched case.

#     # Examples
#     ```julia-repl
#     julia> n=[5]; ϕ1s=[[-1.0,1.0],[-1.0,1.0]]; ϕ2s=[[1.0,1.0],[1.0,1.0]];
#     julia> ∇logZ1s = [CpelTdm.get_∇logZ(n,ϕ1) for ϕ1 in ϕ1s];
#     julia> ∇logZ2s = [CpelTdm.get_∇logZ(n,ϕ2) for ϕ2 in ϕ2s];
#     julia> CpelTdm.comp_unmat_stat_pdm(n,ϕ1s,ϕ2s,∇logZ1s,∇logZ2s)
#     0.7888652058295635
#     ```
# """
# function comp_unmat_stat_pdm(n::Vector{Int64},ϕ1s::Vector{Vector{Float64}},ϕ2s::Vector{Vector{Float64}},
#     ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}})::Float64
    
#     ## Compute stat
#     cmds = []
#     @inbounds for s1=1:length(ϕ1s)
#         @inbounds for s2=1:length(ϕ2s)
#             push!(cmds,comp_cmd(n,ϕ1s[s1],ϕ2s[s2],∇logZ1s[s1],∇logZ2s[s2]))
#         end
#     end

#     # Return sum of GJSD
#     return sum(cmds)/(length(ϕ1s)*length(ϕ2s))

# end
# """
#     `comp_unmat_perm_stat_pdm(n,ϕ1s,ϕ2s,∇logZ1s,∇logZ2s,perm_ids)`

#     Function that produces a permutation statistic for tpdm in unmatched case.

#     # Examples
#     ```julia-repl
#     julia> n=[10]; ϕ1s=fill([0.5,0.5],5); ϕ2s=fill([-0.5,0.5],5); perm_ids=[1,3,5,7,9];
#     julia> ∇logZ1s = [CpelTdm.get_∇logZ(n,ϕ1) for ϕ1 in ϕ1s];
#     julia> ∇logZ2s = [CpelTdm.get_∇logZ(n,ϕ2) for ϕ2 in ϕ2s];
#     julia> CpelTdm.comp_unmat_perm_stat_pdm(n,ϕ1s,ϕ2s,∇logZ1s,∇logZ2s,perm_ids)
#     0.1691594717040136
#     ```
# """
# function comp_unmat_perm_stat_pdm(n::Vector{Int64},ϕ1s::Vector{Vector{Float64}},ϕ2s::Vector{Vector{Float64}},
#     ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}},perm_ids::Vector{Int64})::Float64

#     # Get vectors for each group
#     ϕ2sp = vcat(ϕ1s,ϕ2s)
#     ϕ1sp = ϕ2sp[perm_ids]
#     deleteat!(ϕ2sp,perm_ids)

#     # Get vectors for each group
#     ∇logZ2sp = vcat(∇logZ1s,∇logZ2s)
#     ∇logZ1sp = ∇logZ2sp[perm_ids]
#     deleteat!(∇logZ2sp,perm_ids)  

#     # Return
#     return comp_unmat_stat_pdm(n,ϕ1sp,ϕ2sp,∇logZ1sp,∇logZ2sp)
    
# end
"""
    `unmat_subreg_test(MODELS_G1,MODELS_G2)`

    Function that performs hypothesis testing in unmatched samples group comparison in the k-th α-subregion.

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function unmat_subreg_test(ms_g1::Vector{RegStruct},ms_g2::Vector{RegStruct},k::Int64)::SubregStatTestStruct

    # Initialize output object
    out_struct = SubregStatTestStruct()

    # Get stats for k-th α-subregion
    μs_g1 = [m.proc ? m.mml[k] : NaN for m in ms_g1]
    μs_g2 = [m.proc ? m.mml[k] : NaN for m in ms_g2]
    hs_g1 = [m.proc ? m.nme_vec[k] : NaN for m in ms_g1]
    hs_g2 = [m.proc ? m.nme_vec[k] : NaN for m in ms_g2]

    # Compute observed stats
    tmml_obs = comp_unmat_stat_mml(μs_g1,μs_g2)
    tnme_obs = comp_unmat_stat_nme(hs_g1,hs_g2)
        # tpdm_obs = comp_unmat_stat_pdm(ms_g1,ms_g2)

    # Compute number of possible randomizations
    L = binomial(length(μs_g1)+length(μs_g2),length(μs_g1))
    exact = L<LMAX

    # Create iteratable object with all combinations
    comb_iter = combinations(1:(length(μs_g1)+length(μs_g2)),length(μs_g1))

    # Get group label combinations to use
    comb_iter_used = []
    if exact
        # Use all group assignments
        comb_iter_used = comb_iter
    else
        # Use Lmax group assignments
        ind_subset = rand(1:L,LMAX)
        @inbounds for (ind,comb) in enumerate(comb_iter)
            (ind in ind_subset) && push!(comb_iter_used,comb)
        end
    end

    # Use method for random permutation
    tmml_perms = map(x->comp_unmat_perm_stat_mml(μs_g1,μs_g2,x),comb_iter_used)
    tnme_perms = map(x->comp_unmat_perm_stat_nme(hs_g1,hs_g2,x),comb_iter_used)
        # tpdm_perms = map(x->comp_unmat_perm_stat_pdm(ms_g1,ms_g2,x),comb_iter_used)

    # Compute p-values two-sided test
    if exact
        # Compute exact p-value
        tmml_pval = sum(abs.(tmml_perms).>=abs(tmml_obs))/length(tmml_perms)
        tnme_pval = sum(abs.(tnme_perms).>=abs(tnme_obs))/length(tnme_perms)
            # tpdm_pval = sum(tpdm_perms.>=tpdm_obs)/length(tpdm_perms)
    else
        # Compute p-value in MC setting
        tmml_pval = (1.0+sum(abs.(tmml_perms).>=abs(tmml_obs)))/(1.0+length(tmml_perms))
        tnme_pval = (1.0+sum(abs.(tnme_perms).>=abs(tnme_obs)))/(1.0+length(tnme_perms))
            # tpdm_pval = (1.0+sum(tpdm_perms.>=tpdm_obs))/(1.0+length(tpdm_perms))
    end

    # Fill return object
    out_struct.tmml_test = (tmml_obs,tmml_pval)
    out_struct.tnme_test = (tnme_obs,tnme_pval)
        # out_struct.tpdm_test = (tpdm_obs,tpdm_pval)

    # Set as processed
    out_struct.proc = true

    # Return output structure
    return out_struct
    
end
##################################################################################################
## Matched test
##################################################################################################
"""
    `comp_mat_stat_tpdm(n,ϕ1s,ϕ2s)`

    Function that computes GJSD between pairs (matched case).

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_mat_stat_tpdm(mods_g1::Vector{RegStruct},mods_g2::Vector{RegStruct})::Float64
    
    # Compute stat
    gjsd = 0.0
    @inbounds for s=1:length(mods_g1)
        gjsd += comp_gjsd(mods_g1[s],mods_g2[s])
    end
    
    # Return mean GJSD
    return gjsd/length(mods_g1)

end
"""
    `mat_reg_test_tpdm(mods_g1,mods_g2)`

    Function that performs hypothesis testing in matched samples group comparison for Tpdm.

    # Examples
    ```julia-repl
    julia> using Random; Random.seed!(1234);
    ```
"""
function mat_reg_test_tpdm(mods_g1::Vector{RegStruct},mods_g2::Vector{RegStruct})::NTuple{2,Float64}

    # Return stat-pval pair
    return (comp_mat_stat_tpdm(mods_g1,mods_g2),NaN)
    
end
##################################################################################################
## Unmatched test Tmml & Tnme α-subregion
##################################################################################################
"""
    `comp_mat_diff_mml(MMLs_G1,MMLs_G2)`

    Function that computes MML differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_mat_diff_mml(μs_g1::Vector{Float64},μs_g2::Vector{Float64})::Vector{Float64}
    
    # Compute mean differences
    diffs = []
    @inbounds for s=1:length(μs_g1)
        push!(diffs,μs_g1[s]-μs_g2[s])
    end

    # Return vector of differences
    return diffs

end
"""
    `comp_mat_diff_nme(NMEs_G1,NMEs_G2)`

    Function that computes NME differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_mat_diff_nme(hs_g1::Vector{Float64},hs_g2::Vector{Float64})::Vector{Float64}
    
    # Compute entropy differences
    diffs = []
    @inbounds for s=1:length(hs_g1)
        push!(diffs,hs_g1[s]-hs_g2[s])
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
    julia> 
    ```
"""
function comp_mat_j_stat(diffs::Vector{Float64},j::Int64)::Float64
    
    # Init diff in permutation
    diffs_perm = 0.0
    
    # Change sign if pertinent
    @inbounds for i=1:length(diffs)
        diffs_perm += Bool(j & 1) ? diffs[i] : -diffs[i]
        j >>= 1
    end

    # Return permutation statistic
    return diffs_perm/length(diffs)

end
"""
    `comp_mat_perm_stats(n,ϕ1s,ϕ2s)`

    Function that computes permutation statistics given a vector of differences between pairs (matched case).

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function comp_mat_perm_stats(diffs::Vector{Float64},js::Vector{Int64})::Vector{Float64}
    
    # Return all possible signed sums
    return [comp_mat_j_stat(diffs,j) for j in js]

end
# """
#     `comp_mat_stat_pdm(n,ϕ1s,ϕ2s)`

#     Function that computes GJSD between pairs (matched case).

#     # Examples
#     ```julia-repl
#     julia> n=[5]; ϕ1s=[[-1.0,1.0],[-1.0,1.0]]; ϕ2s=[[1.0,1.0],[1.0,1.0]];
#     julia> ∇logZ1s = [CpelTdm.get_∇logZ(n,ϕ1) for ϕ1 in ϕ1s];
#     julia> ∇logZ2s = [CpelTdm.get_∇logZ(n,ϕ2) for ϕ2 in ϕ2s];
#     julia> CpelTdm.comp_mat_stat_pdm(n,ϕ1s,ϕ2s,∇logZ1s,∇logZ2s)
#     0.7888652058295635
#     ```
# """
# function comp_mat_stat_pdm(n::Vector{Int64},ϕ1s::Vector{Vector{Float64}},ϕ2s::Vector{Vector{Float64}},
#     ∇logZ1s::Vector{Vector{Float64}},∇logZ2s::Vector{Vector{Float64}})::Float64
    
#     # Compute stat
#     cmds = 0.0
#     @inbounds for s=1:length(ϕ1s)
#         cmds += comp_cmd(n,ϕ1s[s],ϕ2s[s],∇logZ1s[s],∇logZ2s[s])
#     end

#     # Return mean GJSD
#     return cmds/length(ϕ1s)

# end
"""
    `mat_subreg_test(MODELS_G1,MODELS_G2,K)`

    Function that performs hypothesis testing in unmatched samples group comparison in the k-th α-subregion.

    # Examples
    ```julia-repl
    julia> 
    ```
"""
function mat_subreg_test(ms_g1::Vector{RegStruct},ms_g2::Vector{RegStruct},k::Int64)::SubregStatTestStruct

    # Initialize output object
    out_struct = SubregStatTestStruct()

    # Get stats for k-th α-subregion
    μs_g1 = [m.mml[k] for m in ms_g1]
    μs_g2 = [m.mml[k] for m in ms_g2]
    hs_g1 = [m.nme_vec[k] for m in ms_g1]
    hs_g2 = [m.nme_vec[k] for m in ms_g2]

    # Compute observed stats
    mml_diffs = comp_mat_diff_mml(μs_g1,μs_g2)
    nme_diffs = comp_mat_diff_nme(hs_g1,hs_g2)
        # tpdm_obs = comp_mat_stat_pdm(ms_g1,ms_g2)

    # Get group label combinations to use
    exact = 2^length(ms_g1)<LMAX
    js = exact ? collect(0:2^length(ms_g1)-1) : rand(0:2^length(ms_g1)-1,LMAX)

    # Compute permutation stats
    tmml_perms = comp_mat_perm_stats(mml_diffs,js)
    tnme_perms = comp_mat_perm_stats(nme_diffs,js)

    # Compute observed stats
    tmml_obs = sum(mml_diffs)/length(ms_g1)
    tnme_obs = sum(nme_diffs)/length(ms_g1)
    
    # Compute p-values
    if exact
        # Compute exact p-value
        tmml_pval = sum(abs.(tmml_perms).>=abs(tmml_obs))/length(tmml_perms)
        tnme_pval = sum(abs.(tnme_perms).>=abs(tnme_obs))/length(tnme_perms)
    else
        # Compute p-value in MC setting
        tmml_pval = (1.0+sum(abs.(tmml_perms).>=abs(tmml_obs)))/(1.0+length(tmml_perms))
        tnme_pval = (1.0+sum(abs.(tnme_perms).>=abs(tnme_obs)))/(1.0+length(tnme_perms))
    end

    # Fill return object
    out_struct.tmml_test = (tmml_obs,tmml_pval)
    out_struct.tnme_test = (tnme_obs,tnme_pval)
        # out_struct.tpdm_test = (tpdm_obs,tpdm_pval)

    # Set as processed
    out_struct.proc = true

    # Return output structure
    return out_struct
    
end