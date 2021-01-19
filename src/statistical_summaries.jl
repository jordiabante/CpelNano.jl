"""
    `get_rs_mats!(REG_DATA)`

    Computes necessary matrices for transfer matrix methods.

    # Examples
    ```julia-repl
    julia> CpelNano.get_rs_mats!(x)
    ```
"""
function get_rs_mats!(rs::RegStruct)::Nothing

    # Get αs & βs
    αs,βs = get_αβ_from_ϕ(rs.ϕhat,rs)

    # Get u's
    rs.tm.u1 = get_u(αs[1])
    rs.tm.uN = get_u(αs[end])

    # Get W's
    Ws = Vector{Array{Float64,2}}()
    @inbounds for n=1:length(βs)
        push!(Ws,get_W(αs[n],αs[n+1],βs[n]))
    end

    # Record W's
    rs.tm.Ws = Ws

    # Return nothing
    return nothing

end
"""
    `get_rs_lgtrck_mats!(REG_DATA)`

    Computes necessary matrices for transfer matrix methods that relies on the log-sum trick.

    # Examples
    ```julia-repl
    julia> CpelNano.get_rs_lgtrck_mats!(x)
    ```
"""
function get_rs_lgtrck_mats!(rs::RegStruct)::Nothing

    # Get αs & βs
    αs,βs = get_αβ_from_ϕ(rs.ϕhat,rs)

    # Get u's
    rs.tm.log_u1 = get_log_u(αs[1])
    rs.tm.log_uN = get_log_u(αs[end])

    # Get W's
    logWs = Vector{Array{Float64,2}}()
    @inbounds for n=1:length(βs)
        push!(logWs,get_log_W(αs[n],αs[n+1],βs[n]))
    end

    # Record W's
    rs.tm.log_Ws = logWs

    # Return nothing
    return nothing

end
"""
    `get_rs_logZ!(REG_DATA)`

    Computes and stores the (log) partition function using log-sum trick transfer matrix.

    # Examples
    ```julia-repl
    julia> CpelNano.get_rs_logZ!(x)
    ```
"""
function get_rs_logZ!(rs::RegStruct)::Nothing

    # Store partition function
    rs.logZ = get_log_Z(rs.tm.log_u1,rs.tm.log_uN,rs.tm.log_Ws)

    # Return nothing
    return nothing

end
"""
    `comp_log_g_p(REG_DATA,P,Xp)`

    Computes log(g_1(xp)) values.

    # Examples
    ```julia-repl
    julia> CpelNano.comp_log_g_p(rs,p,xp)
    ```
"""
function comp_log_g_p(rs::RegStruct,p::Int64,xp::Float64)::Float64

    # If p=1 then return nothing
    p==1 && return 0.00

    # Get αs & βs
    αs,βs = get_αβ_from_ϕ(rs.ϕhat,rs)

    # Get log(uN)
    aux = (αs[p-1]+xp*βs[p-1])/2.0
    log_uNp = [-aux,+aux]

    # If p=2 then return inner product of v's
    p==2 && return log_vec_vec_mult(log_uNp,log_uNp)

    # Get log(W)
    log_Wp = [
        -αs[p-2]/2.0+βs[p-2]-aux -αs[p-2]/2.0-βs[p-2]+aux;
        +αs[p-2]/2.0-βs[p-2]-aux +αs[p-2]/2.0+βs[p-2]+aux
    ]

    # If p=3 then return inner product of v's
    p==3 && return log_vec_vec_mult(log_vec_mat_mult(rs.tm.log_u1,log_Wp),log_uNp)

    # Return log(g1(xp))
    log_Ws = mult_log_mats(vcat(rs.tm.log_Ws[1:(p-3)],[log_Wp]))
    return log_vec_vec_mult(log_vec_mat_mult(rs.tm.log_u1,log_Ws),log_uNp)

end
"""
    `comp_log_g_q(REG_DATA,Q,Xq)`

    Computes log(g_2(xq)) values.

    # Examples
    ```julia-repl
    julia> CpelNano.comp_log_g_q(rs,q,xq)
    ```
"""
function comp_log_g_q(rs::RegStruct,q::Int64,xq::Float64)::Float64

    # If q=N(=L after re-scaling) then return log(1)
    q==rs.L && return 0.00

    # Get αs & βs
    αs,βs = get_αβ_from_ϕ(rs.ϕhat,rs)
    
    # Get log(u1)
    aux = (αs[q+1]+xq*βs[q])/2.0
    log_u1p = [-aux; +aux]

    # If q=N-1 then return inner product of v's
    q==rs.L-1 && return log_vec_vec_mult(log_u1p,log_u1p)

    # Get Wp
    log_Wp = [
        -aux+βs[q+1]-αs[q+2]/2.0 -aux-βs[q+1]+αs[q+2]/2.0;
        +aux-βs[q+1]-αs[q+2]/2.0 +aux+βs[q+1]+αs[q+2]/2.0
    ]

    # If q=N-2 then return inner product of v's
    q==rs.L-2 && return log_vec_vec_mult(log_vec_mat_mult(log_u1p,log_Wp),rs.tm.log_uN)

    # Return log(g2(xq))
    log_Ws = mult_log_mats(vcat([log_Wp],rs.tm.log_Ws[(q+2):end]))
    return log_vec_vec_mult(log_vec_mat_mult(log_u1p,log_Ws),rs.tm.log_uN)

end
"""
    `get_rs_log_gs!(REG_DATA)`

    Computes and stores log(g_i(±⋅)) values.

    # Examples
    ```julia-repl
    julia> CpelNano.get_rs_log_gs!(rs)
    ```
"""
function get_rs_log_gs!(rs::RegStruct)::Nothing

    # Loop over analysis regions
    @inbounds for k=1:rs.nls_rgs.num

        # Get p and q
        p = minimum(rs.nls_rgs.cpg_ind[k])
        q = maximum(rs.nls_rgs.cpg_ind[k])

        if p!=0 && q!=0
            # Compute log g1(±x_p) & log g2(±x_q)
            push!(rs.tm.log_gs.pm,comp_log_g_p(rs,p,-1.0))
            push!(rs.tm.log_gs.pp,comp_log_g_p(rs,p,+1.0))
            push!(rs.tm.log_gs.qm,comp_log_g_q(rs,q,-1.0))
            push!(rs.tm.log_gs.qp,comp_log_g_q(rs,q,+1.0))
        else
            # Push NaN
            push!(rs.tm.log_gs.pm,NaN)
            push!(rs.tm.log_gs.pp,NaN)
            push!(rs.tm.log_gs.qm,NaN)
            push!(rs.tm.log_gs.qp,NaN)
        end

    end

    # Return nothing
    return nothing

end
"""

    `get_rs_exp_log_g1!(RS_DATA)`
    
    Computes E[log(g1(Xp))] for all analysis regions using the transfer matrix method and log-trick.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_rs_exp_log_g1!(rs)
    ```
"""
function get_rs_exp_log_g1!(rs::RegStruct)::Nothing

    # Loop over analysis regions
    rs.exps.log_g1 = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num
        # Get index & next if no CpG sites
        p = minimum(rs.nls_rgs.cpg_ind[k])
        p==0 && continue
        # Compute expectation
        if p==1
            # If first CpG site (can happen k>1)
            rs.exps.log_g1[k] = 0.0
        else
            # Compute expectation using TM (using marginal seems to be unstable)
            logD = log.([rs.tm.log_gs.pm[k] 0.0;0.0 rs.tm.log_gs.pp[k]])
            log_WsDWs = mult_log_mats(vcat(rs.tm.log_Ws[1:(p-1)],[logD],rs.tm.log_Ws[p:end]))
            log_u1WsDWs = log_vec_mat_mult(rs.tm.log_u1,log_WsDWs)
            rs.exps.log_g1[k] = exp(log_vec_vec_mult(log_u1WsDWs,rs.tm.log_uN)-rs.logZ)
        end
    end

    # Return nothing
    return nothing

end
"""

    `comp_crs_exps_log_g1(RS,RS_MIX)`
    
    Computes E[log(g1(Xp;θ̃));θi] for all analysis regions using the transfer matrix method and log-trick.
    
    # Examples
    ```julia-repl
    julia> CpelNano.comp_crs_exps_log_g1(rs,rs_mix)
    ```
"""
function comp_crs_exps_log_g1(rs::RegStruct,rs_mix::RegStruct)::Vector{Float64}

    # Loop over analysis regions
    exps_log_g1 = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num
        # Get index & next if no CpG sites
        p = minimum(rs.nls_rgs.cpg_ind[k])
        p==0 && continue
        # Compute expectation
        if p==1
            # If first CpG site (can happen k>1)
            exps_log_g1[k] = 0.0
        else
            # Compute expectation using TM (using marginal seems to be unstable)
            logD = log.([rs_mix.tm.log_gs.pm[k] 0.0;0.0 rs_mix.tm.log_gs.pp[k]])
            log_WsDWs = mult_log_mats(vcat(rs.tm.log_Ws[1:(p-1)],[logD],rs.tm.log_Ws[p:end]))
            log_u1WsDWs = log_vec_mat_mult(rs.tm.log_u1,log_WsDWs)
            exps_log_g1[k] = exp(log_vec_vec_mult(log_u1WsDWs,rs.tm.log_uN)-rs.logZ)
        end
    end

    # Return vector
    return exps_log_g1

end
"""

    `get_rs_exp_log_g2!(RS_DATA)`
    
    Computes E[log(g2(Xq))] for all analysis regions using the transfer matrix method and log-trick.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_rs_exp_log_g2!(rs)
    ```
"""
function get_rs_exp_log_g2!(rs::RegStruct)::Nothing

    # Loop over analysis regions
    rs.exps.log_g2 = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num
        # Get index & next if no CpG sites
        q = maximum(rs.nls_rgs.cpg_ind[k])
        q==0 && continue
        # Compute expectation
        if q==rs.N
            # If last CpG site (can happen k<K)
            rs.exps.log_g2[k] = 0.0
        else
            # Compute expectation using TM (using marginal seems to be unstable)
            logD = log.([rs.tm.log_gs.qm[k] 0.0;0.0 rs.tm.log_gs.qp[k]])
            log_WsDWs = mult_log_mats(vcat(rs.tm.log_Ws[1:(q-1)],[logD],rs.tm.log_Ws[q:end]))
            log_u1WsDWs = log_vec_mat_mult(rs.tm.log_u1,log_WsDWs)
            rs.exps.log_g2[k] = exp(log_vec_vec_mult(log_u1WsDWs,rs.tm.log_uN)-rs.logZ)
        end
    end

    # Return nothing
    return nothing

end
"""

    `comp_crs_exps_log_g2(RS,RS_MIX)`
    
    Computes E[log(g2(Xq;θ̃));θi] for all analysis regions using the transfer matrix method and log-trick.
    
    # Examples
    ```julia-repl
    julia> CpelNano.comp_crs_exps_log_g2(rs,rs_mix)
    ```
"""
function comp_crs_exps_log_g2(rs::RegStruct,rs_mix::RegStruct)::Vector{Float64}

    # Loop over analysis regions
    exps_log_g2 = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num
        # Get index & next if no CpG sites
        q = maximum(rs.nls_rgs.cpg_ind[k])
        q==0 && continue
        # Compute expectation
        if q==rs.N
            # If last CpG site (can happen k<K)
            exps_log_g2[k] = 0.0
        else
            # Compute expectation using TM (using marginal seems to be unstable)
            logD = log.([rs_mix.tm.log_gs.qm[k] 0.0;0.0 rs_mix.tm.log_gs.qp[k]])
            log_WsDWs = mult_log_mats(vcat(rs.tm.log_Ws[1:(q-1)],[logD],rs.tm.log_Ws[q:end]))
            log_u1WsDWs = log_vec_mat_mult(rs.tm.log_u1,log_WsDWs)
            exps_log_g2[k] = exp(log_vec_vec_mult(log_u1WsDWs,rs.tm.log_uN)-rs.logZ)
        end
    end

    # Return vector
    return exps_log_g2

end
"""
    `comp_mml!(REG_DATA)`

    Function that computes mean methylation level (MML) vector.

    # Examples
    ```julia-repl
    julia> CpelNano.comp_mml!(rs) 
    ```
"""
function comp_mml!(rs::RegStruct)::Nothing

    # Loop over analysis regions
    rs.mml = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num

        # Check analysis region has CpG sites
        p = minimum(rs.nls_rgs.cpg_ind[k])
        p==0 && continue
        
        # Check number of CpG sites
        Nk = length(rs.nls_rgs.cpg_ind[k])
        
        # Set MML vector
        rs.mml[k] = 0.5/Nk * (Nk + sum(rs.exps.ex[rs.nls_rgs.cpg_ind[k]]))

    end

    # Return nothing
    return nothing

end
"""
    `comp_nme!(REG_DATA)`

    Function that computes normalized methylation entropy (NME) vector for each analysis 
    region in the estimation region.

    # Examples
    ```julia-repl
    julia> CpelNano.comp_nme!(rs)
    ```
"""
function comp_nme!(rs::RegStruct)::Nothing

    # Get αs & βs
    αs,βs = get_αβ_from_ϕ(rs.ϕhat,rs)

    # Loop over analysis regions
    rs.nme = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num

        # Get indices
        p = minimum(rs.nls_rgs.cpg_ind[k])
        q = maximum(rs.nls_rgs.cpg_ind[k])
        p==q==0 && continue

        # Add terms and normalize
        nme = rs.logZ - rs.exps.log_g1[k] - rs.exps.log_g2[k] - αs[p:q]'*rs.exps.ex[p:q]
        nme += p==q ? 0.0 : - βs[p:(q-1)]' * rs.exps.exx[p:(q-1)]
        nme /= (q-p+1)*LOG2
        
        # Store
        rs.nme[k] = nme
        
    end

    # Return nothing
    return nothing

end
"""
    `comp_crs_ent(RS,RS_MIX)`

    Function that computes cross-entropy H(θi,θ̃). This function is used in function
    `comp_cmd` to compute the Coefficient of Methylation Divergence (CMD).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_crs_ent(rs,rs_mix)
    ```
"""
function comp_crs_ent(rs::RegStruct,rs_mix::RegStruct)::Vector{Float64}

    # Get αs & βs of geometric mixture of Isings
    α̃s,β̃s = get_αβ_from_ϕ(rs_mix.ϕhat,rs_mix)

    # Compute E[log g1(Xp;θ̃);θi] and E[log g2(Xq;θ̃);θi] for all analysis regions
    exps_log_g1 = comp_crs_exps_log_g1(rs,rs_mix)
    exps_log_g2 = comp_crs_exps_log_g2(rs,rs_mix)
        # print_log("exps_log_g1=$(exps_log_g1)")
        # print_log("exps_log_g2=$(exps_log_g2)")

    # Loop over analysis regions
    crs_ent_vec = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num

        # Get indices
        p = minimum(rs.nls_rgs.cpg_ind[k])
        q = maximum(rs.nls_rgs.cpg_ind[k])
        p==q==0 && continue

        # Compute cross entropy
        crs_ent = rs_mix.logZ - α̃s[p:q]'*rs.exps.ex[p:q] - exps_log_g1[k] - exps_log_g2[k]
        crs_ent += p==q ? 0.0 : - β̃s[p:(q-1)]' * rs.exps.exx[p:(q-1)]
        
        # Store
        crs_ent_vec[k] = crs_ent
        
    end

    # Return vector
    return crs_ent_vec

end
"""
    `comp_cmd(REG_DATA_1,REG_DATA_2)`

    Function that computes the normalized generalised Jensen-Shannon divergence (G-JSD), or
    Coefficient of Methylation Divergence (CMD), over each analysis region.

    # Examples
    ```julia-repl
    julia> comp_cmd(REG_DATA_1,REG_DATA_2)
    ```
"""
function comp_cmd(rs1::RegStruct,rs2::RegStruct)::Vector{Float64}

    ## Struct for geometric mixture
    
    # Create structure
    mix = CpelNano.RegStruct()
    mix.cpg_pos = rs1.cpg_pos; 
    mix.N = rs1.N; 
    mix.ρn = rs1.ρn; 
    mix.dn = rs1.dn; 
    mix.nls_rgs = rs1.nls_rgs;
    mix.Nl = rs1.Nl;
    mix.L = rs1.L;
    mix.ρl = rs1.ρl;
    mix.dl = rs1.dl;
    mix.ϕhat = 0.5.*(rs1.ϕhat+rs2.ϕhat)

    # Required quantities
    get_rs_lgtrck_mats!(mix)
    get_rs_logZ!(mix)
    get_rs_exps!(mix)
    get_rs_log_gs!(mix)

    ## Compute cross-entropies
    
    # Cross-entropy w/ θ1 & θ2
    cross_1 = comp_crs_ent(rs1,mix)
    cross_2 = comp_crs_ent(rs2,mix)
    
    ## Compute CMD for each analysis region
    
    cmd_vec = fill(NaN,rs1.nls_rgs.num)
    @inbounds for k=1:rs1.nls_rgs.num

        # Get indices
        p = minimum(rs1.nls_rgs.cpg_ind[k])
        q = maximum(rs1.nls_rgs.cpg_ind[k])

        # Skip if empty analysis region
        p==q==0 && continue
        
        # Shannon's Entropy
        ent_1 = rs1.nme[k] * (q-p+1) * LOG2
        ent_2 = rs2.nme[k] * (q-p+1) * LOG2
            # print_log("H(θ1)=$(ent_1); H(θ2)=$(ent_2)")
            # print_log("H(θ1,θ̃)=$(cross_1[k]); H(θ2,θ̃)=$(cross_2[k])")

        # If CMD is undefined continue
        cross_1[k]==cross_2[k]==0.0 && continue

        # Store CMD k-th component
        cmd_vec[k] = 1.0 - (ent_1 + ent_2) / (cross_1[k] + cross_2[k]) 

    end

    # Return CMD
    return cmd_vec
    
end
###################################################################################################
# OTHER QUANTITIES
###################################################################################################
"""
    `get_ρ(u1, uN, Ws, Z)`
    
    Computes Pearson correlation between x_n and x_(n+1)
    
    # Examples
    ```julia-repl
    julia> Ngrp = 10; αs = fill(0.0,Ngrp); βs = fill(0.0,Ngrp-1);
    julia> u1 = CpelNano.get_u(αs[1]); uN = CpelNano.get_u(αs[end]);
    julia> Ws = [CpelNano.get_W(αs[n],αs[n+1],βs[n]) for n=1:length(βs)];
    julia> Z = CpelNano.get_Z(u1,uN,Ws);
    julia> CpelNano.get_ρ(u1, uN, Ws, Z)
    9-element Array{Float64,1}:
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
    0.0
 ```
 """
function get_ρ(u1::Vector{Float64},uN::Vector{Float64},Ws::Vector{Array{Float64,2}},Z::Float64)::Vector{Float64}

    # Get E[X] & E[X{i}X{i+1}]
    exs = CpelNano.get_E_X(u1,uN,Ws,Z)
    exxs = CpelNano.get_E_XX(u1,uN,Ws,Z)
    
    # Compute √Var(Xi)
    std_vec = sqrt.(ones((length(Ws)+1))-exs.^2) 

    # Compute covariance
    ex = exs[1:length(exs)-1]
    ex_prime = exs[2:length(exs)]
    cov_vec = exxs - (ex .* ex_prime)
    
    # Return correlation
    return cov_vec ./ (std_vec[1:length(std_vec)-1] .* std_vec[2:length(std_vec)])

end

"""
    `get_E_XiXj(u1,uN,Ws,Z,i,j)`
    
    Computes E[X_iX_j], where j is larger than i
    
    # Examples
    ```julia-repl
    julia> Ngrp = 10; αs = fill(0.0,Ngrp); βs = fill(0.0,Ngrp-1);
    julia> u1 = CpelNano.get_u(αs[1]); uN = CpelNano.get_u(αs[end]);
    julia> Ws = [CpelNano.get_W(αs[n],αs[n+1],βs[n]) for n=1:length(βs)];
    julia> Z = CpelNano.get_Z(u1,uN,Ws);
    julia> CpelNano.get_E_XiXj(u1,uN,Ws,Z,1,10)
    0.0
 ```
 """
function get_E_XiXj(u1::Vector{Float64},uN::Vector{Float64},Ws::Vector{Array{Float64,2}},Z::Float64,i::Int64,j::Int64)::Float64
    
    # Scale Ws for numerical stability
    #Ws = US*Ws

    # Confirm that j is larger than i
    if i >= j
        return NaN
    end

    ## Compute E[X_iX_j]

    # Compute transfer matrix product up to i
    if i == 1
        exixj_temp = u1' * [-1 0; 0 1]
    else
        exixj_temp = u1' * prod(Ws[1:i-1]) * [-1 0; 0 1]
    end

    # Compute transfer matrix product from i to N
    if j == (length(Ws)+1)
        exixj = exixj_temp * prod(Ws[i:end]) * [-1 0; 0 1] * uN
    else
        exixj = exixj_temp * prod(Ws[i:j-1]) * [-1 0; 0 1] * prod(Ws[j:end]) * uN
    end
    
    # Return E[X_iX_j]
    return exixj / (u1' * prod(Ws) * uN)
    
end

"""
    `Σmat(u1,uN,Ws,Z)`
    
    Computes the covariance matrix in upper triangular form, where the diagonal is the variance and the [i,j] element 
    is the covariance between Xi and Xj
    
    # Examples
    ```julia-repl
    julia> Ngrp = 10; αs = fill(0.0,Ngrp); βs = fill(0.0,Ngrp-1);
    julia> u1 = CpelNano.get_u(αs[1]); uN = CpelNano.get_u(αs[end]);
    julia> Ws = [CpelNano.get_W(αs[n],αs[n+1],βs[n]) for n=1:length(βs)];
    julia> Z = CpelNano.get_Z(u1,uN,Ws);
    julia> CpelNano.Σmat(u1,uN,Ws,Z)
    10×10 UpperTriangular{Float64,Array{Float64,2}}:
     1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     ⋅   1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     ⋅    ⋅   1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
     ⋅    ⋅    ⋅   1.0  0.0  0.0  0.0  0.0  0.0  0.0
     ⋅    ⋅    ⋅    ⋅   1.0  0.0  0.0  0.0  0.0  0.0
     ⋅    ⋅    ⋅    ⋅    ⋅   1.0  0.0  0.0  0.0  0.0
     ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0  0.0  0.0  0.0
     ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0  0.0  0.0
     ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0  0.0
     ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅   1.0
 ```
 """
function Σmat(u1::Vector{Float64},uN::Vector{Float64},Ws::Vector{Array{Float64,2}},Z::Float64)::UpperTriangular

    # Num CpG sites
    N = length(Ws)+1

    # Get E[X]
    exs = CpelNano.get_E_X(u1,uN,Ws,Z)
    
    # Compute Var(Xi)
    var_vec = ones(N)-exs.^2 #add square
    
    # Construct covariance matrix
    cov_matrix = Array{Float64}(undef,N,N)

    @inbounds for i in 1:N
        @inbounds for j in i:N
            # Diagonal elements are the variance
            if i == j
                cov_matrix[i,j] = var_vec[i] 
            else
                cov_matrix[i,j] = CpelNano.get_E_XiXj(u1,uN,Ws,Z,i,j)
            end
        end
    end

    # Return Σ
    return UpperTriangular(cov_matrix)

end