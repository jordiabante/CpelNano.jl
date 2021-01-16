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
    `get_Z!(REG_DATA)`

    Computes and stores partition function.

    # Examples
    ```julia-repl
    julia> CpelNano.get_Z!(x)
    ```
"""
function get_Z!(rs::RegStruct)::Nothing

    # Store partition function
    rs.Z = get_Z(rs.tm.u1,rs.tm.uN,rs.tm.Ws)

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
    `get_log_gs!(REG_DATA)`

    Computes and stores log(g_i(±⋅)) values.

    # Examples
    ```julia-repl
    julia> CpelNano.get_log_gs!(rs)
    ```
"""
function get_log_gs!(rs::RegStruct)::Nothing

    # Loop over analysis regions
    for k=1:rs.nls_rgs.num

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

    `get_exp_log_g1!(RS_DATA)`
    
    Computes E[log(g1(Xp))] for all analysis regions using the transfer matrix method and log-trick.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_exp_log_g1!(rs)
    ```
"""
function get_exp_log_g1!(rs::RegStruct)::Nothing

    # Loop over analysis regions
    rs.exps.log_g1 = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num
        # Get index
        p = minimum(rs.nls_rgs.cpg_ind[k])
        # Continue if no CpG sites
        p==0 && continue
        # Compute expectation
        if p==1
            # If first CpG site (can happen k>1)
            rs.exps.log_g1[k] = 0.0
        else
            # Compute
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

    `get_exp_log_g2!(RS_DATA)`
    
    Computes E[log(g2(Xq))] for all analysis regions using the transfer matrix method and log-trick.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_exp_log_g2!(rs)
    ```
"""
function get_exp_log_g2!(rs::RegStruct)::Nothing

    # Loop over analysis regions
    rs.exps.log_g2 = fill(NaN,rs.nls_rgs.num)
    @inbounds for k=1:rs.nls_rgs.num
        # Get index
        q = maximum(rs.nls_rgs.cpg_ind[k])
        # Continue if no CpG sites
        q==0 && continue
        # Compute expectation
        if q==rs.N
            # If last CpG site (can happen k<K)
            rs.exps.log_g2[k] = 0.0
        else
            # Compute
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
    for k=1:rs.nls_rgs.num
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

    Function that computes normalized methylation entropy vector (NMEV).

    # Examples
    ```julia-repl
    julia> CpelNano.comp_nme!(rs)
    ```
"""
function comp_nme!(rs::RegStruct)::Nothing

    # Get αs & βs
    αs,βs = get_αβ_from_ϕ(rs.ϕhat,rs)

    # Loop over subregions
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
    `comp_gjsd(REG_DATA_1,REG_DATA_2)`

    Function that computes the normalized generalised Jensen-Shannon divergence over
    the entire analysis region.

    # Examples
    ```julia-repl
    julia> 𝒜s=[1:20]; ℬs=[1:19]; ϕ=[0.0,0.0]; x=CpelNano.RegStruct([],ϕ,𝒜s,ℬs);
    julia> CpelNano.get_∇logZ!(x); CpelNano.get_rs_mats!(x); CpelNano.get_Z!(x); CpelNano.comp_nme_vec!(x);
    julia> CpelNano.comp_gjsd(x,x)
    0.0
    julia> 𝒜s=[1:20]; ℬs=[1:19]; ϕ1=[-3.75,1.5]; ϕ2=[3.75,1.5];
    julia> x1=CpelNano.RegStruct([],ϕ1,𝒜s,ℬs); x2=CpelNano.RegStruct([],ϕ2,𝒜s,ℬs);
    julia> CpelNano.get_∇logZ!(x1); CpelNano.get_rs_mats!(x1); CpelNano.get_Z!(x1); CpelNano.comp_nme_vec!(x1);
    julia> CpelNano.get_∇logZ!(x2); CpelNano.get_rs_mats!(x2); CpelNano.get_Z!(x2); CpelNano.comp_nme_vec!(x2);
    julia> CpelNano.comp_gjsd(x1,x2)
    0.9993863729172315
    ```
"""
function comp_gjsd(dat1::RegStruct,dat2::RegStruct)::Float64

    # Parameter geometric mixture of Ising
    αs1,βs1 = get_αβ_from_ϕ(dat1.ϕhat,dat1)
    αs2,βs2 = get_αβ_from_ϕ(dat2.ϕhat,dat2)
    ϕγ = 0.5.*(dat1.ϕhat+dat2.ϕhat)
    αsγ = 0.5.*(αs1+αs2)
    βsγ = 0.5.*(βs1+βs2)
    
    # Get partition function
    logu1 = get_log_u(αsγ[1])
    loguN = get_log_u(αsγ[end])
    logWs = [get_log_W(αsγ[n],αsγ[n+1],βsγ[n]) for n=1:length(βsγ)]
    logZγ = get_log_Z(logu1,loguN,logWs)
        # print_log("log(Z1)=$(log(dat1.Z)); log(Z2)=$(log(dat2.Z)); log(Zγ)=$(logZγ)")

    # Compute numerator of 1-GJSD
    gjsd = log(dat1.Z) + log(dat2.Z) - (dat1.ϕhat'*dat1.∇logZ + dat2.ϕhat'*dat2.∇logZ)

    # Normalize 1-GJSD
    gjsd /= 2*logZγ - ϕγ'*(dat1.∇logZ + dat2.∇logZ)
    
    # Return GJSD
    return min(1.0,max(0.0,1.0-gjsd))
    
end
"""
    `comp_gjsd_vec(REG_DATA_1,REG_DATA_2)`

    Function that computes the normalized generalised Jensen-Shannon divergence (G-JSD) over
    each α-subregion.

    # Examples
    ```julia-repl
    julia> 𝒜s=[1:20]; ℬs=[1:19]; ϕ=[0.0,0.0]; x=CpelNano.RegStruct([],ϕ,𝒜s,ℬs);
    julia> CpelNano.get_∇logZ!(x); CpelNano.get_rs_mats!(x); CpelNano.get_Z!(x); CpelNano.comp_nme_vec!(x);
    julia> CpelNano.comp_gjsd_vec(x,x)
    0.0
    julia> 𝒜s=[1:20]; ℬs=[1:19]; ϕ1=[-1.75,1.25]; ϕ2=[1.75,1.25];
    julia> x1=CpelNano.RegStruct([],ϕ1,𝒜s,ℬs); x2=CpelNano.RegStruct([],ϕ2,𝒜s,ℬs);
    julia> CpelNano.get_∇logZ!(x1); CpelNano.get_rs_mats!(x1); CpelNano.get_Z!(x1); CpelNano.comp_nme_vec!(x1);
    julia> CpelNano.get_∇logZ!(x2); CpelNano.get_rs_mats!(x2); CpelNano.get_Z!(x2); CpelNano.comp_nme_vec!(x2);
    julia> CpelNano.comp_gjsd_vec(x1,x2)
    0.9674329433697847
    julia> 𝒜s=[1:5,6:10,11:15,16:20]; ℬs=[1:9,10:19];
    julia> ϕ1=vcat(fill(1.0,4),[1.0,1.0]); ϕ2=vcat(fill(-1.0,4),[1.0,1.0]);
    julia> x1=CpelNano.RegStruct([],ϕ1,𝒜s,ℬs); x2=CpelNano.RegStruct([],ϕ2,𝒜s,ℬs);
    julia> CpelNano.get_∇logZ!(x1); CpelNano.get_rs_mats!(x1); CpelNano.get_Z!(x1); 
    julia> CpelNano.get_log_gs!(x1); CpelNano.comp_mml!(x1); CpelNano.comp_nme_vec!(x1);
    julia> CpelNano.get_∇logZ!(x2); CpelNano.get_rs_mats!(x2); CpelNano.get_Z!(x2); 
    julia> CpelNano.get_log_gs!(x2); CpelNano.comp_mml!(x2); CpelNano.comp_nme_vec!(x2);
    julia> CpelNano.comp_gjsd_vec(x1,x2)
    4-element Array{Float64,1}:
     0.8506403090777066
     0.9156156533053169
     0.9156156533053222
     0.8506403090776973
    ```
"""
function comp_gjsd_vec(dat1::RegStruct,dat2::RegStruct)::Vector{Float64}

    ## Initialize
    
    # Get αs & βs
    K1 = length(dat1.𝒜s)
    gjsd_vec = fill(NaN,K1)
    
    # Parameter geometric mixture of Ising
    αs1,βs1 = get_αβ_from_ϕ(dat1.ϕhat,dat1)
    αs2,βs2 = get_αβ_from_ϕ(dat2.ϕhat,dat2)
    ϕγ = 0.5.*(dat1.ϕhat+dat2.ϕhat)
    αsγ = 0.5.*(αs1+αs2)
    βsγ = 0.5.*(βs1+βs2)

    # Create struct for geometric mixture
    mix = RegStruct([],ϕγ,dat1.𝒜s,dat1.ℬs)
    get_rs_mats!(mix)
    get_Z!(mix)
    get_log_gs!(mix)

    ## Get E[X] and E[XX] under ϕ1 and ϕ2

    # Scale Ws for numerical stability
    Ws1 = US * dat1.tm.Ws
    Ws2 = US * dat2.tm.Ws

    # Get E[X]
    ex1 = get_E_X(dat1.tm.u1,dat1.tm.uN,Ws1,dat1.Z)
    ex2 = get_E_X(dat2.tm.u1,dat2.tm.uN,Ws2,dat2.Z)
    
    # Get E[XX]
    exx1 = get_E_XX(dat1.tm.u1,dat1.tm.uN,Ws1,dat1.Z)
    exx2 = get_E_XX(dat2.tm.u1,dat2.tm.uN,Ws2,dat2.Z)

    # Compute scaled Z's
    Z1s = dat1.tm.u1' * prod(Ws1) * dat1.tm.uN
    Z2s = dat2.tm.u1' * prod(Ws2) * dat2.tm.uN

    # Loop over subregions
    @inbounds for k=1:K1

        # Get p and q
        p = minimum(mix.𝒜s[k])
        q = maximum(mix.𝒜s[k])

        # Create V matrix for ϕγ
        Vp = [mix.tm.log_gs.pm[k] 0.0; 0.0 mix.tm.log_gs.pp[k]]
        Vq = [mix.tm.log_gs.qm[k] 0.0; 0.0 mix.tm.log_gs.qp[k]]

        ## Deal with ϕ1

        # First term
        if p==1
            aux1 = 0.0
        elseif p==dat1.L
            aux1 = dat1.tm.u1' * prod(Ws1) * Vp * dat1.tm.uN / Z1s
        else
            aux1 = dat1.tm.u1' * prod(Ws1[1:(p-1)]) * Vp * prod(Ws1[p:end]) * dat1.tm.uN / Z1s
        end
            # aux1 = p==1 ? 0.0 : dat1.tm.u1' * prod(Ws1[1:(p-1)]) * Vp * prod(Ws1[p:end]) * 
            #     dat1.tm.uN / Z1s

        # Second term
        if q==dat1.L
            aux2 = 0.0
        elseif q==1
            aux2 = dat1.tm.u1' * Vq * prod(Ws1) * dat1.tm.uN / Z1s
        else
            aux2 = dat1.tm.u1' * prod(Ws1[1:(q-1)]) * Vq * prod(Ws1[q:end]) * dat1.tm.uN / Z1s
        end
            # aux2 = q==dat1.L ? 0.0 : dat1.tm.u1' * prod(Ws1[1:(q-1)]) * Vq * prod(Ws1[q:end]) *
            #     dat1.tm.uN / Z1s

        # Cross entropy
        cross = log(mix.Z)-aux1-aux2-αsγ[p:q]'*ex1[p:q]
        cross += p==q ? 0.0 : -βsγ[p:(q-1)]'*exx1[p:(q-1)]

        ## Deal with ϕ2

        # Deal with first term
        if p==1
            aux1 = 0.0
        elseif p==dat2.L
            aux1 = dat2.tm.u1' * prod(Ws2) * Vp * dat2.tm.uN / Z2s
        else
            aux1 = dat2.tm.u1' * prod(Ws2[1:(p-1)]) * Vp * prod(Ws2[p:end]) * dat2.tm.uN / Z2s
        end
            # aux1 = p==1 ? 0.0 : dat2.tm.u1' * prod(Ws2[1:(p-1)]) * Vp * prod(Ws2[p:end]) * 
            #     dat2.tm.uN / Z2s

        # Second term
        if q==dat2.L
            aux2 = 0.0
        elseif q==1
            aux2 = dat2.tm.u1' * Vq * prod(Ws2) * dat2.tm.uN / Z2s
        else
            aux2 = dat2.tm.u1' * prod(Ws2[1:(q-1)]) * Vq * prod(Ws2[q:end]) * dat2.tm.uN / Z2s
        end
            # aux2 = q==dat2.L ? 0.0 : dat2.tm.u1' * prod(Ws2[1:(q-1)]) * Vq * prod(Ws2[q:end]) * 
            #     dat2.tm.uN / Z2s
        
        # Cross-entropy
        cross += log(mix.Z)-aux1-aux2-αsγ[p:q]'*ex2[p:q]
        cross += p==q ? 0.0 : -βsγ[p:(q-1)]'*exx2[p:(q-1)]

        # Add G-JSD component
        gjsd_vec[k] = (dat1.nme_vec[k]+dat2.nme_vec[k])*(q-p+1)*LOG2
        gjsd_vec[k] /= cross
        
    end

    # Return GJSD
    return max.(0.0,min.(1.0,1.0.-gjsd_vec))
    
end
###################################################################################################
# DEVELOPING QUANTITIES
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