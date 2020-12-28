##################################################################################################
## Transfer matrix method functions for conditional expectations E[⋅|y]
##################################################################################################
"""

    `get_uc(α,call)`
    
    Generates vector u1 and uN used in transfer matrix method for p(x̄|y).
    
    # Examples
    ```julia-repl
    julia> call = CpelNano.MethCallCpgGrp(-50.0,-45.0);
    julia> CpelNano.get_uc(0.0,call)
    2-element Array{Float64,1}:
     1.3887943864964021e-11
     1.6918979226151304e-10
    ```
"""
function get_uc(α::Float64,call::MethCallCpgGrp)::Vector{Float64}
    
    # Return u for CpG passed
    return [exp(get_γc(false,α,call)); exp(get_γc(true,α,call))]

end
"""

    `get_Wc(α1,α2,β,call1,call2)`
    
    Generates transition matrix W used in transfer matrix method for p(x̄|y).
    
    # Examples
    ```julia-repl
    julia> call = CpelNano.MethCallCpgGrp(-50.0,-45.0);
    julia> Wc = CpelNano.get_Wc(0.0,0.0,0.0,call,call)
    2×2 Array{Float64,2}:
     1.92875e-22  2.3497e-21 
     2.3497e-21   2.86252e-20
    ```
"""
function get_Wc(α1::Float64,α2::Float64,β::Float64,call1::MethCallCpgGrp,call2::MethCallCpgGrp)::Array{Float64,2}
    
    # Store intermediate quantities
    γ1_f = get_γc(false,α1,call1)
    γ1_t = get_γc(true,α1,call1)
    γ2_f = get_γc(false,α2,call2)
    γ2_t = get_γc(true,α2,call2)
    exp_γ1_f = exp(γ1_f)
    exp_γ1_t = exp(γ1_t)
    exp_γ2_f = exp(γ2_f)
    exp_γ2_t = exp(γ2_t)
    exp_β = exp(β)

    # Define transition matrix
    W = zeros(2,2)
    W[1,1] = exp_γ1_f * exp_γ2_f * exp_β
    W[2,1] = exp_γ1_t * exp_γ2_f / exp_β
    W[1,2] = exp_γ1_f * exp_γ2_t / exp_β
    W[2,2] = exp_γ1_t * exp_γ2_t * exp_β

    # Return W
    return W

end
"""

    `get_Zc(u1,uL,Ws)`
    
    Computes partition function of p(x̄_m|y_m;ϕ) using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 2; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-100.0,-100.0),L);
    julia> u1 = CpelNano.get_uc(αs[1],obs[1]); uL = CpelNano.get_uc(αs[end],obs[end]);
    julia> Ws = [CpelNano.get_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> CpelNano.get_Zc(u1,uL,Ws)
    5.53558610694695e-87
    ```
"""
function get_Zc(u1::Vector{Float64},uL::Vector{Float64},Ws::Vector{Array{Float64,2}})::Float64
   
    # Return u_1⊺ W_1 ⋯ W_L u_L for m-th observation
    return u1' * prod(Ws) * uL 

end
"""

    `get_Ec_X(u1,uL,Ws)`
    
    Computes E[X̄|y_m;ϕ] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 2; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-175.0,-170.0),L);
    julia> u1 = CpelNano.get_uc(αs[1],obs[1]); uL = CpelNano.get_uc(αs[end],obs[end]);
    julia> Ws = [CpelNano.get_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> Zc = CpelNano.get_Zc(u1,uL,Ws);
    julia> CpelNano.get_Ec_X(u1,uL,Ws,Zc)
    2-element Array{Float64,1}:
     0.9866142981514304
     0.9866142981514304
    ```
"""
function get_Ec_X(u1::Vector{Float64},uL::Vector{Float64},Ws::Vector{Array{Float64,2}},Zc::Float64)::Vector{Float64}
   
    # Compute E[X̄_l|y_m] for l=1,2,…,L
    prod_Ws = prod(Ws)
    exs = fill(NaN,length(Ws)+1)
    @inbounds for l=1:length(exs)
        if l==1
            exs[1] = u1' * D_EX * prod_Ws * uL
        elseif l==length(exs)
            exs[end] = u1' * prod_Ws * D_EX * uL
        else
            exs[l] = u1' * prod(Ws[1:(l-1)]) * D_EX * prod(Ws[l:end]) * uL
        end
    end

    # Return vector of Ec[X̄_l|y_m] for m-th observation
    return exs ./ Zc 

end
"""

    `get_Ec_XX(u1,uL,Ws)`
    
    Computes E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-50.0,-50.0),L);
    julia> u1 = CpelNano.get_uc(αs[1],obs[1]); uL = CpelNano.get_uc(αs[end],obs[end]);
    julia> Ws = [CpelNano.get_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> Zc = CpelNano.get_Zc(u1,uL,Ws);
    julia> CpelNano.get_Ec_XX(u1,uL,Ws,Zc)
    9-element Array{Float64,1}:
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
function get_Ec_XX(u1::Vector{Float64},uL::Vector{Float64},Ws::Vector{Array{Float64,2}},Zc::Float64)::Vector{Float64}
   
    # Init output
    exxs = fill(NaN,length(Ws))
    
    # If only one W
    if length(Ws)==1 
        exxs[1] = u1'*(D_EXX.* Ws[1]) * uL
        return exxs ./ Zc 
    end

    # Compute E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] for l=1,2,…,L-1
    @inbounds for l=1:length(exxs)
        if l==1
            exxs[1] = u1' * (D_EXX .* Ws[1]) * prod(Ws[2:end]) * uL
        elseif l==length(Ws)
            exxs[end] = u1' * prod(Ws[1:(end-1)]) * (Ws[end] .* D_EXX) * uL
        else
            exxs[l] = u1' * prod(Ws[1:(l-1)]) * (Ws[l] .* D_EXX) * prod(Ws[(l+1):end]) * uL
        end
    end

    # Return vector of E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] for m-th observation
    return exxs ./ Zc 

end
"""

    `get_cond_exs(αs,βs,obs)`
    
    Computes E[X̄_{m,l}|Y_m=y_m;ϕ] & E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] using the transfer matrix method. This function 
    does not use the log-sum trick, however. Thus, it is not as numerically stable as `get_cond_exs_log()`.
    
    # Examples
    ```julia-repl
    julia> L = 5; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-50.0,-45.0),L);
    julia> ex,exx = CpelNano.get_cond_exs(αs,βs,obs);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-45.0,-50.0),L);
    julia> ex,exx = CpelNano.get_cond_exs(αs,βs,obs);
    ```
"""
function get_cond_exs(αs::Vector{Float64},βs::Vector{Float64},obs::Vector{MethCallCpgGrp})::NTuple{2,Vector{Float64}}

    # Get u's
    u1 = get_uc(αs[1],obs[1])
    uL = get_uc(αs[end],obs[end])

    # Get W's
    Ws = [get_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:length(βs)]
    # print_log("$(prod(Ws))")

    # Scale Ws for numerical stability
    Ws = UniformScaling(1.0/maximum(maximum.(Ws))) * Ws
    # print_log("Scaled: $(prod(Ws))")
    # readline()

    # Get Zc
    Zc = get_Zc(u1,uL,Ws)

    # Ex
    ex = get_Ec_X(u1,uL,Ws,Zc)

    # Exx
    exx = get_Ec_XX(u1,uL,Ws,Zc)

    # Return
    return ex,exx
    
end
"""

    `get_Ec_logpyx(log_u1,log_uN,log_Ws,log_Zc)`
    
    Computes E[log p(y|x)|ỹ_m] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-50.0,-50.0),L);
    julia> u1 = CpelNano.get_uc(αs[1],obs[1]); 
    julia> uL = CpelNano.get_uc(αs[end],obs[end]);
    julia> Ws = [CpelNano.get_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> Zc = CpelNano.get_Zc(u1,uL,Ws);
    julia> CpelNano.get_Ec_logpyx_log(logu1,loguL,logWs,logZc,obs)
    
    ```
"""
function get_Ec_logpyx(u1::Vector{Float64},uL::Vector{Float64},Ws::Vector{Array{Float64,2}},Zc::Float64,call_m::Vector{MethCallCpgGrp})::Vector{Float64}
    
    # Compute E[log p(y|x)|ỹ]
    prod_Ws = prod(Ws)
    ex_log_pyxs = zeros(length(call_m))
    @inbounds for l=1:length(ex_log_pyxs)

        # Skip if no data
        call_m[l].obs || continue

        # Diagonal matrix to introduce
        Dl = [call_m[l].log_pyx_u 0.0; 0.0 call_m[l].log_pyx_m]

        # Get l-th term
        if l==1
            ex_log_pyxs[1] = u1' * Dl * prod_Ws * uL
        elseif l==length(ex_log_pyxs)
            ex_log_pyxs[end] = u1' * prod_Ws * Dl * uL
        else
            ex_log_pyxs[l] = u1' * prod(Ws[1:(l-1)]) * Dl * prod(Ws[l:end]) * uL
        end

    end

    # Return E[log p(y|x)|ỹ]
    return ex_log_pyxs ./ Zc 

end
##################################################################################################
## Transfer matrix method functions for expectations E[X̄] && E[X̄X̄]
##################################################################################################
"""

    `get_u(α)`
    
    Generates vector u1 and uL used in transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_u(0.0)
    2-element Array{Float64,1}:
     1.0
     1.0
    ```
"""
function get_u(α::Float64)::Vector{Float64}
    
    # Compute intermediate
    aux = exp(0.5 * α)

    # Return u for CpG passed
    return [1.0/aux; aux]

end
"""

    `get_W(α1,α2,β)`
    
    Generates transition matrix W used in transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_W(0.0,0.0,0.0)
    2×2 Array{Float64,2}:
     1.0  1.0
     1.0  1.0
    ```
"""
function get_W(α1::Float64,α2::Float64,β::Float64)::Array{Float64,2}
    
    # Intermediate quantities
    exp_α1 = exp(0.5*α1)
    exp_α2 = exp(0.5*α2)
    exp_β = exp(β)

    # Define transition matrix
    W = zeros(2,2)
    W[1,1] = exp_β / (exp_α1 * exp_α2)
    W[2,1] = exp_α1 / (exp_α2 * exp_β)
    W[1,2] = exp_α2 / (exp_α1 * exp_β)
    W[2,2] = exp_α1 * exp_α2 * exp_β

    # Return W
    return W

end
"""

    `get_Z(u1,uL,Ws)`
    
    Computes partition function of p(x̄;ϕ) using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> u1 = CpelNano.get_u(αs[1]); uL = CpelNano.get_u(αs[end]);
    julia> Ws = [CpelNano.get_W(αs[l],αs[l+1],βs[l]) for l=1:length(βs)];
    julia> CpelNano.get_Z(u1,uL,Ws)
    1024.0
    ```
"""
function get_Z(u1::Vector{Float64},uL::Vector{Float64},Ws::Vector{Array{Float64,2}})::Float64
   
    # Return u_1⊺ W_1 ⋯ W_L uL of p(x̄;ϕ)
    return u1' * prod(Ws) * uL 

end
"""

    `get_E_X(u1,uN,Ws,Z)`
    
    Computes E[X̄;ϕ] using the transfer matrix method. CAUTION: This function is not as numerically 
    stable as `get_E_X_log(logu1,loguN,logWs,logZ)`, albeit significantly faster.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> u1 = CpelNano.get_u(αs[1]); uL = CpelNano.get_u(αs[end]);
    julia> Ws = [CpelNano.get_W(αs[l],αs[l+1],βs[l]) for l=1:length(βs)];
    julia> Ws = UniformScaling(1.0/maximum(maximum.(Ws))) * Ws
    julia> Z = CpelNano.get_Z(u1,uL,Ws);
    julia> CpelNano.get_E_X(u1,uL,Ws,Z)
    10-element Array{Float64,1}:
     0.0
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
function get_E_X(u1::Vector{Float64},uL::Vector{Float64},Ws::Vector{Array{Float64,2}},Z::Float64)::Vector{Float64}

    # Compute E[X̄_l;ϕ] for l=1,2,…,L
    prod_Ws = prod(Ws)
    exs = fill(NaN,length(Ws)+1)
    @inbounds for l=1:length(exs)
        if l==1
            exs[1] = u1' * D_EX * prod_Ws * uL
        elseif l==length(exs)
            exs[end] = u1' * prod_Ws * D_EX * uL
        else
            exs[l] = u1' * prod(Ws[1:(l-1)]) * D_EX * prod(Ws[l:end]) * uL
        end
    end

    # Return vector of E[X̄_l;ϕ] for l=1,2,…,L
    return exs ./ Z

end
"""

    `get_E_XX(u1,uL,Ws,Z)`
    
    Computes E[X̄_{l}X̄_{l+1};ϕ] using the transfer matrix method. CAUTION: This function is not as numerically 
    stable as `get_E_XX_log(logu1,loguN,logWs,logZ)`, albeit significantly faster.
    
    # Examples
    ```julia-repl
    julia> using LinearAlgebra
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> u1 = CpelNano.get_u(αs[1]); uL = CpelNano.get_u(αs[end]);
    julia> Ws = [CpelNano.get_W(αs[l],αs[l+1],βs[l]) for l=1:length(βs)];
    julia> Z = CpelNano.get_Z(u1,uL,Ws);
    julia> Ws = UniformScaling(1.0/maximum(maximum.(Ws))) * Ws
    julia> CpelNano.get_E_XX(u1,uL,Ws,Z)
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
function get_E_XX(u1::Vector{Float64},uL::Vector{Float64},Ws::Vector{Array{Float64,2}},Z::Float64)::Vector{Float64}

    # Init output
    exxs = fill(NaN,length(Ws))
    
    # If only one W
    if length(Ws)==1 
        exxs[1] = u1'*(D_EXX.* Ws[1]) * uL
        return exxs ./ Z
    end

    # Compute E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] for l=1,2,…,L-1
    @inbounds for l=1:length(exxs)
        if l==1
            exxs[1] = u1' * (D_EXX .* Ws[1]) * prod(Ws[2:end]) * uL
        elseif l==length(Ws)
            exxs[end] = u1' * prod(Ws[1:(end-1)]) * (Ws[end] .* D_EXX) * uL
        else
            exxs[l] = u1' * prod(Ws[1:(l-1)]) * (Ws[l] .* D_EXX) * prod(Ws[(l+1):end]) * uL
        end
    end

    # Return vector of E[X̄_{l}X̄_{l+1};ϕ] for l=1,2,…,L-1
    return exxs ./ Z

end
"""

    `get_marg_px(u1,uL,Ws,Z)`
    
    Computes marginal probability P(X̄_l=1;ϕ) using the transfer matrix method. CAUTION: This function is not as 
    numerically stable as `get_marg_px_log(logu1,loguN,logWs,logZ)`, albeit significantly faster.
    
    # Examples
    ```julia-repl
    julia> using LinearAlgebra
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> u1 = CpelNano.get_u(αs[1]); 
    julia> uL = CpelNano.get_u(αs[end]);
    julia> Ws = [CpelNano.get_W(αs[l],αs[l+1],βs[l]) for l=1:(L-1)];
    julia> Ws = UniformScaling(1.0/maximum(maximum.(Ws))) * Ws
    julia> Z = CpelNano.get_Z(u1,uL,Ws);
    julia> CpelNano.get_marg_px(u1,uL,Ws,Z)
    10-element Array{Float64,1}:
     0.5
     0.5
     0.5
     0.5
     0.5
     0.5
     0.5
     0.5
     0.5
    ```
"""
function get_marg_px(u1::Vector{Float64},uL::Vector{Float64},Ws::Vector{Array{Float64,2}},Z::Float64)::Vector{Float64}

    # Compute P(X̄_l=1;ϕ) for l=1,2,…,L
    prod_Ws = prod(Ws)
    exs = fill(NaN,length(Ws)+1)
    @inbounds for l=1:length(exs)
        if l==1
            exs[1] = u1' * D_PX * prod_Ws * uL
        elseif l==length(exs)
            exs[end] = u1' * prod_Ws * D_PX * uL
        else
            exs[l] = u1' * prod(Ws[1:(l-1)]) * D_PX * prod(Ws[l:end]) * uL
        end
    end

    # Return vector of P(X̄_l=1;ϕ) for l=1,2,…,L
    return exs ./ Z

end
##################################################################################################
## Transfer matrix base functions
##################################################################################################
"""

    `log_sum_exp(x)`
    
    Computes log sum using vector x as log(1+sum(exp(x-x_max)))+xmax. 
    
    # Examples
    ```julia-repl
    julia> CpelNano.log_sum_exp([1.0e-2,1.0e-5,1.0e-10])
    0.6981596805576123
    ```
"""
function log_sum_exp(x::Vector{Float64})::Float64

    # Get max and min
    x_max = maximum(x)
    x_min = minimum(x)

    # Return
    return log(1.0+exp(x_min-x_max))+x_max

end
"""

    `log_mat_mult(logW1,logW2)`
    
    Computes log(W1*W2) in a numeric stable way.
    
    # Examples
    ```julia-repl
    julia> CpelNano.log_mat_mult(rand(2,2),rand(2,2))
    ```
"""
function log_mat_mult(logW1::Array{Float64,2},logW2::Array{Float64,2})::Array{Float64,2}

    # Init output
    logW1W2 = fill(NaN,(2,2))

    # Compute each entry
    logW1W2[1,1] = log_sum_exp(logW1[1,:]+logW2[:,1])
    logW1W2[1,2] = log_sum_exp(logW1[1,:]+logW2[:,2])
    logW1W2[2,1] = log_sum_exp(logW1[2,:]+logW2[:,1])
    logW1W2[2,2] = log_sum_exp(logW1[2,:]+logW2[:,2])

    # Return log(W1*W2)
    return logW1W2

end
"""

    `log_vec_mat_mult(logu,logW)`
    
    Computes log(u'*W) in a numeric stable way.
    
    # Examples
    ```julia-repl
    julia> CpelNano.log_vec_mat_mult(rand(2),rand(2,2))
    ```
"""
function log_vec_mat_mult(logu::Vector{Float64},logW::Array{Float64,2})::Vector{Float64}

    # Init output
    loguW = fill(NaN,2)

    # Compute each entry
    loguW[1] = log_sum_exp(logu+logW[:,1])
    loguW[2] = log_sum_exp(logu+logW[:,2])

    # Return log(u'*W)
    return loguW

end
"""

    `log_vec_vec_mult(logu1,logu2)`
    
    Computes log(u1'*u2) in a numeric stable way.
    
    # Examples
    ```julia-repl
    julia> CpelNano.log_vec_vec_mult(rand(2),rand(2))
    ```
"""
function log_vec_vec_mult(logu1::Vector{Float64},logu2::Vector{Float64})::Float64

    # Return log(u1'*u2)
    return log_sum_exp(logu1+logu2)

end
"""

    `mult_log_mats([log(M1),log(M2),...])`
    
    Computes log(M1*M2*...) in a numeric stable way.
    
    # Examples
    ```julia-repl
    julia> CpelNano.mult_log_mats(fill(rand(2,2),10))
    ```
"""
function mult_log_mats(logMs::Vector{Array{Float64,2}})::Array{Float64,2}

    # Check if only one matrix
    length(logMs)==1 && return logMs[1]
    
    # Perform first multiplication
    mult_logMs = log_mat_mult(logMs[1],logMs[2])
    length(logMs)==2 && return mult_logMs

    # Do the rest
    @inbounds for i=3:length(logMs)
        mult_logMs = log_mat_mult(mult_logMs,logMs[i])
    end

    # Return log of the multiplication
    return mult_logMs

end
##################################################################################################
## Transfer matrix method functions for conditional expectations E[⋅|y] using log trick
##################################################################################################
"""

    `get_γc(x,α,call)`
    
    Computes function γ in transfer matrix method for p(x|y).
    
    # Examples
    ```julia-repl
    julia> call = CpelNano.MethCallCpgGrp(-50.0,-45.0);
    julia> CpelNano.get_γc(false,0.0,call)
    -25.0
    julia> CpelNano.get_γc(true,0.0,call)
    -22.5
    ```
"""
function get_γc(x::Bool,α::Float64,call::MethCallCpgGrp)::Float64

    # Compute γ_{m,l} = 0.5⋅(α⋅x + 1(obs)⋅log p(y|x))
    γ = x ? α : - α
    if call.obs
        γ += x ? call.log_pyx_m : call.log_pyx_u
    end

    # Return
    return 0.5 * γ

end
"""

    `get_log_uc(α,call)`
    
    Generates vector log(u1) and log(uN) used in transfer matrix method for p(x|y).
    
    # Examples
    ```julia-repl
    julia> call = CpelNano.MethCallCpgGrp(-50.0,-45.0);
    julia> CpelNano.get_log_uc(0.0,call)
     -25.0
     -22.5
    ```
"""
function get_log_uc(α::Float64,call::MethCallCpgGrp)::Vector{Float64}
    
    # Return u for CpG group passed
    return [get_γc(false,α,call); get_γc(true,α,call)]

end
"""

    `get_log_Wc(α1,α2,β,call1,call2)`
    
    Generates transition matrix log(W) used in transfer matrix method for p(x|y).
    
    # Examples
    ```julia-repl
    julia> call = CpelNano.MethCallCpgGrp(-50.0,-45.0);
    julia> Wc = CpelNano.get_log_Wc(0.0,0.0,0.0,call,call)
    2×2 Array{Float64,2}:
     -50.0  -47.5
     -47.5  -45.0
    ```
"""
function get_log_Wc(α1::Float64,α2::Float64,β::Float64,call1::MethCallCpgGrp,call2::MethCallCpgGrp)::Array{Float64,2}
    
    # Define transition matrix
    log_W = zeros(2,2)
    log_W[1,1] = get_γc(false,α1,call1) + get_γc(false,α2,call2) + β
    log_W[2,1] = get_γc(true,α1,call1) + get_γc(false,α2,call2) - β
    log_W[1,2] = get_γc(false,α1,call1) + get_γc(true,α2,call2) - β
    log_W[2,2] = get_γc(true,α1,call1) + get_γc(true,α2,call2) + β

    # Return log_W
    return log_W

end
"""

    `get_log_Zc(log_u1,log_uN,log_Ws)`
    
    Computes log of partition function of P(X|Y;ϕ) using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-100.0,-100.0),L);
    julia> logu1 = CpelNano.get_log_uc(αs[1],obs[1]); 
    julia> loguL = CpelNano.get_log_uc(αs[end],obs[end]);
    julia> logWs = [CpelNano.get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> logZc = CpelNano.get_log_Zc(logu1,loguL,logWs)
    -993.06
    ```
"""
function get_log_Zc(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}})::Float64
   
    # Return log u_1⊺ W_1 ⋯ W_N u_N for m-th observation
    return log_vec_vec_mult(log_vec_mat_mult(logu1,mult_log_mats(logWs)),loguL)

end
"""

    `get_Ec_X_log(log_u1,log_uN,log_Ws,log_Zc)`
    
    Computes E[X̄_l|ỹ_m] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-45.0,-45.0),L);
    julia> logu1 = CpelNano.get_log_uc(αs[1],obs[1]); 
    julia> loguL = CpelNano.get_log_uc(αs[end],obs[end]);
    julia> logWs = [CpelNano.get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> logZc = CpelNano.get_log_Zc(logu1,loguL,logWs);
    julia> CpelNano.get_Ec_X_log(logu1,loguL,logWs,logZc)
    10-element Array{Float64,1}:
     3.552713678800501e-15
     3.552713678800501e-15
     3.552713678800501e-15
     3.552713678800501e-15
     3.552713678800501e-15
     3.552713678800501e-15
     3.552713678800501e-15
     3.552713678800501e-15
     3.552713678800501e-15
     3.552713678800501e-15
    ```
"""
function get_Ec_X_log(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}},logZc::Float64)::Vector{Float64}

    # Compute E[X̄_{m,l}|ỹ_m] for l=1,2,…,L
    log_exs = Vector{Float64}()
    @inbounds for l=1:(length(logWs)+1)
        if l==1
            log_DWs = mult_log_mats(vcat([logD1],logWs))
            log_u1DWs = log_vec_mat_mult(logu1,log_DWs)
            push!(log_exs,log_vec_vec_mult(log_u1DWs,loguL))
        elseif l==(length(logWs)+1)
            log_WsD = mult_log_mats(vcat(logWs,[logD1]))
            log_u1WsD = log_vec_mat_mult(logu1,log_WsD)
            push!(log_exs,log_vec_vec_mult(log_u1WsD,loguL))
        else
            log_WsDWs = mult_log_mats(vcat(logWs[1:(l-1)],[logD1],logWs[l:end]))
            log_u1WsDWs = log_vec_mat_mult(logu1,log_WsDWs)
            push!(log_exs,log_vec_vec_mult(log_u1WsDWs,loguL))
        end
    end

    # Return E[X̄_m|ỹ_m]
    return exp.(log_exs .- logZc) .- 2.0

end
"""

    `get_Ec_logpyx_log(log_u1,log_uN,log_Ws,log_Zc)`
    
    Computes E[log p(y|x)|ỹ_m] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(-1.0,L); βs = fill(1.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-105.0,-95.0),L);
    julia> logu1 = CpelNano.get_log_uc(αs[1],obs[1]); 
    julia> loguL = CpelNano.get_log_uc(αs[end],obs[end]);
    julia> logWs = [CpelNano.get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> logZc = CpelNano.get_log_Zc(logu1,loguL,logWs);
    julia> CpelNano.get_Ec_logpyx_log(logu1,loguL,logWs,logZc,obs)
    10-element Array{Float64,1}:
     -95.00045412823272
     -95.00006161154741
     -95.00006148227914
     -95.00006148224674
     -95.00006148224674
     -95.00006148223594
     -95.00006148223594
     -95.00006148226834
     -95.00006161153661
     -95.00045412823272
    ```
"""
function get_Ec_logpyx_log(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}},logZc::Float64,call_m::Vector{MethCallCpgGrp})::Vector{Float64}

    # Compute E[log p(y|x)|ỹ]
    log_pyxs = zeros(length(call_m))
    @inbounds for l=1:(length(logWs)+1)
        call_m[l].obs || continue
        logD = log.(max.([call_m[l].log_pyx_u+250.0 0.0; 0.0 call_m[l].log_pyx_m+250.0],1.0))
        if l==1
            log_DWs = mult_log_mats(vcat([logD],logWs))
            log_u1DWs = log_vec_mat_mult(logu1,log_DWs)
            log_pyxs[l] = log_vec_vec_mult(log_u1DWs,loguL)
        elseif l==(length(logWs)+1)
            log_WsD = mult_log_mats(vcat(logWs,[logD]))
            log_u1WsD = log_vec_mat_mult(logu1,log_WsD)
            log_pyxs[l] = log_vec_vec_mult(log_u1WsD,loguL)
        else
            log_WsDWs = mult_log_mats(vcat(logWs[1:(l-1)],[logD],logWs[l:end]))
            log_u1WsDWs = log_vec_mat_mult(logu1,log_WsDWs)
            log_pyxs[l] = log_vec_vec_mult(log_u1WsDWs,loguL)
        end
        log_pyxs[l] = exp(log_pyxs[l] - logZc) - 250.0
    end

    # Return E[log p(y|x)|ỹ]
    return log_pyxs

end
"""

    `get_Ec_XX_log(log_u1,log_uN,log_Ws,log_Zc)`
    
    Computes vector of E[X̄_{m,l}X̄_{m,l+1}|ỹ_m] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-100.0,-100.0),L);
    julia> logu1 = CpelNano.get_log_uc(αs[1],obs[1]);
    julia> loguL = CpelNano.get_log_uc(αs[end],obs[end]);
    julia> logWs = [CpelNano.get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> logZc = CpelNano.get_log_Zc(logu1,loguL,logWs);
    julia> CpelNano.get_Ec_XX_log(logu1,loguL,logWs,logZc)
    9-element Array{Float64,1}:
     -1.099120794378905e-13 
     -1.099120794378905e-13 
     -1.099120794378905e-13 
     -1.099120794378905e-13 
     -1.099120794378905e-13 
     1.1723955140041653e-13
     1.1723955140041653e-13
     1.1723955140041653e-13
     1.1723955140041653e-13
    ```
"""
function get_Ec_XX_log(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}},logZc::Float64)::Vector{Float64}

    # Compute E[X̄_{m,l}X̄_{m,l+1}|ỹ_m] for l=1,2,…,L-1
    log_exxs = Vector{Float64}()
    @inbounds for l=1:length(logWs)
        if l==1
            aux = logWs[1]+logD2
            log_DWs = mult_log_mats(vcat([aux],logWs[2:end]))
            log_u1DWs = log_vec_mat_mult(logu1,log_DWs)
            push!(log_exxs,log_vec_vec_mult(log_u1DWs,loguL))
        elseif l==length(logWs)
            aux = logWs[end]+logD2
            log_WsD = mult_log_mats(vcat(logWs[1:(end-1)],[aux]))
            log_u1WsD = log_vec_mat_mult(logu1,log_WsD)
            push!(log_exxs,log_vec_vec_mult(log_u1WsD,loguL))
        else
            aux = logWs[l]+logD2
            log_WsDWs = mult_log_mats(vcat(logWs[1:(l-1)],[aux],logWs[(l+1):end]))
            log_u1WsDWs = log_vec_mat_mult(logu1,log_WsDWs)
            push!(log_exxs,log_vec_vec_mult(log_u1WsDWs,loguL))
        end
        
    end

    # Return vector of E[X̄_{m,l}X̄_{m,l+1}|ỹ_m]
    return exp.(log_exxs .- logZc) .- 2.0

end
"""

    `get_cond_exs_log(αs,βs,obs)`
    
    Computes conditional expectations using the transfer matrix method. This function uses the
    log-sum trick to cope with the numerical instability.
    
    # Examples
    ```julia-repl
    julia> L = 2; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-50.0,-45.0),L);
    julia> ex,exx = CpelNano.get_cond_exs_log(αs,βs,obs)
    ([0.9866142981513866, 0.9866142981513866], [0.9734077733168003])
    julia> obs = fill(CpelNano.MethCallCpgGrp(-200.0,-200.0),L);
    julia> ex,exx = CpelNano.get_cond_exs_log(αs,βs,obs)
    ([3.552713678800501e-15, 3.552713678800501e-15], [3.552713678800501e-15])
    ```
"""
function get_cond_exs_log(αs::Vector{Float64},βs::Vector{Float64},obs::Vector{MethCallCpgGrp})::NTuple{2,Vector{Float64}}

    # Get log u's
    logu1 = get_log_uc(αs[1],obs[1])
    loguL = get_log_uc(αs[end],obs[end])

    # Get log W's
    logWs = [get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:length(βs)]
    
    # Get Zc
    logZc = get_log_Zc(logu1,loguL,logWs)

    # Ex
    ex = get_Ec_X_log(logu1,loguL,logWs,logZc)

    # Exx
    exx = get_Ec_XX_log(logu1,loguL,logWs,logZc)

    # Return E[X̄] and E[X̄X̄]
    return ex,exx
    
end
##################################################################################################
## Transfer matrix method functions for WGBS inference using log sum trick (not used)
##################################################################################################
"""

    `get_Ec_X_wgbs_log(u1,uL,Ws,Zc,obs)`
    
    Computes E[X̄|y_m;ϕ] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 2; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-2.0,0.0),L);
    julia> logu1 = CpelNano.get_log_uc(αs[1],obs[1]); loguL = CpelNano.get_log_uc(αs[end],obs[end]);
    julia> logWs = [CpelNano.get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> logZc = CpelNano.get_log_Zc(logu1,loguL,logWs);
    julia> CpelNano.get_Ec_X_wgbs_log(logu1,loguL,logWs,logZc,obs)
    2-element Array{Float64,1}:
     1.0
     1.0
    ```
"""
function get_Ec_X_wgbs_log(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}},logZc::Float64,obs::Vector{MethCallCpgGrp})::Vector{Float64}
   
    # Compute E[X̄_l|y_m] for l=1,2,…,L
    log_exs = fill(NaN,length(obs))
    @inbounds for l=1:length(obs)

        # Check if observed
        if obs[l].obs
            
            # If observed just take observed value & continue
            log_exs[l] = obs[l].log_pyx_m > obs[l].log_pyx_u ? log(3.0) : log(1.0)
            
        else
        
            # Compute E[X̄_{m,l}|ỹ_m] for l=1,2,…,L
            if l==1
                log_DWs = mult_log_mats(vcat([logD1],logWs))
                log_u1DWs = log_vec_mat_mult(logu1,log_DWs)
                log_exs[l] = log_vec_vec_mult(log_u1DWs,loguL) - logZc
            elseif l==(length(logWs)+1)
                log_WsD = mult_log_mats(vcat(logWs,[logD1]))
                log_u1WsD = log_vec_mat_mult(logu1,log_WsD)
                log_exs[l] = log_vec_vec_mult(log_u1WsD,loguL) - logZc
            else
                log_WsDWs = mult_log_mats(vcat(logWs[1:(l-1)],[logD1],logWs[l:end]))
                log_u1WsDWs = log_vec_mat_mult(logu1,log_WsDWs)
                log_exs[l] = log_vec_vec_mult(log_u1WsDWs,loguL) - logZc
            end
        
        end

    end

    # Return vector of Ec[X̄_l|y_m] for m-th observation
    return exp.(log_exs) .- 2.0
    # return exs

end
"""

    `get_Ec_XX_wgbs_log(u1,uL,Ws,obs)`
    
    Computes E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 3; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-5.0,0.0),L);
    julia> logu1 = CpelNano.get_log_uc(αs[1],obs[1]); loguL = CpelNano.get_log_uc(αs[end],obs[end]);
    julia> logWs = [CpelNano.get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:(L-1)];
    julia> logZc = CpelNano.get_log_Zc(logu1,loguL,logWs);
    julia> CpelNano.get_Ec_XX_wgbs_log(logu1,loguL,logWs,logZc,obs)
    2-element Array{Float64,1}:
     1.0
     1.0
    ```
"""
function get_Ec_XX_wgbs_log(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}},logZc::Float64,obs::Vector{MethCallCpgGrp})::Vector{Float64}
    
    # Init output
    log_exxs = fill(NaN,length(logWs))
    
    # Compute E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] for l=1,2,…,L-1
    @inbounds for l=1:length(logWs)

        if obs[l].obs && obs[l+1].obs
        
            # If both are observed just take observed value & continue
            xl = obs[l].log_pyx_m > obs[l].log_pyx_u ? 1.0 : -1.0
            xlp1 = obs[l+1].log_pyx_m > obs[l+1].log_pyx_u ? 1.0 : -1.0
            log_exxs[l] = xl==xlp1 ? log(3.0) : log(1.0)

        else

            # If not observed compute E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ]
            if l==1
                aux = logWs[1]+logD2
                log_DWs = mult_log_mats(vcat([aux],logWs[2:end]))
                log_u1DWs = log_vec_mat_mult(logu1,log_DWs)
                log_exxs[l] = log_vec_vec_mult(log_u1DWs,loguL) - logZc
            elseif l==length(logWs)
                aux = logWs[end]+logD2
                log_WsD = mult_log_mats(vcat(logWs[1:(end-1)],[aux]))
                log_u1WsD = log_vec_mat_mult(logu1,log_WsD)
                log_exxs[l] = log_vec_vec_mult(log_u1WsD,loguL) - logZc
            else
                aux = logWs[l]+logD2
                log_WsDWs = mult_log_mats(vcat(logWs[1:(l-1)],[aux],logWs[(l+1):end]))
                log_u1WsDWs = log_vec_mat_mult(logu1,log_WsDWs)
                log_exxs[l] = log_vec_vec_mult(log_u1WsDWs,loguL) - logZc
            end
            
        end

    end

    # Return vector of E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] for m-th observation
    return exp.(log_exxs) .- 2.0

end
"""

    `get_cond_exs_wgbs_log(αs,βs,obs)`
    
    Computes E[X̄_{m,l}|Y_m=y_m;ϕ] & E[X̄_{m,l}X̄_{m,l+1}|Y_m=y_m;ϕ] using the transfer matrix method. This function 
    does not use the log-sum trick, however. Thus, it is not as numerically stable as `get_cond_exs_log()`.
    
    # Examples
    ```julia-repl
    julia> L = 5; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-2.0,0.0),L);
    julia> ex,exx = CpelNano.get_cond_exs_wgbs_log(αs,βs,obs);
    julia> obs = fill(CpelNano.MethCallCpgGrp(0.0,-2.0),L);
    julia> ex,exx = CpelNano.get_cond_exs_wgbs_log(αs,βs,obs);
    ```
"""
function get_cond_exs_wgbs_log(αs::Vector{Float64},βs::Vector{Float64},obs::Vector{MethCallCpgGrp})::NTuple{2,Vector{Float64}}

    ## Get matrices

    # Get u's
    logu1 = get_log_uc(αs[1],obs[1])
    loguL = get_log_uc(αs[end],obs[end])

    # Get W's
    logWs = [get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:length(βs)]

    # Get Zc
    logZc = get_log_Zc(logu1,loguL,logWs)

    ## Compute expectations

    # Ex
    ex = get_Ec_X_wgbs_log(logu1,loguL,logWs,logZc,obs)

    # Exx
    exx = get_Ec_XX_wgbs_log(logu1,loguL,logWs,logZc,obs)

    # Return
    return ex,exx
    
end
##################################################################################################
## Transfer matrix method functions for expectations E[X̄] && E[X̄X̄] using log-sum trick
##################################################################################################
"""

    `get_log_u(α)`
    
    Generates vector log(u) used in transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_log_u(0.0)
    2-element Array{Float64,1}:
     0.0
     0.0
    ```
"""
function get_log_u(α::Float64)::Vector{Float64}
    
    # Compute intermediate quantity
    aux = 0.5 * α
    
    # Return u for CpG passed
    return [-aux;aux]

end
"""

    `get_log_W(α1,α2,β)`
    
    Generates transition matrix log(W) used in transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_log_W(0.0,0.0,0.0)
    2×2 Array{Float64,2}:
     0.0  0.0
     0.0  0.0
    ```
"""
function get_log_W(α1::Float64,α2::Float64,β::Float64)::Array{Float64,2}
    
    # Define transition matrix
    logW = zeros(2,2)
    logW[1,1] = -0.5*(α1+α2) + β
    logW[2,1] = 0.5*(α1-α2) - β
    logW[1,2] = 0.5*(α2-α1) - β
    logW[2,2] = 0.5*(α1+α2) + β

    # Return logW
    return logW

end
"""

    `get_log_Z(u1,uN,Ws)`
    
    Computes log of partition function of q(x̄;ϕ) using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> logu1 = CpelNano.get_log_u(αs[1]); 
    julia> loguL = CpelNano.get_log_u(αs[end]);
    julia> logWs = [CpelNano.get_log_W(αs[l],αs[l+1],βs[l]) for l=1:length(βs)];
    julia> logZ = CpelNano.get_log_Z(logu1,loguL,logWs)
    6.931471805599453
    ```
"""
function get_log_Z(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}})::Float64
   
    # Compute intermediate quantities
    aux = mult_log_mats(logWs)
    aux = log_vec_mat_mult(logu1,aux)
    
    # Return log u_1⊺ W_1 ⋯ W_L u_L
    return log_vec_vec_mult(aux,loguL)

end
"""
    `get_∇logZ(REG_ST,ϕhat)`

    Numerically computes the gradient of the log partition function of a model with estimated parameter 
    vector ϕhat=[a,b,c]. This is equivalent to computing the expected value of sufficient statistis (SS).

    # Examples
    ```julia-repl
    julia> x=CpelNano.RegStruct(); x.L=10; x.Nl=fill(5.0,x.L); x.ρl=fill(0.1,x.L); x.dl=fill(10.0,x.L);
    julia> CpelNano.get_∇logZ(x,[0.0,0.0,0.0])
    3-element Array{Float64,1}:
     0.0
     0.0
     0.0 
    ```
"""
function get_∇logZ(rg::RegStruct,ϕhat::Vector{Float64})::Vector{Float64}

    # Define function
    function f(ϕ::Vector{Float64})

        # Get αs & βs
        αs,βs = get_αβ_from_ϕ(ϕ,rg)

        # Get matrices
        logu1 = get_log_u(αs[1])
        loguL = get_log_u(αs[end])
        logWs = [get_log_W(αs[l],αs[l+1],βs[l]) for l=1:length(βs)]
        
        # Return log ζ(ϕ)
        return get_log_Z(logu1,loguL,logWs)

    end

    # Return ∇logζ(ϕ)
    return Calculus.gradient(f,ϕhat)

end
"""
    `get_∇logZ!(REG_ST)`

    Numerically computes the gradient of the log partition function of a model and stores it in the passed
    RegStruct REG_DATA.

    # Examples
    ```julia-repl
    julia> x=CpelNano.RegStruct(); x.L=10; x.Nl=fill(5.0,x.L); x.ρl=fill(0.1,x.L); x.dl=fill(10.0,x.L);
    julia> x.ϕhat=[0.0,0.0,0.0]; CpelNano.get_∇logZ!(x); x.∇logZ
    3-element Array{Float64,1}:
     0.0
     0.0
     0.0
    ```
"""
function get_∇logZ!(rg::RegStruct)::Nothing

    # Set ∇logZ
    rg.∇logZ  = get_∇logZ(rg,rg.ϕhat)

    # Return nothing
    return nothing

end
"""

    `get_E_X_log(log_u1,log_uL,log_Ws,log_Z)`
    
    Computes E[X̄_l;ϕ] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> logu1 = CpelNano.get_log_u(αs[1]); 
    julia> loguL = CpelNano.get_log_u(αs[end]);
    julia> logWs = [CpelNano.get_log_W(αs[l],αs[l+1],βs[l]) for l=1:(L-1)];
    julia> logZ = CpelNano.get_log_Z(logu1,loguL,logWs);
    julia> CpelNano.get_E_X_log(logu1,loguL,logWs,logZ)
    10-element Array{Float64,1}:
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
function get_E_X_log(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}},logZ::Float64)::Vector{Float64}

    # Compute E[X̄_l;ϕ] for l=1,2,…,L
    log_exs = Vector{Float64}()
    @inbounds for l=1:(length(logWs)+1)
        if l==1
            log_DWs = mult_log_mats(vcat([logD1],logWs))
            log_u1DWs = log_vec_mat_mult(logu1,log_DWs)
            push!(log_exs,log_vec_vec_mult(log_u1DWs,loguL))
        elseif l==(length(logWs)+1)
            log_WsD = mult_log_mats(vcat(logWs,[logD1]))
            log_u1WsD = log_vec_mat_mult(logu1,log_WsD)
            push!(log_exs,log_vec_vec_mult(log_u1WsD,loguL))
        else
            log_WsDWs = mult_log_mats(vcat(logWs[1:(l-1)],[logD1],logWs[l:end]))
            log_u1WsDWs = log_vec_mat_mult(logu1,log_WsDWs)
            push!(log_exs,log_vec_vec_mult(log_u1WsDWs,loguL))
        end
    end

    # Return vector of E[X̄_l;ϕ] for l=1,2,…,L
    return exp.(log_exs .- logZ) .- 2.0

end
"""

    `get_E_XX_log(log_u1,log_uL,log_Ws,log_Z)`
    
    Computes E[X̄_{l}X̄_{l+1};ϕ] using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> logu1 = CpelNano.get_log_u(αs[1]);
    julia> loguL = CpelNano.get_log_u(αs[end]);
    julia> logWs = [CpelNano.get_log_W(αs[l],αs[l+1],βs[l]) for l=1:(L-1)];
    julia> logZ = CpelNano.get_log_Z(logu1,loguL,logWs);
    julia> CpelNano.get_E_XX_log(logu1,loguL,logWs,logZ)
    9-element Array{Float64,1}:
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
function get_E_XX_log(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}},logZ::Float64)::Vector{Float64}

    # Compute E[X̄_{l}X̄_{l+1};ϕ] for l=1,2,…,L-1
    log_exxs = Vector{Float64}()
    @inbounds for l=1:length(logWs)
        if l==1
            aux = logWs[1]+logD2
            log_DWs = mult_log_mats(vcat([aux],logWs[2:end]))
            log_u1DWs = log_vec_mat_mult(logu1,log_DWs)
            push!(log_exxs,log_vec_vec_mult(log_u1DWs,loguL))
        elseif l==length(logWs)
            aux = logWs[end]+logD2
            log_WsD = mult_log_mats(vcat(logWs[1:(end-1)],[aux]))
            log_u1WsD = log_vec_mat_mult(logu1,log_WsD)
            push!(log_exxs,log_vec_vec_mult(log_u1WsD,loguL))
        else
            aux = logWs[l]+logD2
            log_WsDWs = mult_log_mats(vcat(logWs[1:(l-1)],[aux],logWs[(l+1):end]))
            log_u1WsDWs = log_vec_mat_mult(logu1,log_WsDWs)
            push!(log_exxs,log_vec_vec_mult(log_u1WsDWs,loguL))
        end
        
    end

    # Return vector of E[X̄_{l}X̄_{l+1};ϕ] for l=1,2,…,L-1
    return exp.(log_exxs .- logZ) .- 2.0

end
"""

    `get_marg_px_log(log_u1,log_uL,log_Ws,log_Z)`
    
    Computes marginal probability p(x_l) using the transfer matrix method.
    
    # Examples
    ```julia-repl
    julia> L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> logu1 = CpelNano.get_log_u(αs[1]); 
    julia> loguL = CpelNano.get_log_u(αs[end]);
    julia> logWs = [CpelNano.get_log_W(αs[l],αs[l+1],βs[l]) for l=1:(L-1)];
    julia> logZ = CpelNano.get_log_Z(logu1,loguL,logWs);
    julia> CpelNano.get_marg_px_log(logu1,loguL,logWs,logZ)
    10-element Array{Float64,1}:
     0.5
     0.5
     0.5
     0.5
     0.5
     0.5
     0.5
     0.5
     0.5
    ```
"""
function get_marg_px_log(logu1::Vector{Float64},loguL::Vector{Float64},logWs::Vector{Array{Float64,2}},logZ::Float64)::Vector{Float64}

    # Compute p(xl) for l=1,2,…,L
    log_exs = Vector{Float64}()
    @inbounds for l=1:(length(logWs)+1)
        if l==1
            log_DWs = mult_log_mats(vcat([logD3],logWs))
            log_u1DWs = log_vec_mat_mult(logu1,log_DWs)
            push!(log_exs,log_vec_vec_mult(log_u1DWs,loguL))
        elseif l==(length(logWs)+1)
            log_WsD = mult_log_mats(vcat(logWs,[logD3]))
            log_u1WsD = log_vec_mat_mult(logu1,log_WsD)
            push!(log_exs,log_vec_vec_mult(log_u1WsD,loguL))
        else
            log_WsDWs = mult_log_mats(vcat(logWs[1:(l-1)],[logD3],logWs[l:end]))
            log_u1WsDWs = log_vec_mat_mult(logu1,log_WsDWs)
            push!(log_exs,log_vec_vec_mult(log_u1WsDWs,loguL))
        end
    end

    # Return vector of p(xl) for l=1,2,…,L
    return exp.(log_exs .- logZ) .- 2.0

end