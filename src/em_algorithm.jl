##################################################################################################
## Common functions
##################################################################################################
"""

    `conv_check(ϕ1,ϕ2)`
    
    Check for convergence in the parameter vector estimate. 
    
    # Examples
    ```julia-repl
    julia> CpelNano.conv_check(zeros(10),0.1.+zeros(10))
    false
    julia> CpelNano.conv_check(zeros(10),0.01.+zeros(10))
    true
    ```
"""
function conv_check(ϕ1::Vector{Float64},ϕ2::Vector{Float64})::Bool
    
    # Return
    return (sqrt(sum((ϕ1-ϕ2).^2))<=3.0e-2)

end
"""

    `pick_ϕ(pairs)`
    
    Pick ϕhat that maximizes objective function Q.
    
    # Examples
    ```julia-repl
    julia> CpelNano.pick_ϕ(pairs)
    ```
"""
function pick_ϕ(pairs::Vector{Tuple{Vector{Float64},Float64}})::Tuple{Vector{Float64},Float64}
    
    # Init
    Qhat = -Inf
    ϕhat = Vector{Float64}()

    # Choose ϕ that maximizes E[p(X,y)|Y=y;ϕ]
    @inbounds for pair in pairs
        # Get ϕ and Q
        ϕ = pair[1]
        Q = pair[2]
        # Update pair if improvement
        ϕhat,Qhat = Q>Qhat ? (ϕ,Q) : (ϕhat,Qhat)
    end

    # Return optimal pair
    return ϕhat,Qhat

end
"""

    `gen_ϕ0s(max_init)`
    
    Function that returns a set of ϕ0 to initialize with the EM algorithm.
    
    # Examples
    ```julia-repl
    julia> CpelNano.gen_ϕ0s(10)
    ```
"""
function gen_ϕ0s(max_init::Int64)::Vector{Vector{Float64}}
    
    # Init
    ϕ0s = Vector{Vector{Float64}}()

    # Generate 
    @inbounds for _ in 1:max_init
        push!(ϕ0s,vcat(rand(-5.0:0.01:5.0),rand(-50.0:0.01:50.0),rand(0.0:0.01:10.0)))
    end

    # Return parameter estimates
    return ϕ0s

end
"""

    `get_ess(exs)`
    
    Computes E[Sa(X̄)|y_1,…,y_m;ϕ], E[Sb(X̄)|y_1,…,y_m;ϕ], and E[Sc(X̄)|y_1,…,y_m;ϕ] from first and second 
    order stats E[X̄|y_m] and E[X̄X̄|y_m].
    
    # Examples
    ```julia-repl
    julia> rs = CpelNano.RegStruct(); rs.L = 10; rs.Nl = fill(5.0,rs.L); 
    julia> rs.ρl = fill(0.1,rs.L); rs.dl = fill(10.0,rs.L-1); 
    julia> αs = fill(0.0,rs.L); βs = fill(0.0,rs.L-1);
    julia> obs = fill(CpelNano.MethCallCpgGrp(-50.0,-50.0),rs.L);
    julia> exs1 = CpelNano.get_cond_exs_log(αs,βs,obs);
    julia> exs2 = CpelNano.get_cond_exs_log(αs,βs,obs);
    julia> exs = [exs1,exs2];
    julia> ess = CpelNano.get_ess([exs1,exs2],rs)
    (3.552713678800501e-13, 3.552713678800501e-13, 6.394884621840901e-15)
    ```
"""
function get_ess(exs::Vector{NTuple{2,Vector{Float64}}},rs::RegStruct)::Vector{Float64}

    # Compute E[S(X̄)|y_1,…,y_m]
    ES = fill(0.0,3)
    @inbounds for m=1:length(exs)
        EX,EXX = exs[m]
        ES[1] += dot(rs.Nl,EX)
        ES[2] += dot(rs.Nl.*rs.ρl,EX)
        ES[3] += dot(1.0./rs.dl,EXX)
    end

    # Return E[Sa(X̄)|y_1,…,y_m;ϕ], E[Sb(X̄)|y_1,…,y_m;ϕ], and E[Sc(X̄)|y_1,…,y_m;ϕ]
    return ES

end
##################################################################################################
## EM algorithm for nanopore sequencing
##################################################################################################
"""

    `em_iter(αs,βs,data,config)`
    
    A single iteration of the EM algorithm.
    
    # Examples
    ```julia-repl
    julia> M = 10; L = 10; αs = fill(0.0,L); βs = fill(0.05,L-1);
    julia> pobs = 1.0; μ_m = 40.0; σ_m = 2.0; μ_u = 80.0; σ_u = 2.0;
    julia> rs = CpelNano.cpel_samp_ont(M,αs,βs,pobs,μ_m,σ_m,μ_u,σ_u);
    julia> rs.L = 10; rs.Nl = fill(5.0,rs.L); rs.ρl = fill(0.1,rs.L); rs.dl = fill(10.0,rs.L-1);
    julia> CpelNano.em_iter(zeros(3),rs)
    ```
"""
function em_iter(ϕinit::Vector{Float64},rs::RegStruct)::Vector{Float64}

    # Get α and β
    αs,βs = CpelNano.get_αβ_from_ϕ(ϕinit,rs)
        # print_log("αs: $(αs)")
        # print_log("βs: $(βs)")

    ## E-step 

    # Compute E[X̄|y_m] & E[X̄X̄|y_m] for m=1,2,…,M
    map_out = map(x->CpelNano.get_cond_exs_log(αs,βs,x),rs.calls)
        # print_log("map_out: $(map_out)")

    # Compute E[Sa(X̄)|y_1,…,y_m;ϕ], E[Sb(X̄)|y_1,…,y_m;ϕ], and E[Sc(X̄)|y_1,…,y_m;ϕ]
    ES = CpelNano.get_ess(map_out,rs)
        # print_log("ES: $(ES)")

    ## M-step

    # Set up system (requires to work w/ vectors)
    function f!(U::Vector{Float64},ϕ::Vector{Float64})
        ∇logZ = CpelNano.get_∇logZ(rs,ϕ)
        U[1] =  ∇logZ[1] - ES[1] / length(rs.calls)
        U[2] =  ∇logZ[2] - ES[2] / length(rs.calls)
        U[3] =  ∇logZ[3] - ES[3] / length(rs.calls)
    end
    
    # Solve NL system
    sol =
    try
        # Call solver
        NLsolve.nlsolve(f!,ϕinit;iterations=20,ftol=1e-3)
    catch x
        # Report x error if found
        print_log("An instance of EM algorithm did not converge.")
    end

    ## Convergence
    
    # Leave if no convergence with previous estimate
    (isnothing(sol) || !converged(sol)) && return ϕinit

    # Project vector to allowed hypercube
    θproj = sol.zero
    θproj[1] = max(-AMAX,min(θproj[1],AMAX))
    θproj[2] = max(-BMAX,min(θproj[2],BMAX))
    θproj[3] = max(-CMAX,min(θproj[3],CMAX))

    # Return new ϕ
    return θproj
    
end
"""
    `comp_Q(REG_ST,ϕ)`

    Function that computes E[log p(X,y)|Y=y;ϕ].

    # Examples
    ```julia-repl
    julia> M = 50; L = 80; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> pobs = 1.0; μ_m = 40.0; σ_m = 2.0; μ_u = 80.0; σ_u = 2.0;
    julia> rs = CpelNano.cpel_samp_ont(M,αs,βs,pobs,μ_m,σ_m,μ_u,σ_u);
    julia> rs.L = L; rs.Nl = fill(5.0,rs.L); rs.ρl = fill(0.1,rs.L); rs.dl = fill(10.0,rs.L-1); 
    julia> CpelNano.comp_Q(rs,zeros(3))
    -2.8009282290116913
    ```
"""
function comp_Q(rs::RegStruct,ϕ::Vector{Float64})::Float64
    
    # Get αs & βs
    αs,βs = get_αβ_from_ϕ(ϕ,rs)

    # Compute Ec[X|y_m] & Ec[XX|y_m] for m=1,2,…,M
    map_out = map(x->get_cond_exs_log(αs,βs,x),rs.calls)

    # Compute Ec[S(X)] & Ec[SS(X)]
    ES = get_ess(map_out,rs)

    # Add <ϕ,E[S]>
    Q = rs.m * dot(ϕ,ES)

    # Matrices
    u1 = CpelNano.get_log_u(αs[1])
    uL = CpelNano.get_log_u(αs[end])
    Ws = [CpelNano.get_log_W(αs[l],αs[l+1],βs[l]) for l=1:length(βs)]

    # Add - m⋅log[Z(ϕ)]
    Q -= rs.m * CpelNano.get_log_Z(u1,uL,Ws)
    # print_log("Q before: $(Q)")

    # Add observational part
    @inbounds for m=1:rs.m

        # Get m-th observation
        obs = rs.calls[m]
    
        # Get u's
        u1 = get_log_uc(αs[1],obs[1])
        uL = get_log_uc(αs[end],obs[end])

        # Get W's
        Ws = [get_log_Wc(αs[l],αs[l+1],βs[l],obs[l],obs[l+1]) for l=1:length(βs)]

        # Get Zc
        Zc = get_log_Zc(u1,uL,Ws)
        # print_log("Zc: $(Zc)")

        # Scale Ws for numerical stability
        Q += sum(get_Ec_logpyx_log(u1,uL,Ws,Zc,obs))
        
    end
    # print_log("Q after: $(Q)")
    # readline()

    # Return value
    return log(Q)/rs.m

end
"""

    `em_inst(data,ϕ0,config)`
    
    Single instance of the EM algorithm.
    
    # Examples
    ```julia-repl
    julia> M = 10; L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> pobs = 1.0; μ_m = 2.0; σ_m = 2.0; μ_u = -2.0; σ_u = 2.0;
    julia> reg_dat = CpelNano.cpel_samp_ont(M,αs,βs,pobs,μ_m,σ_m,μ_u,σ_u);
    julia> ϕhat,Q = CpelNano.em_inst(reg_dat,Float64[],CpelNano.CpelNanoConfig())
    ([0.5210509995207954, 0.6242225806898237, 0.05607186596907174, 0.1085379770548646], -0.26316332446017654)
    ```
"""
function em_inst(rs::RegStruct,ϕ0::Vector{Float64},config::CpelNanoConfig)::Tuple{Vector{Float64},Float64}

    # Initialize
    i = 1
    ϕt = ϕ0
    conv = false
    max_iters_hit = false

    config.verbose && print_log("Starting EM algorithm instance...")
    config.verbose && print_log("ϕt=$(ϕt); Q(ϕ)=$(comp_Q(rs,ϕt))")

    # Iterate
    while !conv && !max_iters_hit

        # Find next ϕ update
        ϕtp1 = em_iter(ϕt,rs)
        
        # Print out Q
        config.verbose && print_log("ϕtp1=$(ϕtp1); Q(ϕ)=$(comp_Q(rs,ϕtp1))")

        # Check for convergence
        conv = conv_check(ϕt,ϕtp1)
        
        # Check if too many iters
        max_iters_hit = i>=config.max_em_iters

        # Update ϕ
        ϕt = ϕtp1

        # Increase counter
        i += 1

    end
    
    # Check if more than 2 iterations (sometimes it gets stuck right at start)
    conv = conv && i>2

    # Get final pair
    Q = conv ? comp_Q(rs,ϕt) : -Inf

    # Return pair (ϕ,Q)
    return ϕt,Q

end
"""

    `get_ϕhat!(REG_ST,CONFIG)`
    
    Run multiple instances of the EM algorithm by randomly initializing the algorithm max_init times.
    
    # Examples
    ```julia-repl
    julia> a = -2.5; b = 20.0; c = 5.0; ϕ = [a,b,c]; L = 50;
    julia> x=CpelNano.RegStruct(); x.L=L; x.Nl=fill(1.0,x.L);
    julia> x.ρl=rand(0.001:0.001:0.1,x.L); x.dl=fill(10.0,x.L-1);
    julia> α,β = CpelNano.get_αβ_from_ϕ(ϕ,x);
    julia> M = 20; pobs = 1.0; μ_m = 40.0; σ_m = 2.0; μ_u = 45.0; σ_u = 2.0;
    julia> rs = CpelNano.cpel_samp_ont(M,α,β,pobs,μ_m,σ_m,μ_u,σ_u);
    julia> rs.L=x.L; rs.Nl=x.Nl; rs.ρl=x.ρl; rs.dl=x.dl; 
    julia> config = CpelNano.CpelNanoConfig(); config.verbose=true;
    julia> config.max_em_iters = 50; config.max_em_init = 5;
    julia> CpelNano.get_ϕhat!(rs,config); rs.ϕhat
    ```
"""
function get_ϕhat!(rs::RegStruct,config::CpelNanoConfig)::Nothing
    
    # Initialize vector
    ϕ0s = gen_ϕ0s(config.max_em_init)

    # Launch multiple instances of EM in parallel
    config.verbose && print_log("Starting EM algorithm...")
    em_out = map(ϕ0->em_inst(rs,ϕ0,config),ϕ0s)

    # Choose ϕ with maximum expected loglike
    ϕhat,Qhat = pick_ϕ(em_out)
    config.verbose && print_log("Converged to Q=$(Qhat) and ϕ̂=$(ϕhat).")

    # Set ϕhat
    rs.ϕhat = Qhat>-Inf ? ϕhat : Vector{Float64}()

    # Return nothing
    return nothing

end
##################################################################################################
## EM algorithm for WGBS sequencing
##################################################################################################
"""

    `em_iter_wgbs(αs,βs,data,config)`
    
    A single iteration of the EM algorithm in WGBS case.
    
    # Examples
    ```julia-repl
    julia> M = 10; L = 10; αs = fill(0.0,L); βs = fill(0.05,L-1); pobs = 1.0;
    julia> rs = CpelNano.cpel_samp_wgbs(M,αs,βs,pobs);
    julia> rs.L = 10; rs.m = M; rs.Nl = fill(5.0,rs.L); rs.ρl = fill(0.1,rs.L); rs.dl = fill(10.0,rs.L-1);
    julia> CpelNano.em_iter_wgbs(zeros(3),rs)
    ```
"""
function em_iter_wgbs(ϕinit::Vector{Float64},rs::RegStruct)::Vector{Float64}

    # Get α and β
    αs,βs = CpelNano.get_αβ_from_ϕ(ϕinit,rs)

    ## E-step 

    # Compute E[X̄|y_m] & E[X̄X̄|y_m] for m=1,2,…,M
    map_out = map(x->CpelNano.get_cond_exs(αs,βs,x),rs.calls)

    # Compute E[Sa(X̄)|y_1,…,y_m;ϕ], E[Sb(X̄)|y_1,…,y_m;ϕ], and E[Sc(X̄)|y_1,…,y_m;ϕ]
    ES = CpelNano.get_ess(map_out,rs)

    ## M-step

    # Set up system (requires to work w/ vectors)
    function f!(U::Vector{Float64},ϕ::Vector{Float64})
        ∇logZ = CpelNano.get_∇logZ(rs,ϕ)
        U[1] =  ∇logZ[1] - ES[1] / rs.m
        U[2] =  ∇logZ[2] - ES[2] / rs.m
        U[3] =  ∇logZ[3] - ES[3] / rs.m
    end
    
    # Solve NL system
    sol =
    try
        # Call solver
        NLsolve.nlsolve(f!,ϕinit;iterations=20,ftol=1e-3)
    catch x
        # Report x error if found
        print_log("An instance of EM algorithm did not converge.")
    end

    ## Convergence
    
    # Leave if no convergence with previous estimate
    (isnothing(sol) || !converged(sol)) && return ϕinit

    ## Project
    
    # Project vector to allowed hypercube
    θproj = sol.zero
    θproj[1] = max(-AMAX,min(θproj[1],AMAX))
    θproj[2] = max(-BMAX,min(θproj[2],BMAX))
    θproj[3] = max(-CMAX,min(θproj[3],CMAX))

    # Return new ϕ
    return θproj
    
end
"""
    `comp_Q_wgbs(REG_ST,ϕ)`

    Function that computes E[log p(X,y)|Y=y;ϕ] in WGBS case.

    # Examples
    ```julia-repl
    julia> M = 50; L = 80; αs = fill(0.0,L); βs = fill(0.0,L-1); pobs = 1.0;
    julia> rs = CpelNano.cpel_samp_wgbs(M,αs,βs,pobs);
    julia> rs.L = 10; rs.m = M; rs.Nl = fill(5.0,rs.L); rs.ρl = fill(0.1,rs.L); rs.dl = fill(10.0,rs.L-1);
    julia> CpelNano.comp_Q_wgbs(rs,zeros(3))
    -2.747243306144796
    ```
"""
function comp_Q_wgbs(rs::RegStruct,ϕ::Vector{Float64})::Float64
    
    # Get αs & βs
    αs,βs = CpelNano.get_αβ_from_ϕ(ϕ,rs)

    # Compute Ec[X|y_m] & Ec[XX|y_m] for m=1,2,…,M
    map_out = map(x->CpelNano.get_cond_exs(αs,βs,x),rs.calls)

    # Compute Ec[S(X)] & Ec[SS(X)]
    ES = CpelNano.get_ess(map_out,rs)

    # Add <ϕ,E[S]>
    Q = dot(ϕ,ES)

    # Matrices
    u1 = CpelNano.get_u(αs[1])
    uL = CpelNano.get_u(αs[end])
    Ws = [CpelNano.get_W(αs[l],αs[l+1],βs[l]) for l=1:length(βs)]

    # Scale Ws for numerical stability
    scaling = maximum(maximum.(Ws))
    Ws = UniformScaling(1.0/scaling) * Ws

    # Add - m⋅log[Z(ϕ)]
    Q -= rs.m * (log(CpelNano.get_Z(u1,uL,Ws)) + length(Ws)*log(scaling))

    # Return normalized Q
    return Q/(rs.L*rs.m)

end
"""

    `em_inst_wgbs(data,ϕ0,config)`
    
    Single instance of the EM algorithm when working with WGBS data.
    
    # Examples
    ```julia-repl
    julia> M = 10; L = 10; αs = fill(0.0,L); βs = fill(0.0,L-1); pobs = 1.0;
    julia> reg_dat = CpelNano.cpel_samp_wgbs(M,αs,βs,pobs);
    julia> ϕhat,Q = CpelNano.em_inst_wgbs(reg_dat,zeros(3),CpelNano.CpelNanoConfig())
    ([0.5210509995207954, 0.6242225806898237, 0.05607186596907174, 0.1085379770548646], -0.26316332446017654)
    ```
"""
function em_inst_wgbs(rs::RegStruct,ϕ0::Vector{Float64},config::CpelNanoConfig)::Tuple{Vector{Float64},Float64}

    # Initialize
    i = 1
    ϕt = ϕ0
    ϕtp1 = NaN
    conv = false
    max_iters_hit = false

    # Iterate
    while !conv && !max_iters_hit

        # Find next ϕ update
        ϕtp1 = em_iter_wgbs(ϕt,rs)
        
        # Print out Q
        config.verbose && print_log("ϕtp1=$(ϕtp1); Q(ϕ)=$(comp_Q_wgbs(rs,ϕtp1))")

        # Check for convergence
        conv = conv_check(ϕt,ϕtp1)
        
        # Check if too many iters
        max_iters_hit = i>=config.max_em_iters

        # Update ϕ
        ϕt = ϕtp1

        # Increase counter
        i += 1

    end
    
    # Check if more than 2 iterations (sometimes it gets stuck right at start)
    conv = conv && i>2

    # Get final pair
    Q = conv ? comp_Q_wgbs(rs,ϕtp1) : -Inf

    # Return pair (ϕ,Q)
    return ϕtp1,Q

end
"""

    `get_ϕhat_wgbs!(REG_ST,CONFIG)`
    
    Run multiple instances of the EM algorithm by randomly initializing the algorithm max_init times 
    when working with WGBS data.
    
    # Examples
    ```julia-repl
    julia> a = -0.5; b = 50.0; c = 5.0; ϕ = [a,b,c]; L = 100;
    julia> x=CpelNano.RegStruct(); x.L=L; x.Nl=fill(1.0,x.L);
    julia> x.ρl=rand(0.001:0.001:0.1,x.L); x.dl=fill(10.0,x.L-1); 
    julia> α,β = CpelNano.get_αβ_from_ϕ(ϕ,x);
    julia> M = 50; pobs = 0.25; rs = CpelNano.cpel_samp_wgbs(M,α,β,pobs); 
    julia> rs.L=x.L; rs.Nl=x.Nl; rs.ρl=x.ρl; rs.dl=x.dl; 
    julia> config = CpelNano.CpelNanoConfig(); config.max_em_iters = 100; 
    julia> config.max_em_init = 1; config.verbose = true;
    julia> CpelNano.get_ϕhat_wgbs!(rs,config); rs.ϕhat
    ```
"""
function get_ϕhat_wgbs!(rs::RegStruct,config::CpelNanoConfig)::Nothing
    
    # Initialize vector
    ϕ0s = gen_ϕ0s(config.max_em_init)

    # Launch multiple instances of EM in parallel
    em_out = map(ϕ0->em_inst_wgbs(rs,ϕ0,config),ϕ0s)

    # Choose ϕ with maximum expected loglike
    ϕhat,Qhat = pick_ϕ(em_out)
    config.verbose && print_log("Converged to Q=$(Qhat) and ϕ̂=$(ϕhat).")

    # Set ϕhat
    rs.ϕhat = Qhat>-Inf ? ϕhat : Vector{Float64}()

    # Return nothing
    return nothing

end