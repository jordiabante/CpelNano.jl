"""
    `est_marg_model(calls)`
    
    Function to estimate marginal model from without accounting for error.
    
    # Examples
    ```julia-repl
    julia>
    ```
"""
function est_marg_model(calls::Vector{Vector{MethCallCpgGrp}})::Vector{Float64}
   
    # Get total number of CpG sites
    N = length(calls[1])

    # For each call 
    pvec = fill(0.0, N)
    nvec = fill(0.0, N)
    for call in calls
        for i = 1:length(call)
            if call[i].obs
                pvec[i] += call[i].log_pyx_m - call[i].log_pyx_u >= 0 ? 1 : 0
                nvec[i] += 1
            end
        end
    end

    # Return phat vector
    return pvec ./ nvec
 
end
"""
    `comp_mml_marg_model(pvec,𝒜s)`
    
    Function to compute average MML within analysis regions.
    
    # Examples
    ```julia-repl
    julia>
    ```
"""
function comp_mml_marg_model(pvec::Vector{Float64}, 𝒜s::Vector{UnitRange{Int64}})::Vector{Float64}
   
    # Return MML
    return [sum(pvec[x]) / length(x) for x in 𝒜s]
 
end
"""
    `comp_nme_marg_model(pvec,𝒜s)`
    
    Function to compute average NME within analysis regions.
    
    # Examples
    ```julia-repl
    julia>
    ```
"""
function comp_nme_marg_model(pvec::Vector{Float64}, 𝒜s::Vector{UnitRange{Int64}})::Vector{Float64}
   
    # Compute H(X) = -p⋅log(p) + (1-p)⋅log(p)
    cpg_nme = - 1.0 / (LOG2) * [x * log(max(x, eps(Float64))) + (1.0 - x) * log(max(1.0 - x, eps(Float64))) for x in pvec]

    # Return NME
    return [sum(cpg_nme[x]) / length(x) for x in 𝒜s]
 
end