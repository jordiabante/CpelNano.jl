"""
    `get_αβ_from_ϕ(ϕ,REG_STRUCT)`

    Get full α and β vectors from vector ϕ=[a,b,c] using Nl, ρ̄l, and dl, assuming energy function
    from Jenkinson et al 2017 and 2018.

        α_l = (a + b⋅ρ̄_l)⋅N_l , β_l = c/d_l

    # Examples
    ```julia-repl
    julia> rs=CpelNano.RegStruct(); rs.L=10; rs.Nl=fill(5.0,rs.L); rs.ρl=fill(0.1,rs.L); rs.dl=fill(10.0,rs.L-1); 
    julia> CpelNano.get_αβ_from_ϕ([0.0,0.0,0.0],rs)
    ([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]) 
    ```
"""
function get_αβ_from_ϕ(ϕ::Vector{Float64},rs::RegStruct)::NTuple{2,Vector{Float64}}

    # Obtain full α & β vectors
    α = [rs.Nl[l]*(ϕ[1]+ϕ[2]*rs.ρl[l]) for l=1:rs.L]
    β = [ϕ[3]/rs.dl[l] for l=1:(rs.L-1)]

    # Return α and β vector
    return α,β
    
end