"""

    `comp_smp_ex!(RS)`
    
    Stores sample average of E[X]. 
    
    # Examples
    ```julia-repl
    julia> CpelNano.comp_smp_ex!(rs)
    ```
"""
function comp_smp_ex!(rs::RegStruct)::Nothing
    
    # Loop over groups
    ex = fill(NaN,rs.N)
    @inbounds for l=1:rs.L
        aux = 0.0
        n_data = 0.0
        # Loop over calls
        @inbounds for m=1:length(rs.calls)
            # Get call
            call = rs.calls[m][l]
            if call.obs
                # If observed add contribution
                aux += call.log_pyx_m>call.log_pyx_u ? 1.0 : -1.0 
                n_data += 1.0
            end
        end
        # Store if any data point
        if n_data>0.0
            ex[rs.cpg_grps[l].cpg_ind] .= aux/n_data
        end
    end

    # Store in struct
    rs.exps.ex = ex

    # Return
    return nothing

end
"""

    `comp_smp_exx!(RS)`
    
    Stores sample average of E[XX]. 
    
    # Examples
    ```julia-repl
    julia> CpelNano.comp_smp_exx!(rs)
    ```
"""
function comp_smp_exx!(rs::RegStruct)::Nothing
    
    # Loop over groups
    exx = fill(NaN,rs.N-1)
    @inbounds for l=1:(rs.L-1)
        aux = 0.0
        n_data = 0.0
        # Loop over calls
        @inbounds for m=1:length(rs.calls)
            # Get call
            call1 = rs.calls[m][l]
            call2 = rs.calls[m][l+1]
            # Interaction between groups
            if call1.obs && call2.obs
                # If observed add contribution
                aux1 = call1.log_pyx_m>call1.log_pyx_u ? 1.0 : -1.0
                aux2 = call2.log_pyx_m>call2.log_pyx_u ? 1.0 : -1.0
                aux += aux1*aux2
                # Increase counter
                n_data += 1.0
            end
        end
        # Store between groups if any data point
        if n_data>0.0
            exx[rs.cpg_grps[l].cpg_ind[end]] = aux/n_data
        end
        # Store within groups if not singleton
        if length(rs.cpg_grps[l].cpg_ind)>1
            exx[rs.cpg_grps[l].cpg_ind[1:(end-1)]] .= 1.0
        end
    end

    # Take care of last group if not singleton
    if length(rs.cpg_grps[end].cpg_ind)>1
        exx[rs.cpg_grps[end].cpg_ind[1:(end-1)]] .= 1.0
    end

    # Store in struct
    rs.exps.exx = exx

    # Return
    return nothing

end
"""

    `comp_smp_exps!(RS)`
    
    Stores sample averages of E[X] and E[XX]. 
    
    # Examples
    ```julia-repl
    julia> CpelNano.comp_smp_exps!(rs)
    ```
"""
function comp_smp_exps!(rs::RegStruct)::Nothing
    
    # Compute expectations
    comp_smp_ex!(rs)
    comp_smp_exx!(rs)

    # Return
    return nothing

end