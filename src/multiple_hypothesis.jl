###################################################################################################
# MULTIPLE HYPOTHESIS CORRECTION
###################################################################################################
"""
    `add_qvalues_to_path(PATH)`

    Function that takes in the differential output and adds BH q-value.

    # Examples
    ```julia-repl
    julia> CpelNano.add_qvalues_to_path(path)
    ```
"""
function add_qvalues_to_path(path::String)::Nothing

    # Leave if no data
    filesize(path) > 0 || return nothing
    
    # Get data
    all_data = readdlm(path, '\t', Any)
    qvals = fill(NaN, size(all_data)[1])
    
    # Multiple hypothesis testing correction
    ind = .!isnan.(all_data[:,5])
    if sum(ind) > 0 
        qvals[ind] = MultipleTesting.adjust(convert(Vector{Float64}, all_data[ind,5]), BenjaminiHochberg())
    end
    
    # Append to output matrix
    all_data = hcat(all_data, qvals)
    
    # Write to temp output
    temp_path = path * ".tmp"
    open(temp_path, "w") do io
        writedlm(io, all_data, '\t')
    end

    # Move to original file
    mv(temp_path, path, force=true)

    # Return
    return nothing

end
"""
    `mult_hyp_corr(CONFIG)`

    Function that takes in all the differential output and adds BH q-value in each one.

    # Examples
    ```julia-repl
    julia> CpelNano.mult_hyp_corr(config)
    ```
"""
function mult_hyp_corr(config::CpelNanoConfig)::Nothing

    # Add q-values
    add_qvalues_to_path(config.out_diff_files.tmml_file)
    add_qvalues_to_path(config.out_diff_files.tnme_file)
    config.matched || add_qvalues_to_path(config.out_diff_files.tcmd_file)

    # Return
    return nothing

end