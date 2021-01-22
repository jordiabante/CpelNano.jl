##################################################################################################
# Differential analysis step 1 - Tpdm
##################################################################################################
"""
    `get_all_models_chr(MODEL_FILES,CHR)`

    Function that takes reads in all models in MODEL_FILES found in chr CHR.

    # Examples
    ```julia-repl
    julia> CpelNano.get_all_models_chr(mod_files,chr)
    ```
"""
function get_all_models_chr(mod_files::Vector{String}, chr::String)::Vector{Dict{String,RegStruct}}

    # Get all data from chromosome
    mods = Vector{Dict{String,RegStruct}}()
    @inbounds for f in mod_files
        push!(mods, read_model_file_chr(f, chr))
    end

    # Return models
    return mods

end
"""
    `get_unique_ids(MODELS)`

    Function that returns unique analysis regions IDs.

    # Examples
    ```julia-repl
    julia> CpelNano.get_unique_ids(mods)
    ```
"""
function get_unique_ids(mods::Vector{Dict{String,RegStruct}})::Vector{String}

    # Loop over all samples to get list of unique IDs
    uniq_ids = Vector{String}()
    @inbounds for s in mods
        @inbounds for key in keys(s)
            found = findfirst(isequal(key), uniq_ids)
            isnothing(found) && push!(uniq_ids, key)
        end
    end

    # Return models
    return uniq_ids

end
"""
    `get_mods_reg(REG_ID,MODELS)`

    Function that takes model files and returns vector with models from REG_ID.

    # Examples
    ```julia-repl
    julia> CpelNano.get_mods_reg(reg_id,mods_g1)
    ```
"""
function get_mods_reg(reg_id::String, mods::Vector{Dict{String,RegStruct}})::Vector{RegStruct}

    # Find models with data
    mods_reg = Vector{RegStruct}()
    @inbounds for mod in mods
        
        # Check if sample has region analyzed
        if reg_id in keys(mod) 
            # Push model
            push!(mods_reg, mod[reg_id])
        else
            # Push empty model
            push!(mods_reg, RegStruct())
        end

    end

    # Return models from region
    return mods_reg

end
"""
    `check_samples_for_testing(MODELS_G1,MODELS_G2,MATCHED)`

    Function that takes looks at the samples, or pairs in the matched setting, available for 
    doing the statistical test.

    # Examples
    ```julia-repl
    julia> CpelNano.check_samples_for_testing(ms_g1,ms_g2,matched)
    ```
"""
function check_samples_for_testing(ms_g1::Vector{RegStruct},ms_g2::Vector{RegStruct},
    matched::Bool)::Tuple{Vector{RegStruct},Vector{RegStruct},Bool}

    # Check which models are available
    proc_1 = [m.proc for m in ms_g1]
    proc_2 = [m.proc for m in ms_g2]

    # Get samples/pairs
    if matched
        
        # Find pairs with data
        ind = proc_1 .& proc_2

        # Do test if p-values can be below 0.05
        test_go = (0.5^(sum(ind) - 1)) < 0.05 
        
        # Get sample pairs with data
        ms_g1 = ms_g1[ind]
        ms_g2 = ms_g2[ind]
        
    else
        
        # Check if sufficient data
        s1 = sum(proc_1)
        s2 = sum(proc_2)

        # Do test if p-values can be below 0.05
        test_go = 2.0 / binomial(s1 + s2, s1) < 0.05
        
        # Get samples with data
        ms_g1 = ms_g1[proc_1]
        ms_g2 = ms_g2[proc_2]
        
    end

    # Return tuple with models
    return ms_g1, ms_g2, test_go
    
end
"""
    `pmap_diff_analysis_reg(REG_ID,MODELS_G1,MODELS_G2,CONFIG)`

    Function that takes model files for groups 1 and 2 and performs differential analysis at the
    analysis region level in analysis region with REG_ID.

    # Examples
    ```julia-repl
    julia> CpelNano.differential_analysis(reg_id,mods_g1,mods_g2,config)
    ```
"""
function pmap_diff_analysis_reg(reg_id::String,mods_g1::Vector{Dict{String,RegStruct}},
    mods_g2::Vector{Dict{String,RegStruct}},config::CpelNanoConfig)::NTuple{2,Float64}

    # Find models with data for both groups
    ms_g1 = get_mods_reg(reg_id, mods_g1)
    ms_g2 = get_mods_reg(reg_id, mods_g2)

    # Check which models are available
    ms_g1, ms_g2, test_go = check_samples_for_testing(ms_g1, ms_g2, config.matched)
    config.min_pval && !test_go && return (NaN, NaN)

    # Return struct with test results
    return config.matched ? mat_reg_test_tpdm(ms_g1, ms_g2) : unmat_reg_test_tpdm(ms_g1, ms_g2)

end
"""
    `test_diff_analysis_subreg!(TEST_OUT,MODELS_G1,MODELS_G2,MATCHED)`

    Function that performs testing in each subregion.

    # Examples
    ```julia-repl
    julia> CpelNano.test_diff_analysis_subreg!(test_out,ms_g1,ms_g2,matched)
    ```
"""
function test_diff_analysis_subreg!(test_out::RegStatTestStruct,ms_g1::Vector{RegStruct},
    ms_g2::Vector{RegStruct},matched::Bool)::Nothing

    # For each subregion, do test
    k = 1
    @inbounds for i = 1:length(test_out.nls_reg_cpg_occ)
        
        # If no occupancy continue to next α-subregion
        test_out.nls_reg_cpg_occ[i] || continue
        
        # Perform hypothesis testing
        if matched
            # Push matched test
            push!(test_out.nls_reg_tests, mat_nls_reg_test(ms_g1, ms_g2, k))
        else
            # Push unmatched test
            push!(test_out.nls_reg_tests, unmat_nls_reg_test(ms_g1, ms_g2, k))
        end

        # Increase counter
        k += 1
        
    end
    
    # Return nothing
    return nothing

end
"""
    `pmap_diff_analysis_subreg(REG_ID,MODELS_G1,MODELS_G2,CONFIG)`

    Function that takes model files for groups 1 and 2 and performs differential analysis at the
    analysis region level in analysis region with REG_ID.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_diff_analysis_subreg(reg_id,mods_g1,mods_g2,config)
    ```
"""
function pmap_diff_analysis_subreg(reg_id::String,mods_g1::Vector{Dict{String,RegStruct}},
    mods_g2::Vector{Dict{String,RegStruct}},config::CpelNanoConfig)::RegStatTestStruct

    # Find models with data for both groups
    ms_g1 = get_mods_reg(reg_id, mods_g1)
    ms_g2 = get_mods_reg(reg_id, mods_g2)

    # Check which models are available
    ms_g1, ms_g2, test_go = check_samples_for_testing(ms_g1, ms_g2, config.matched)

    # Init output structure
    test_out = RegStatTestStruct(reg_id, ms_g1[1].cpg_occ)

    # Return empty test if not enough samples
    config.min_pval && !test_go && return test_out
    
    # For each subregion, do test
    test_diff_analysis_subreg!(test_out, ms_g1, ms_g2, config.matched)
    
    # Return permutation test output
    return test_out

end
"""
    `get_out_paths(CONFIG)`

    Get output paths for differential files.

    # Examples
    ```julia-repl
    julia> CpelNano.get_out_paths(config)
    ```
"""
function get_out_paths(config::CpelNanoConfig)::NTuple{3,String}

    # Output files
    tpdm_file = "$(config.out_dir)/$(config.out_prefix)_tpdm_reg.bedGraph"
    tmml_file = "$(config.out_dir)/$(config.out_prefix)_tmml_subreg.bedGraph"
    tnme_file = "$(config.out_dir)/$(config.out_prefix)_tnme_subreg.bedGraph"
    
    # Return permutation test output
    return tpdm_file, tmml_file, tnme_file

end
"""
    `get_filt_regs(PMAP_OUT,UNIQUE_IDS,MATCHED)`

    Select analysis regions for further testing.

    # Examples
    ```julia-repl
    julia> CpelNano.get_filt_regs(pmap_out,uniq_ids,matched)
    ```
"""
function get_filt_regs(pmap_out::Vector{NTuple{2,Float64}}, uniq_ids::Vector{String}, matched::Bool)::Vector{String}

    # Get significant regions
    pass_regs = Vector{String}()
    @inbounds for i = 1:length(uniq_ids)

        # Push if necessary
        if matched
            # Push if Tpdm ≥ Tmin
            pmap_out[i][1] >= 0.05 && push!(pass_regs, uniq_ids[i])
        else
            # Push if significant
            pmap_out[i][2] <= 0.05 && push!(pass_regs, uniq_ids[i])
        end
        
    end
    
    # Return significant regions ID
    return pass_regs

end
"""
    `add_qvalues_to_path(PATH)`

    Function that takes in the differential output and adds BH q-value.

    # Examples
    ```julia-repl
    julia> CpelAsm.add_qvalues_to_path(path)
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
    `mult_hyp_corr(DIFF_PATHS)`

    Function that takes in all the differential output and adds BH q-value in each one.

    # Examples
    ```julia-repl
    julia> CpelAsm.mult_hyp_corr(diff_paths)
    ```
"""
function mult_hyp_corr(diff_paths::Vector{String})

    # Add q-values
    add_qvalues_to_path(diff_paths[1])
    add_qvalues_to_path(diff_paths[2])
    # add_qvalues_to_path(diff_paths[3])

    # Return
    return nothing

end
"""
    `differential_analysis(MODEL_FILES_G1,MODEL_FILES_G2,FASTA,CONFIG)`

    Function that takes model files for groups 1 and 2 and performs differential analysis at the
    analysis region level.

    # Examples
    ```julia-repl
    julia> CpelNano.differential_analysis(mod_files_g1,mod_files_g2,fasta,config)
    ```
"""
function differential_analysis(mod_files_g1::Vector{String}, mod_files_g2::Vector{String}, fasta::String, config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names, chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)
    
    # Output files
    tpdm_file, tmml_file, tnme_file = get_out_paths(config)

    # Loop over chromosomes
    for i = 1:length(chr_names)

        # Get chromosome name and size
        chr = chr_names[i]
        chr_size = chr_sizes[i]

        # Print info
        print_log("Processing chr: $(chr) ...")
        print_log("Reading in data from chr $(chr) ...")

        ## Deal with analysis regions
        
        # Get all data from chromosome
        mods_g1 = get_all_models_chr(mod_files_g1, chr)
        mods_g2 = get_all_models_chr(mod_files_g2, chr)

        # Find unique IDs
        print_log("Finding unique region IDs in chr $(chr) ...")
        uniq_ids = unique(vcat(get_unique_ids(mods_g1), get_unique_ids(mods_g2)))

        # Process each analysis region in chr
        print_log("Performing differential PDM analysis in chr $(chr) ...")
        pmap_out = pmap(x -> pmap_diff_analysis_reg(x, mods_g1, mods_g2, config), uniq_ids)
            # print_log("$(pmap_out)")

        # Add data to respective bedGraph file
        print_log("Writing differential PDM analysis results from chr $(chr) ...")
        write_reg_diff_output(pmap_out, uniq_ids, config)

        # Keep only IDs that pass filter if desired
        print_log("Identifying significant PDM regions in chr $(chr) ...")
        uniq_ids = config.filter ? get_filt_regs(pmap_out, uniq_ids, config.matched) : uniq_ids
            # print_log("$(uniq_ids)")
        
        # Clean pmap_out
        pmap_out = Vector{RegStatTestStruct}()

        ## Deal with subregion differential analysis

        # Process each subregion within sig analysis region
        if length(uniq_ids) > 0
            print_log("Performing differential MML & NME analysis in chr $(chr) ...")
            pmap_out = pmap(x -> pmap_diff_analysis_subreg(x, mods_g1, mods_g2, config), uniq_ids)
                # print_log("$(pmap_out)")
        else
            print_log("No analysis regions to look do MML & NME analysis in chr $(chr) ...")
        end

        # Add last to respective bedGraph file
        print_log("Writing differential MML & NME analysis results from chr $(chr) ...")
        write_nls_reg_diff_output(pmap_out, uniq_ids, config)

    end

    # Multiple hypothesis testing
    print_log("Performing multiple hypothesis correction ...")
    mult_hyp_corr([tmml_file,tnme_file])

    # Sort bedGraph
    # sort_bedGraph(tpdm_file)
    # sort_bedGraph(tmml_file)
    # sort_bedGraph(tnme_file)
    
    # Return nothing
    return nothing

end