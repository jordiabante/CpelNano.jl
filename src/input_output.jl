##################################################################################################
# INPUT/OUTPUT
##################################################################################################
"""
    `print_log(MESSAGE)`

    Function that prints MESSAGE to stderr.

    # Examples
    ```julia-repl
    julia> CpelNano.print_log("Hello")
    [2020-03-30 16:24:18]: Hello
    ```
"""
function print_log(mess::String)

    # Right format for date
    date = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println(stderr,"[$(date)]: "  * mess)
    flush(stderr)

    # Return
    return nothing

end
"""
    `check_output_exists(CONFIG)`

    Function that checks if at least one output file exists.

    # Examples
    ```julia-repl
    julia> CpelNano.check_output_exists(config)
    ```
"""
function check_output_exists(config::CpelNanoConfig)::Bool
    
    # Init
    exists = isfile(config.out_files.mml_file)
    exists |= isfile(config.out_files.nme_file)
    exists |= isfile(config.out_files.theta_file)

    # Print out message if output exists
    if exists
        message = "At least an output file exists. Please, change the output directory and try again."
        print_log(message)
        sleep(5)
    end
    
    # Return bool
    return exists

end
"""
    `check_diff_output_exists(CONFIG)`

    Function that checks if at least one differential output file exists.

    # Examples
    ```julia-repl
    julia> CpelNano.check_diff_output_exists(config)
    ```
"""
function check_diff_output_exists(config::CpelNanoConfig)::Bool
    
    # Init
    exists = isfile(config.out_diff_files.tmml_file)
    exists |= isfile(config.out_diff_files.tnme_file)
    exists |= isfile(config.out_diff_files.tcmd_file)

    # Print out message if output exists
    if exists
        message = "At least an output file exists. Please, change the output directory and try again."
        print_log(message)
        sleep(5)
    end
    
    # Return bool
    return exists

end
##################################################################################################
# ESTIMATION OUTPUT
##################################################################################################
"""
    `write_output_ex(PATH,RS_VEC)`

    Function that stores file with E[X].

    # Examples
    ```julia-repl
    julia> CpelNano.write_output_ex(path,rs_vec)
    ```
"""
function write_output_ex(path::String,rs_vec::Vector{RegStruct})::Nothing

    # Write E[X]
    open(path,"a") do io
        for rs in rs_vec
            # Check if region analyzed
            rs.proc || continue
            # Loop over CpG groups
            @inbounds for l=1:rs.L
                isnan(rs.exps.ex[l]) && continue
                write(io,"$(rs.chr)\t$(rs.cpg_pos[l])\t$(rs.cpg_pos[l])\t$(rs.exps.ex[l])\n")
            end
        end
    end

    # Return
    return nothing

end
"""
    `write_output_exx(PATH,RS_VEC)`

    Function that stores file with E[XX].

    # Examples
    ```julia-repl
    julia> CpelNano.write_output_exx(path,rs_vec)
    ```
"""
function write_output_exx(path::String,rs_vec::Vector{RegStruct})::Nothing

    # Write E[XX]
    open(path,"a") do io
        for rs in rs_vec
            # Check if region analyzed
            rs.proc || continue
            # Loop over CpG groups
            @inbounds for l=1:(rs.L-1)
                isnan(rs.exps.exx[l]) && continue
                write(io,"$(rs.chr)\t$(rs.cpg_pos[l])\t$(rs.cpg_pos[l+1])\t$(rs.exps.exx[l])\n")
            end
        end
    end

    # Return
    return nothing

end
"""
    `write_output_mml(PATH,REGS_DATA)`

    Function that stores file with analysis region information μ(X)=[μ1(X),μ2(X),…,μK1(X)].

    # Examples
    ```julia-repl
    julia> CpelNano.write_output_mml(path,regs_data)
    ```
"""
function write_output_mml(path::String,regs_data::Vector{RegStruct})::Nothing

    # Write μ(X)
    open(path,"a") do io
        for rs in regs_data

            # Check if region analyzed
            rs.proc || continue

            # Loop over analysis region
            @inbounds for k=1:rs.nls_rgs.num

                # Check if data 
                isnan(rs.mml[k]) && continue

                # Get coordinates
                reg_st = minimum(rs.nls_rgs.chr_int[k])
                reg_end = maximum(rs.nls_rgs.chr_int[k])

                # Round mml to 4 decimals
                mml = round(rs.mml[k],digits=4)

                # Write μk(X) for k-th analysis region
                write(io,"$(rs.chr)\t$(reg_st)\t$(reg_end)\t$(mml)\n")

            end

        end
    end

    # Return
    return nothing

end
"""
    `write_output_nme(PATH,REGS_DATA)`

    Function that stores file with analysis region information h(X)=[h1(X),h2(X),…,hK1(X)].

    # Examples
    ```julia-repl
    julia> CpelNano.write_output_nme(path,regs_data)
    ```
"""
function write_output_nme(path::String,regs_data::Vector{RegStruct})::Nothing

    # Write h(X)
    open(path,"a") do io
        for rs in regs_data

            # Check if region analyzed
            rs.proc || continue

            # Loop over analysis region
            @inbounds for k=1:rs.nls_rgs.num

                # Check if data 
                isnan(rs.nme[k]) && continue 

                # Get coordinates
                reg_st = minimum(rs.nls_rgs.chr_int[k])
                reg_end = maximum(rs.nls_rgs.chr_int[k])

                # Round nme to 4 decimals
                nme = round(rs.nme[k],digits=4)

                # Write hk(X) for k-th analysis region
                write(io,"$(rs.chr)\t$(reg_st)\t$(reg_end)\t$(nme)\n")

            end

        end
    end

    # Return
    return nothing

end
"""
    `write_output_ϕ(PATH,RS_VEC)`

    Function that stores parameter vector ϕ as well as α & β.

    # Examples
    ```julia-repl
    julia> CpelNano.write_output_ϕ(path,rs_vec)
    ```
"""
function write_output_ϕ(path::String,rs_vec::Vector{RegStruct})::Nothing

    # Write ϕ
    open(path,"a") do io
        for rs in rs_vec

            # Continue if no data
            rs.proc || continue

            # Obtain each field
            ϕ = join(rs.ϕhat,',')

            # Obtain α & β
            α,β = get_αβ_from_ϕ(rs.ϕhat,rs)
            α = join(α,',')
            β = join(β,',')

            # Write to file
            write(io,"$(rs.chr)\t$(rs.chrst)\t$(rs.chrend)\t$(ϕ)\t$(α)\t$(β)\n")

        end
    end

    # Return
    return nothing

end
"""
    `write_output(REGS_DATA,CONFIG)`

    Function that writes output of `analyze_nano` to output files.

    # Examples
    ```julia-repl
    julia> CpelNano.write_output(regs_data,config)
    ```
"""
function write_output(regs_data::Vector{RegStruct},config::CpelNanoConfig)::Nothing

    # Write
    write_output_ϕ(config.out_files.theta_file,regs_data)
    write_output_exx(config.out_files.exx_file,regs_data)
    write_output_ex(config.out_files.ex_file,regs_data)
    write_output_mml(config.out_files.mml_file,regs_data)
    write_output_nme(config.out_files.nme_file,regs_data)
    
    # Return nothing
    return nothing
    
end
##################################################################################################
# DIFFERENTIAL ANALYSIS
##################################################################################################
"""
    `write_output_tmml(PATH,TEST_DATA)`

    Function that writes differential analysis output for statistic Tmml to output files.

    # Examples
    ```julia-repl
    julia> CpelNano.write_output_tmml(path,test_data)
    ```
"""
function write_output_tmml(path::String,test_data::Vector{RegStatTestStruct})::Nothing

    # Write
    open(path,"a") do io
        for ts in test_data
            # Loop over analysis region
            @inbounds for k=1:length(ts.coords)
                # Gather data
                chrst = minimum(ts.coords[k])
                chrend = maximum(ts.coords[k])
                tmml = ts.tests.tmml_test[k][1]
                pmml = ts.tests.tmml_test[k][2]
                # Skip if no CGs
                isnan(tmml) && continue
                # Round Tmml to 4 decimals
                tmml = round(tmml,digits=4)
                # Write for k-th analysis region
                write(io,"$(ts.chr)\t$(chrst)\t$(chrend)\t$(tmml)\t$(pmml)\n")
            end
        end
    end
    
    # Return nothing
    return nothing

end
"""
    `write_output_tnme(PATH,TEST_DATA)`

    Function that writes differential analysis output for statistic Tnme to output files.

    # Examples
    ```julia-repl
    julia> CpelNano.write_output_tnme(path,test_data)
    ```
"""
function write_output_tnme(path::String,test_data::Vector{RegStatTestStruct})::Nothing

    # Write
    open(path,"a") do io
        for ts in test_data
            # Loop over analysis region
            @inbounds for k=1:length(ts.coords)
                # Gather data
                chrst = minimum(ts.coords[k])
                chrend = maximum(ts.coords[k])
                tnme = ts.tests.tnme_test[k][1]
                pnme = ts.tests.tnme_test[k][2]
                # Skip if no CGs
                isnan(tnme) && continue
                # Round Tnme to 4 decimals
                tnme = round(tnme,digits=4)
                # Write for k-th analysis region
                write(io,"$(ts.chr)\t$(chrst)\t$(chrend)\t$(tnme)\t$(pnme)\n")
            end
        end
    end
    
    # Return nothing
    return nothing

end
"""
    `write_output_tcmd(PATH,TEST_DATA)`

    Function that writes differential analysis output for statistic Tcmd to output files.

    # Examples
    ```julia-repl
    julia> CpelNano.write_output_tcmd(path,test_data)
    ```
"""
function write_output_tcmd(path::String,test_data::Vector{RegStatTestStruct})::Nothing

    # Write
    open(path,"a") do io
        for ts in test_data
            # Loop over analysis region
            @inbounds for k=1:length(ts.coords)
                # Gather data
                chrst = minimum(ts.coords[k])
                chrend = maximum(ts.coords[k])
                tcmd = ts.tests.tcmd_test[k][1]
                pcmd = ts.tests.tcmd_test[k][2]
                # Skip if no CGs
                isnan(tcmd) && continue
                # Round Tcmd to 4 decimals
                tcmd = round(max(0.0,tcmd),digits=4)
                # Write for k-th analysis region
                write(io,"$(ts.chr)\t$(chrst)\t$(chrend)\t$(tcmd)\t$(pcmd)\n")
            end
        end
    end
    
    # Return nothing
    return nothing

end
"""
    `write_diff_output(REGS_DATA,CONFIG)`

    Function that writes differential analysis output to output files.

    # Examples
    ```julia-repl
    julia> CpelNano.write_diff_output(regs_data,config)
    ```
"""
function write_diff_output(test_data::Vector{RegStatTestStruct},config::CpelNanoConfig)::Nothing

    # Write
    write_output_tmml(config.out_diff_files.tmml_file,test_data)
    write_output_tnme(config.out_diff_files.tnme_file,test_data)
    write_output_tcmd(config.out_diff_files.tcmd_file,test_data)
    
    # Return nothing
    return nothing

end
##################################################################################################
# READ IN
##################################################################################################
"""
    `parse_unit_range(STRING)`

    Function that parses a unit range in the form of a string.

    # Examples
    ```julia-repl
    julia> CpelNano.parse_unit_range("1:2")
    1:2
    ```
"""
function parse_unit_range(x::String)
    
    # Split string using the colon
    y = split(x,':')

    # Return unit range
    return parse(Int,String(y[1])):parse(Int,String(y[2]))

end
"""
    `read_bed_reg(CONFIG)`

    Function that returns a dictionary with regions of interest from BED file provided. We expect a BED 
    file with chr-reg_st-reg_end as fields, separated by tabulations. We do not expect a header in the
    file.

    # Examples
    ```julia-repl
    julia> CpelNano.read_bed_reg(config)
    ```
"""
function read_bed_reg(config::CpelNanoConfig)::Dict{String,Vector{UnitRange{Int64}}}

    # Initizalize
    targ_regs = Dict{String,Vector{UnitRange{Int64}}}()
    open(config.bed_reg) do f
        for line in enumerate(eachline(f))
            # Get line info
            line_vec = split(line[2],"\t")
            # Get chromosome
            line_chr = String(line_vec[1])
            # Get start/end
            line_chrst = parse(Int64,line_vec[2])
            line_chrend = parse(Int64,line_vec[3])
            line_chrend >= line_chrst || continue
            # print_log("line_chr: $(line_chr); line_chrst: $(line_chrst); line_chrend: $(line_chrend)")
            # Add to dictionary
            if haskey(targ_regs,line_chr)
                push!(targ_regs[line_chr],line_chrst:line_chrend)
            else
                targ_regs[line_chr] = [line_chrst:line_chrend]
            end
            # print_log("targ_regs: $(targ_regs)")
        end
    end

    # Return dictionary
    return targ_regs

end
"""
    `get_dic_nanopolish(REG,CALL_FILE)`

    Function that returns reads in methylation calls by nanopolish and creates a dictionary of data. This
    function assumes that the file is sorted using the chromosome and start coordinate of the CpG groups, and
    it assumes that the file contains a header. Furthermore, this function assumes the following fields: 
    1. chromosome; 3. start; 4. end; 5. read_name.

    # Examples
    ```julia-repl
    julia> CpelNano.get_dic_nanopolish(reg,call_file)
    ```
"""
function get_dic_nanopolish(reg_data::RegStruct,call_file::String)::Dict{String,Vector{Vector{SubString}}}

    # Read in (assume it's sorted)
    chr_visited = false
    line_dic = Dict{String,Vector{Vector{SubString}}}()
    open(call_file) do f
        for line in enumerate(eachline(f))
            # Skip header
            line[1]>1 || continue
            # Get line info
            line_vec = split(line[2],"\t")
            # Get chromosome
            line_chr = line_vec[1]
            line_chr == reg_data.chr || continue
            # Set flag
            chr_visited = true
            # Get start
            line_chrst = parse(Int64,line_vec[3])
            line_chrst <= reg_data.chrend || continue
            # Get end
            line_chrend = parse(Int64,line_vec[4])
            line_chrend >= reg_data.chrst || continue
            # Get ID
            line_id = line_vec[5]
            # Add to dictionary
            haskey(line_dic,line_id) ? append!(line_dic[line_id],[line_vec]) : line_dic[line_id] = [line_vec]
            # Break if done: checked the right chromosome and beyond end position
            chr_visited && line_chrst>reg_data.chrend && break 
        end
    end
    
    # Return call dictionary
    return line_dic

end
"""
    `read_mod_file_chr(FILE,CHR)`

    Function that reads in the output of `analyze_nano` containing models in chromosome CHR.

    # Examples
    ```julia-repl
    julia> CpelNano.read_mod_file_chr(file,chr)
    ```
"""
function read_mod_file_chr(file::String,chr::String,fa_rec::FASTA.Record,config::CpelNanoConfig)::Dict{String,RegStruct}

    # Go over all lines in file
    out_dic = Dict{String,RegStruct}()
    open(file) do f
        for line in enumerate(eachline(f))

            # Get line data
            line_vec = split(line[2],"\t")

            # Get chromosome info
            reg_chr = line_vec[1]
            reg_chr == chr || continue

            # Init struct
            rs = RegStruct()
            rs.chr = reg_chr
            
            # Get start or end
            rs.chrst = parse(Int64,line_vec[2])
            rs.chrend = parse(Int64,line_vec[3])
            
            # Get a,b,c
            rs.ϕhat = parse.(Float64,String.(split(line_vec[4],",")))

            # Get genomic info
            get_grp_info!(rs,fa_rec,config.min_grp_dist)

            # Re-scale model
            rscle_grp_mod!(rs)
            
            # Statistical summaries
            get_stat_sums!(rs,config)

            # Mark region as processed
            rs.proc = true

            # Set region ID based on chr_st_end
            reg_id = join([rs.chr,rs.chrst,rs.chrend],"_")

            # Add to dictionary 
            out_dic[reg_id] = rs
            
        end
    end
    
    # Return dictionary
    return out_dic

end
##################################################################################################
# SORT BEDGRAPH
##################################################################################################
"""
    `sort_bedGraph(BEDGRAPH)`

    Function that sorts bedGraph file.

    # Examples
    ```julia-repl
    julia> CpelNano.sort_bedGraph(bedGraph)
    ```
"""
function sort_bedGraph(bedGraph::String)::Nothing

    # Return nothing
    return nothing

end