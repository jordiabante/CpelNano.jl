
"""
    `get_nano_pore_mod(NANO_MOD_FILE)`

    Function that stores the signal to noise ratio for each k-mer in a tabulated file with columns

        | Unmethylated k-mer | Signal-to-noise Ratio |

    # Examples
    ```julia-repl
    julia> CpelNano.get_nano_pore_mod(nano_mod_file)
    ```
"""
function get_nano_pore_mod(nano_mod_file::String)::Dict{String,NTuple{4,Float64}}

    # Read in model file
    nano_mod = readdlm(nano_mod_file)

    # Create dictionary
    nano_dic = Dict{String,NTuple{4,Float64}}()
    for i = 2:size(nano_mod)[1]
        kmer = nano_mod[i,1]
        occursin('C', kmer) || occursin('M', kmer) || continue
        if occursin('M', kmer)
            kmer_mod = replace(kmer, "M" => "C")
            # print_log("kmer: $(kmer) => kmer_mod: $(kmer_mod)")
            if kmer_mod in keys(nano_dic)
                μu = nano_dic[kmer_mod][1]
                σu = nano_dic[kmer_mod][2]
                nano_dic[kmer_mod] = (μu, σu, nano_mod[i,2], nano_mod[i,3])
            else
                nano_dic[kmer_mod] = (NaN, NaN, nano_mod[i,2], nano_mod[i,3])
            end
        else
            if kmer in keys(nano_dic)
                μm = nano_dic[kmer][3]
                σm = nano_dic[kmer][4]
                nano_dic[kmer] = (nano_mod[i,2], nano_mod[i,3], μm, σm)
            else
                nano_dic[kmer] = (nano_mod[i,2], nano_mod[i,3], NaN, NaN)
            end
        end
    end

    # Return dictionary
    return nano_dic

end
"""
    `pmap_split_nanopolish_file(NANO_FILE,N_FILES,HEADER,PREF,READ_FILE_IND,READ_NAME_VEC,IND_FILE)`

    Function writes one of the files in `split_nanopolish_file`.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_split_nanopolish_file(nano_file,n_files,header,pref,read_file_ind,read_name_vec,ind_file)
    ```
"""
function pmap_split_nanopolish_file(nano_file::String, n_files::Int64, header::String, pref::SubString{String}, read_file_ind::Vector{Int64}, read_nm_vec::Vector{String}, ind_file::Int64)::Nothing

    # Name of file
    nth_file = "$(pref).rep_$(ind_file).tsv"
    print_log("Writing file $(ind_file) out of $(n_files)...")
    # Write individual file
    open(nth_file, "w") do io
        # Add header
        write(io, "$(header)\n")
        # Add reads
        open(nano_file) do f
            for line in enumerate(eachline(f))
                # Skip header
                line[1] == 1 && continue
                # Get ID
                line_id = split(line[2], "\t")[5]
                # Get index of read in permuted vector
                read_ind_vec = findfirst(x -> x == line_id, read_nm_vec)
                # If read belongs to this file print
                if read_file_ind[read_ind_vec] == ind_file
                    write(io, "$(line[2])\n")
                end
            end
        end
    end

    # Return nothing
    return nothing

end
"""
    `split_nanopolish_file(NANO_FILE,N_FILES)`

    Function that splits nanopolish file into `N_FILES` files.

    # Examples
    ```julia-repl
    julia> CpelNano.split_nanopolish_file(nano_file,n_files)
    ```
"""
function split_nanopolish_file(nano_file::String, n_files::Int64)::Nothing

    # Prefix output file
    pref = split(nano_file, ".tsv")[1]

    print_log("Getting read IDs...")

    # Obtain all read names
    header = ""
    read_nm_vec = fill("", 200000)
    open(nano_file) do f
        next_entry = 1
        for line in enumerate(eachline(f))
            # Record and skip header
            if line[1] == 1
                header =  line[2]
                continue
            end
            # Get ID
            line_id = split(line[2], "\t")[5]
            # Add if necessary
            min_ind = max(1, next_entry - 100)
            ind_eq = findfirst(isequal(line_id), read_nm_vec[min_ind:next_entry])
            if isnothing(ind_eq)
                read_nm_vec[next_entry] = line_id
                mod(next_entry, 1000) == 0 && println(next_entry)
                next_entry += 1
            end
        end
    end
    
    # Remove empty
    deleteat!(read_nm_vec, read_nm_vec .== "")

    # Remove duplicates
    read_nm_vec = unique(read_nm_vec)

    # Number of reads per file
    base_num = div(length(read_nm_vec), n_files)
    extra_num = mod(length(read_nm_vec), n_files)
    
    print_log("Permuting reads...")

    # Permute order of reads so to not create bias
    permute!(read_nm_vec, sample(1:length(read_nm_vec), length(read_nm_vec);replace=false))

    print_log("Randomly assigning read IDs to files...")

    # Assign read to each file
    read_file_ind = Vector{Int64}()
    file_bnds = [i * base_num + extra_num for i = 1:n_files]
    @inbounds for i = 1:length(read_nm_vec)
        push!(read_file_ind, findfirst(x -> i <= x, file_bnds))
    end

    # Write output files
    pmap(x -> pmap_split_nanopolish_file(nano_file, n_files, header, pref, read_file_ind, read_nm_vec, x), 1:n_files)

    # Return nothing
    return nothing

end
"""
    `get_calls_nanopolish(CALL_DIC,REG_DATA)`

    Function that using the takes overlapping the analysis region and extracts information for each CpG 
    group.

    # Examples
    ```julia-repl
    julia> CpelNano.get_calls_nanopolish(call_dic,rs)
    ```
"""
function get_calls_nanopolish!(call_dic::Dict{String,Vector{Vector{SubString}}}, rs::RegStruct, mod_nse::Bool)::Nothing

    # For each read in call_dic
    calls = Vector{Vector{MethCallCpgGrp}}()
    @inbounds for kk in keys(call_dic)

        # Initizalize vector of group calls for the read
        call = [MethCallCpgGrp() for l = 1:length(rs.cpg_grps)]
        
        @inbounds for cpg_grp_dat in call_dic[kk]
            
            # Get chr_start (0-based in nanopolish)
            # Note: nanopolish does not consider the ±5 flanking bp
            cpg_grp_st = parse(Int64, cpg_grp_dat[3]) + 1 - 5
            cpg_grp_end = parse(Int64, cpg_grp_dat[4]) + 1 + 5
            
            # Find correct index of CpG group in region
            x_index_st = findfirst(x -> minimum(x.grp_int) == cpg_grp_st, rs.cpg_grps)
            x_index_end = findfirst(x -> maximum(x.grp_int) == cpg_grp_end, rs.cpg_grps)
            (isnothing(x_index_st) || isnothing(x_index_end)) && continue
                # print_log("x_index_st: $(x_index_st)")
                
            # Check if there is a discrepancy in CpG group coordinates
            if x_index_st != x_index_end
                print_log("[WARNING] CpG group coordinates do not match nanopolish output ... ")
                continue
            end
            
            # Record log p(y|x̄)
            log_pyx_m = parse(Float64, cpg_grp_dat[7])
            log_pyx_u = parse(Float64, cpg_grp_dat[8])
            if mod_nse
                # Modeling noise uses log p's given by nanopolish
                call[x_index_st].log_pyx_m = log_pyx_m
                call[x_index_st].log_pyx_u = log_pyx_u
            else
            # Not modeling noise does not use log p's given by nanopolish
                call[x_index_st].log_pyx_m = log_pyx_m > log_pyx_u ? log_pyx_right_x : log_pyx_wrong_x
        call[x_index_st].log_pyx_u = log_pyx_m > log_pyx_u ? log_pyx_wrong_x : log_pyx_right_x
            end

            # CpG group is observed in read
            call[x_index_st].obs = true

        end
        
        # Push call
        push!(calls, call)

    end
    
    # Store call info
    rs.calls = calls
rs.m = length(calls)

# Return nothing
    return nothing

end
"""
    `get_calls!(REG,CALL_FILE,CONFIG)`

    Function that stores observations in REG from calls in CALL_FILE.

    # Examples
    ```julia-repl
    julia> CpelNano.get_calls!(reg,call_file,config)
    ```
"""
function get_calls!(rs::RegStruct, call_file::String, config::CpelNanoConfig)::Nothing

    # Get calls
    if config.caller == "nanopolish"
        # nanopolish (sorted)
        data_dic = get_dic_nanopolish(rs, call_file)
        get_calls_nanopolish!(data_dic, rs, config.mod_nse)
else
        # exit if not valid caller
        print_log("Invalid option for config.caller ...")
        sleep(5)
        exit(0)
    end

    # Return nothing
    return nothing

end
"""
    `rscle_grp_mod!(RS)`

    Function that rescales group model (coarse) to per cpg model (fine).

    # Examples
    ```julia-repl
    julia> CpelNano.rscle_grp_mod!(rs)
    ```
"""
function rscle_grp_mod!(rs::RegStruct)::Nothing

    # Redefine groups as singletons
    rs.L = rs.N
    rs.dl = rs.dn
    rs.ρl = rs.ρn
    rs.Nl = fill(1, rs.N)

    # Return nothing
    return nothing

end
"""
    `pmap_analyze_reg(REG_INT,CHR,NANO,FA_REC,CONFIG)`

    Function that estimates ϕ for a given analysis region starting at REG_ST and stores MML, NME, as well as the 
    parameter vector ϕ.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_analyze_reg(reg_int,chr,nano,fa_rec,config)
    ```
"""
function pmap_analyze_reg(reg_int::UnitRange{Int64}, chr::String, nano::String, fa_rec::FASTA.Record, config::CpelNanoConfig)::RegStruct

    ## Get info about analysis region

    # Initizalize empty region structure
    rs = RegStruct()

    # Store chromosomal info (1-based)
    rs.chr = chr
    rs.chrst = minimum(reg_int)
    rs.chrend = maximum(reg_int)

    # Print out region being processed
    print_log("Processing $(rs.chr):$(rs.chrst)-$(rs.chrend)")

    # Get CpG site information (1-based)
    get_grp_info!(rs, fa_rec, config.min_grp_dist)
    length(rs.cpg_pos) > 9 || return rs
    config.verbose && print_log("CpG sites: $(rs.cpg_pos)")
    config.verbose && print_log("CpG groups: $(rs.cpg_grps)")

    # Retrieve region calls
    get_calls!(rs, nano, config)
    config.verbose && print_log("Number of reads: $(rs.m)")
    config.verbose && print_log("Average depth: $(get_ave_depth(rs))")
    config.verbose && print_log("Percentage groups observed: $(perc_gprs_obs(rs) * 100)%")
        # config.verbose && print_log("Calls: $(rs.calls)")
        # readline()

    # Skip region if not enough data
    get_ave_depth(rs) >= config.min_cov || return rs
    perc_gprs_obs(rs) >= 0.66 || return rs
    
    # Skip region if only one group
    rs.L > 1 || return rs

    ## Estimate ϕs

    # Estimate parameter vector ϕ
    get_ϕhat!(rs, config)

    # Leave if we couldn't estimate ϕ
    length(rs.ϕhat) > 0 || return rs
    config.verbose && print_log("ϕhat: $(rs.ϕhat)")
        # readline()

    ## Re-scale to per cpg site resolution
    rscle_grp_mod!(rs)
    config.verbose && print_log("rs.Nl: $(rs.Nl)")
    
    ## Statistical summaries
    get_stat_sums!(rs, config)

    # Set estimation region as processed
    rs.proc = true

    # Return
    return rs

end
"""
    `analyze_nano(NANO,FASTA,CONFIG)`

    Function that estimates ϕ for all analysis regions and stores MML, NME, as well as the parameter vector ϕ.

    # Examples
    ```julia-repl
    julia> CpelNano.analyze_nano(nano,fasta,config)
    ```
"""
function analyze_nano(nano::String, fasta::String, config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names, chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

    # Get output file names
    config.out_files = OutputFiles(config.out_dir, config.out_prefix)

    # Check if output file exists
    check_output_exists(config) && return nothing
    
    # If BED file provided, identify intervals of interest
    print_log("Checking for target regions ...")
    targ_regs = isempty(config.bed_reg) ? nothing : read_bed_reg(config)

    # Loop over chromosomes
    for i = 1:length(chr_names)
        
        # Get chromosome name and size
        chr = chr_names[i]
        chr_size = chr_sizes[i]
        
        # Print info
        print_log("Processing chr $(chr) of length $(chr_size) ...")
        
        # Get analysis region start positions
        if isnothing(targ_regs)
            # If no BED file provided partition entire chromosome
            print_log("Partitioning chr $(chr) ...")
            chr_part = get_chr_part(chr, config, fasta)
        else
            # If BED file provided, partition regions of interest. NOTE: ensure no split CpG groups
            print_log("Partitioning target regions in chr $(chr) ...")
            haskey(targ_regs, chr) || continue
            chr_part = get_chr_part(chr, config, fasta, targ_regs[chr])
        end
        
        # Get FASTA record for chr
        fa_rec = get_chr_fa_rec(chr, fasta)
        
        # Process each analysis region in chr
        out_pmap = pmap(x -> pmap_analyze_reg(x, chr, nano, fa_rec, config), chr_part)

        # Write chromosome
        write_output(out_pmap, config)

    end

    # Return nothing
    return nothing

end
"""
    `pmap_nano_smp_est(REG_INT,CHR,NANO,FA_REC,CONFIG)`

    Function that estimates computes sample statistics E[X] and E[XX] for a given analysis region.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_nano_smp_est(reg_int,chr,nano,fa_rec,config)
    ```
"""
function pmap_nano_smp_est(reg_int::UnitRange{Int64}, chr::String, nano::String, fa_rec::FASTA.Record, config::CpelNanoConfig)::RegStruct

    ## Get info about analysis region

    # Initizalize empty region structure
    rs = RegStruct()

    # Store chromosomal info (1-based)
    rs.chr = chr
    rs.chrst = minimum(reg_int)
    rs.chrend = maximum(reg_int)

    # Print out region being processed
    print_log("Processing $(rs.chr):$(rs.chrst)-$(rs.chrend)")

    # Get CpG site information (1-based)
    get_grp_info!(rs, fa_rec, config.min_grp_dist)
    length(rs.cpg_pos) > 9 || return rs
    config.verbose && print_log("CpG sites: $(rs.cpg_pos)")
    config.verbose && print_log("CpG groups: $(rs.cpg_grps)")

    # Retrieve region calls
    get_calls!(rs, nano, config)
    config.verbose && print_log("Number of reads: $(rs.m)")
    config.verbose && print_log("Average depth: $(get_ave_depth(rs))")
    config.verbose && print_log("Percentage groups observed: $(perc_gprs_obs(rs) * 100)%")
        # config.verbose && print_log("Calls: $(rs.calls)")
        # readline()

    # Skip region if not enough data
    get_ave_depth(rs) >= config.min_cov || return rs
    perc_gprs_obs(rs) >= 0.66 || return rs
    
    # Skip region if only one group
    rs.L > 1 || return rs

    ## Compute sample averages
    comp_smp_exps!(rs)
    config.verbose && print_log("E[X]: $(rs.exps.ex)")
    config.verbose && print_log("E[XX]: $(rs.exps.exx)")
        # config.verbose && readline()

    # Set estimation region as processed
    rs.proc = true

    # Return
    return rs

end
"""
    `nano_smp_est(NANO,FASTA,CONFIG)`

    Function that computes sample statistics E[X] and E[XX].

    # Examples
    ```julia-repl
    julia> CpelNano.nano_smp_est(nano,fasta,config)
    ```
"""
function nano_smp_est(nano::String, fasta::String, config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names, chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

    # Get output file names
    config.out_files = OutputFiles(config.out_dir, config.out_prefix)

    # Check if output file exists
    check_output_exists(config) && return nothing

    # If BED file provided, identify intervals of interest
    print_log("Checking for target regions ...")
    targ_regs = isempty(config.bed_reg) ? nothing : read_bed_reg(config)

    # Loop over chromosomes
    for i = 1:length(chr_names)
        
        # Get chromosome name and size
        chr = chr_names[i]
        chr_size = chr_sizes[i]
        
        # Print info
        print_log("Processing chr $(chr) of length $(chr_size) ...")
        
        # Get analysis region start positions
        print_log("Partitioning chr $(chr) ...")
        chr_part = get_chr_part(chr, config, fasta)
        
        # Get FASTA record for chr
        fa_rec = get_chr_fa_rec(chr, fasta)
        
        # Process each analysis region in chr
        out_pmap = pmap(x -> pmap_nano_smp_est(x, chr, nano, fa_rec, config), chr_part)

        # Write output
        write_output_ex(config.out_files.ex_file, out_pmap)
        write_output_exx(config.out_files.exx_file, out_pmap)
        write_output_ϕ_marg(config.out_files.theta_file, out_pmap)

    end

    # Return nothing
    return nothing

end