"""
    `get_calls_nanopolish(CALL_DIC,REG_DATA)`

    Function that using the takes overlapping the analysis region and extracts information for each CpG 
    group.

    # Examples
    ```julia-repl
    julia> CpelNano.get_calls_nanopolish(call_dic,rs)
    ```
"""
function get_calls_nanopolish!(call_dic::Dict{String,Vector{Vector{SubString}}},rs::RegStruct)::Nothing

    # For each read in call_dic
    calls = Vector{Vector{MethCallCpgGrp}}()
    @inbounds for kk in keys(call_dic)
        # Initizalize vector of group calls for the read
        call = [MethCallCpgGrp() for l=1:length(rs.cpg_grps)]
        @inbounds for cpg_grp_dat in call_dic[kk]
            
            # Get chr_start (0-based in nanopolish)
            # Note: nanopolish does not consider the ±5 flanking bp
            cpg_grp_st = parse(Int64,cpg_grp_dat[3]) + 1 - 5
            cpg_grp_end = parse(Int64,cpg_grp_dat[4]) + 1 + 5
            
            # Find correct index of CpG group in region
            x_index_st = findfirst(x->minimum(x.grp_int)==cpg_grp_st,rs.cpg_grps)
            x_index_end = findfirst(x->maximum(x.grp_int)==cpg_grp_end,rs.cpg_grps)
            (isnothing(x_index_st) || isnothing(x_index_end)) && continue
                # print_log("x_index_st: $(x_index_st)")

            # Check if there is a discrepancy in CpG group coordinates
            if x_index_st != x_index_end
                print_log("[WARNING] CpG group coordinates do not match nanopolish output ... ")
                continue
            end
            
            # Record log p(y|x̄)
            call[x_index_st].log_pyx_m = parse(Float64,cpg_grp_dat[7])
            call[x_index_st].log_pyx_u = parse(Float64,cpg_grp_dat[8])
            
            # CpG group is observed in read
            call[x_index_st].obs = true
            
        end
        push!(calls,call)
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
function get_calls!(rs::RegStruct,call_file::String,config::CpelNanoConfig)::Nothing

    # Get calls
    if config.caller=="nanopolish"
        # nanopolish (sorted)
        data_dic = get_dic_nanopolish(rs,call_file)
        get_calls_nanopolish!(data_dic,rs)
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
    rs.Nl = fill(1,rs.N)

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
function pmap_analyze_reg(reg_int::UnitRange{Int64},chr::String,nano::String,fa_rec::FASTA.Record,config::CpelNanoConfig)::RegStruct

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
    get_grp_info!(rs,fa_rec,config.min_grp_dist)
    length(rs.cpg_pos)>9 || return rs
    config.verbose && print_log("CpG sites: $(rs.cpg_pos)")
    config.verbose && print_log("CpG groups: $(rs.cpg_grps)")

    # Retrieve region calls
    get_calls!(rs,nano,config)
    config.verbose && print_log("Number of reads: $(rs.m)")
    config.verbose && print_log("Average depth: $(get_ave_depth(rs))")
    config.verbose && print_log("Percentage groups observed: $(perc_gprs_obs(rs)*100)%")
        # config.verbose && print_log("Calls: $(rs.calls)")
        # readline()

    # Skip region if not enough data
    get_ave_depth(rs)>=config.min_cov || return rs
    perc_gprs_obs(rs)>=0.66 || return rs

    # Skip region if only one group
    rs.L>1 || return rs

    ## Estimate ϕs

    # Estimate parameter vector ϕ
    get_ϕhat!(rs,config)

    # Leave if we couldn't estimate ϕ
    length(rs.ϕhat)>0 || return rs
    config.verbose && print_log("ϕhat: $(rs.ϕhat)")
        # readline()

    ## Re-scale to per cpg site resolution
    rscle_grp_mod!(rs)
    config.verbose && print_log("rs.Nl: $(rs.Nl)")
    
    ## Estimate intermediate quantities
    
    # Get matrices for transfer matrix method with logtrick
    get_rs_lgtrck_mats!(rs)
    
    # Store partition function
    get_rs_logZ!(rs)

    # Get E[X] & E[XX]
    get_rs_exps!(rs)

    # Analysis region information
    get_nls_reg_info!(rs,config)

    # Store log(g_i(±⋅))
    get_rs_log_gs!(rs)

    # Compute expectations E[log g(X)]
    get_rs_exp_log_g1!(rs)
    get_rs_exp_log_g2!(rs)

    ## Compute output

    # Compute μ(X)
    comp_mml!(rs)
        # print_log("μ: $(rs.mml)")

    # Compute h(X)
    comp_nme!(rs)
        # print_log("h: $(rs.nme)")
        # readline()

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
function analyze_nano(nano::String,fasta::String,config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names,chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

    # Get output file names
    config.out_files = OutputFiles(config.out_dir,config.out_prefix)

    # Check if output file exists
    check_output_exists(config) && return nothing

    # If BED file provided, identify intervals of interest
    print_log("Checking for target regions ...")
    targ_regs = isempty(config.bed_reg) ? nothing : read_bed_reg(config)

    # Loop over chromosomes
    for i=1:length(chr_names)
        
        # Get chromosome name and size
        chr = chr_names[i]
        chr_size = chr_sizes[i]
        
        # Print info
        print_log("Processing chr $(chr) of length $(chr_size) ...")
        
        # Get analysis region start positions
        if isnothing(targ_regs)
            # If no BED file provided partition entire chromosome
            print_log("Partitioning chr $(chr) ...")
            chr_part = get_chr_part(chr,config,fasta)
        else
            # If BED file provided, partition regions of interest. NOTE: ensure no split CpG groups
            print_log("Partitioning target regions in chr $(chr) ...")
            haskey(targ_regs,chr) || continue
            chr_part = get_chr_part(chr,config,fasta,targ_regs[chr])
        end
        
        # Get FASTA record for chr
        fa_rec = get_chr_fa_rec(chr,fasta)
        
        # Process each analysis region in chr
        out_pmap = pmap(x->pmap_analyze_reg(x,chr,nano,fa_rec,config),chr_part)

        # Write chromosome
        write_output(out_pmap,config)

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
function pmap_nano_smp_est(reg_int::UnitRange{Int64},chr::String,nano::String,fa_rec::FASTA.Record,config::CpelNanoConfig)::RegStruct

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
    get_grp_info!(rs,fa_rec,config.min_grp_dist)
    length(rs.cpg_pos)>9 || return rs
    config.verbose && print_log("CpG sites: $(rs.cpg_pos)")
    config.verbose && print_log("CpG groups: $(rs.cpg_grps)")

    # Retrieve region calls
    get_calls!(rs,nano,config)
    config.verbose && print_log("Number of reads: $(rs.m)")
    config.verbose && print_log("Average depth: $(get_ave_depth(rs))")
    config.verbose && print_log("Percentage groups observed: $(perc_gprs_obs(rs)*100)%")
        # config.verbose && print_log("Calls: $(rs.calls)")
        # readline()

    # Skip region if not enough data
    get_ave_depth(rs)>=config.min_cov || return rs
    perc_gprs_obs(rs)>=0.66 || return rs

    # Skip region if only one group
    rs.L>1 || return rs

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
function nano_smp_est(nano::String,fasta::String,config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names,chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

    # Get output file names
    config.out_files = OutputFiles(config.out_dir,config.out_prefix)

    # Check if output file exists
    check_output_exists(config) && return nothing

    # If BED file provided, identify intervals of interest
    print_log("Checking for target regions ...")
    targ_regs = isempty(config.bed_reg) ? nothing : read_bed_reg(config)

    # Loop over chromosomes
    for i=1:length(chr_names)
        
        # Get chromosome name and size
        chr = chr_names[i]
        chr_size = chr_sizes[i]
        
        # Print info
        print_log("Processing chr $(chr) of length $(chr_size) ...")
        
        # Get analysis region start positions
        print_log("Partitioning chr $(chr) ...")
        chr_part = get_chr_part(chr,config,fasta)
        
        # Get FASTA record for chr
        fa_rec = get_chr_fa_rec(chr,fasta)
        
        # Process each analysis region in chr
        out_pmap = pmap(x->pmap_nano_smp_est(x,chr,nano,fa_rec,config),chr_part)

        # Write output
        write_output_ex(config.out_files.ex_file,out_pmap)
        write_output_exx(config.out_files.exx_file,out_pmap)

    end

    # Return nothing
    return nothing

end