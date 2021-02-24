"""
    `get_align_strand(PAIRED_END,FLAG1,FLAG2)`

    Function that returns the strand of methylation call based on Bismark's logic. In single-end mode,
    OT and CTOT reads will both receive a FLAG of 0 while OB and CTOB both receive a FLAG of 16. In
    paired end mode:

                    Read 1       Read 2
        OT:         99          147
        OB:         83          163
        CTOT:      147           99
        CTOB:      163           83

    This table was extracted from
        https://github.com/FelixKrueger/Bismark/issues/151

    # Examples
    ```julia-repl
    julia> CpelTdm.get_align_strand(true,UInt16(99),UInt16(147))
    "OT"
    ```
"""
function get_align_strand(pe::Bool, flag1::UInt16, flag2::UInt16)::String

    # Treat SE and PE separately
    if !pe
        if flag1 == 0
            s = "OT"
        elseif flag1 == 16
            s = "OB"
        else
            print_log("SE: was expecting FLAG 0 or and encountered $(flag1) instead.")
            print_log("Exiting julia ...")
            exit(1)
        end
    else
        if flag1 == 99 && flag2 == 147
            s = "OT"
        elseif flag1 == 83 && flag2 == 163
            s = "OB"
        elseif flag1 == 147 && flag2 == 99
            s = "CTOT"
        elseif flag1 == 163 && flag2 == 83
            s = "CTOB"
        else
            print_log("PE: unexpected flag combination. Expected 99/147, 147/99, 83/163, 163/83.")
            print_log("Exiting julia ...")
            exit(1)
        end
    end

    # Return strand
    return s

end
"""
    `order_bams(PAIRED_END,RECORDS)`

    Function that returns an AlignTemp object with R1 as the first record in the forward strand and
    with the methylation call strand taken from `get_align_strand`.

    # Examples
    ```julia-repl
    julia> CpelTdm.order_bams(true,RECORDS)
    ```
"""
function order_bams(pe::Bool, records::Vector{BAM.Record})::AlignTemp

    # Check which record appears first in forward strand
    if BAM.position(records[1]) <= BAM.position(records[2])
        s = get_align_strand(pe, BAM.flag(records[1]), BAM.flag(records[2]))
        return AlignTemp(s, records[1], records[2])
    else
        s = get_align_strand(pe, BAM.flag(records[2]), BAM.flag(records[1]))
        return AlignTemp(s, records[2], records[1])
    end

end
"""
    `clean_records(PAIRED_END,RECORDS)`

    Function that takes a set of records and returns an AllAlignTemps object that contains all the
    properly aligned templates as an array of AlignTemp, which contains information about the
    methylation call strand as well as the relevant BAM records. In the PE case, R1 corresponds to
    the BAM record that appears before in the forward strand.

    # Examples
    ```julia-repl
    julia> CpelTdm.clean_records(true,RECORDS)
    ```
"""
function clean_records(pe::Bool, records::Vector{BAM.Record})::AllAlignTemps

    # Initialize struct
    out = AllAlignTemps(pe, [])

    # Consider PE vs. SE
    if pe
        # Handle PE case, first retrieve unique template names
        temp_names = unique([BAM.tempname(x) for x in records])
        for temp_name in temp_names
            # Get records with same template name
            temp_recs = filter(x -> BAM.tempname(x) == temp_name, records)

            # There should only be two records with the same template name
            length(temp_recs) == 2 && push!(out.templates, order_bams(pe, temp_recs))
        end
    else
        # Handle SE case
        out.templates = [AlignTemp(get_align_strand(pe, BAM.flag(x), UInt16(0)), x, BAM.Record()) for x in records]
    end

    # Return struct
    return out

end
"""
    `try_olaps(READER,CHR,WINDOW)`

    Function that tries to find BAM records overlaping with `CHR` at positions `WINDOW`.

    # Examples
    ```julia-repl
    julia> CpelTdm.try_olaps(reader,chr,win)
    ```
"""
function try_olaps(reader::BAM.Reader, chr::String, win::UnitRange{Int64})::Vector{BAM.Record}

    # NOTE: known bug that has to do w/ BioAlignments when close to end?
    records_olap =
    try
        collect(eachoverlap(reader, chr, minimum(win):maximum(win)))
    catch
        Array{BAM.Record,1}()
    end

    return records_olap

end
"""
    `split_bam(BAM,FASTA,WIN_SIZE)`

    Function that splits BAM file in equally size windows of size WIN_SIZE.

    # Examples
    ```julia-repl
    julia> CpelTdm.split_bam(bam,fasta,win_size)
    ```
"""
function split_bam(bam::String, fasta::String, win_size::Int64)::Nothing

    # Output BAM files prefix
    bam_prefix = splitext(bam)[1]

    # Get FASTA info
    chr_names, chr_sizes = get_genome_info(fasta)

    # Open reader
    j = 1
    reader = open(BAM.Reader, bam, index=bam * ".bai")
    header = BAM.header(reader)

    for (i, chr) in enumerate(chr_names)

        # Print info
        print_log("Splitting chromosome $(chr)...")

        # Get partitions
        chr_part = 1:win_size:chr_sizes[i]

        # Get reads overlapping with each region
        for chr_st in chr_part

            # Get end
            chr_end = min(chr_sizes[i] - 200, chr_st + win_size)

            # Get overlap reads
            recs = try_olaps(reader, chr, chr_st:chr_end)

            # Open BAM writer with header
            bamw = BAM.Writer(BGZFStream(open(bam_prefix * "_$(j).bam", "w"), "w"), header)

            # Write records
            map(rec -> write(bamw, rec), recs)

            # Close writer
            close(bamw)

            # Increase counter
            j += 1

        end

    end

    # Close reader
    close(reader)

    # Return nothing
    return nothing

end
"""
    `get_calls_bam(BAM_PATH,CHR,FEAT_ST,FEAT_END,CPG_POS,PE,TRIM)`

    Function that reads in BAM file in `BAM_PATH` and returns methylation vectors for those records
    that overlap with (1-based) genomic coordinates `chr:chrSt-chrEnd` at `cpg_pos`. A trimming given
    by TRIM=(Rf_5p,Rf_3p,Rr_5p,Rr_3p) is applied to the reads. The information was taken from:

        XM: meth calls (https://github.com/FelixKrueger/Bismark/tree/master/Docs)
        XS: uniqueness (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

    For info on OT, CTOT, OB, CTOB nomenclature see

        http://software.broadinstitute.org/software/igv/book/export/html/37.

    # Examples
    ```julia-repl
    julia> CpelTdm.get_calls_bam(BAM_PATH,"chr1",30,80,[40,60],false,(0,0,0,0))
    ```
"""
function get_calls_bam(bam::String,chr::String,roi_st::Int64,roi_end::Int64,cpg_pos::Vector{Int64},
    pe::Bool,trim::NTuple{4,Int64})::Tuple{Int64,Vector{Vector{MethCallCpgGrp}}}

    # Init return object
    calls = Vector{Vector{MethCallCpgGrp}}()

    # Get records overlapping window.
    reader = open(BAM.Reader, bam, index=bam * ".bai")
    chr_size = reader.refseqlens[findfirst(x -> x == chr, reader.refseqnames)]
    roi_end = min(roi_end, chr_size - 151)
    records_olap = try_olaps(reader, chr, roi_st:roi_end)
    length(records_olap) > 0 || return 0, calls
    close(reader)

    # Relevant flags in BAM file (both Bismark and Arioc)
    filter!(x -> (BAM.ismapped(x)) && (haskey(x, "XM")) && (!haskey(x, "XS")) && (BAM.flag(x) in
           FLAGS_ALLOWED) && (BAM.mappingquality(x) > THRESH_MAPQ),records_olap)

    # Clean records & organize them
    records_org = clean_records(pe, records_olap)

    # Loop over organized records
    for record in records_org.templates

        # Get methylation call and offset (depends on strand where call is made)
        if !(records_org.paired)
            
            # Obtain methylation call from single end
            meth_call = record.R1[:"XM"]
            meth_call = SubString(meth_call, (1 + trim[1]):(length(meth_call) - trim[2]))
            OFFSET = record.strand in ["OT","CTOT"] ? 1 : 2
            OFFSET -= BAM.position(record.R1) + trim[1]

        else
            
            # Obtain methylation call
            R1_call = record.R1[:"XM"]
            R2_call = record.R2[:"XM"]
            R1_call = SubString(R1_call, (1 + trim[1]):(length(R1_call) - trim[2]))
            R2_call = SubString(R2_call, (1 + trim[4]):(length(R2_call) - trim[3]))
            dist_st = abs(BAM.position(record.R2) + trim[4] - BAM.position(record.R1) - trim[1])
            meth_call = R1_call * "."^max(0, dist_st - length(R1_call))
            meth_call *= R2_call[max(1, (length(R1_call) + 1 - dist_st)):end]
            OFFSET = record.strand in ["OT","CTOT"] ? 1 : 2
            OFFSET -= BAM.position(record.R1) + trim[1]

        end

        # Cross positions of CpG sites if template contains CpG sites
        obs_cpgs = [x.offset for x in eachmatch(r"[zZ]", meth_call)] .- OFFSET
        length(obs_cpgs) > 0 || continue

        # If overlapping CpG sites then store and add to xobs
        olap_cpgs = findall(x -> x in obs_cpgs, cpg_pos)
        length(olap_cpgs) > 0 || continue
        x = zeros(Int64, length(cpg_pos))
        obs_cpgs = obs_cpgs[findall(x -> x in cpg_pos, obs_cpgs)] .+ OFFSET
        x[olap_cpgs] = reduce(replace, ["Z" => 1,"z" => -1], init=split(meth_call[obs_cpgs], ""))

        # Define call for read pair
        call = Vector{MethCallCpgGrp}() 
        @inbounds for n = 1:length(x)

            # Init call for CpG site
            cpg_call = MethCallCpgGrp()

            if x[n] == -1 
                # If call is unmethylated assign p(y=-1|x=-1)≈1
                cpg_call.obs = true
                cpg_call.log_pyx_u = log_pyx_right_x
                cpg_call.log_pyx_m = log_pyx_wrong_x
            elseif x[n] == 1 
                # If call is methylated assign p(y=1|x=1)≈1
                cpg_call.obs = true
                cpg_call.log_pyx_u = log_pyx_wrong_x
                cpg_call.log_pyx_m = log_pyx_right_x
            else    
                # nothing
            end

            # Push call
            push!(call, cpg_call)

        end

        # Add to set of observations
        push!(calls, call)

    end # end loop over templates sequenced

    # Return
    return length(calls), calls

end
"""
    `pmap_analyze_reg_bam(REG_INT,CHR,NANO,FA_REC,CONFIG)`

    Function that estimates ϕ for a given analysis region starting at REG_ST and stores MML, NME, as well as the 
    parameter vector ϕ.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_analyze_reg_bam(reg_int,chr,nano,fa_rec,config)
    ```
"""
function pmap_analyze_reg_bam(reg_int::UnitRange{Int64}, chr::String, bam::String, fa_rec::FASTA.Record, config::CpelNanoConfig)::RegStruct

    ## Get info about analysis region

    # Initizalize empty region structure
    rs = RegStruct()

    # Store chromosomal info (1-based)
    rs.chr = chr
    rs.chrst = minimum(reg_int)
    rs.chrend = maximum(reg_int)

    # Print out region being processed
    print_log("Processing $(rs.chr):$(rs.chrst)-$(rs.chrend)")

    # Get CpG site positions (1-based)
    get_grp_info!(rs, fa_rec, 1)
    length(rs.cpg_pos) > 9 || return rs
    config.verbose && print_log("CpG sites: $(length(rs.cpg_pos))")
    config.verbose && print_log("CpG groups: $(length(rs.cpg_grps))")

    # Retrieve region calls
    rs.m, rs.calls = get_calls_bam(bam, chr, rs.chrst, rs.chrend, rs.cpg_pos, config.pe, config.trim)
    config.verbose && print_log("Number of reads: $(rs.m)")
    config.verbose && print_log("Average depth: $(get_ave_depth(rs))")
    config.verbose && print_log("Percentage groups observed: $(perc_gprs_obs(rs) * 100)%")
    # config.verbose && length(rs.calls)>0 && print_log("Calls: $(rs.calls[1])")
    # config.verbose && readline()

    # Skip region if not enough data
    get_ave_depth(rs) >= config.min_cov || return rs
    perc_gprs_obs(rs) >= 0.66 || return rs

    # Skip region if only one group
    rs.L > 1 || return rs

    ## Estimate ϕs

    # Estimate parameter vector ϕs
    get_ϕhat_wgbs!(rs, config)

    # Leave if we couldn't estimate ϕ
    length(rs.ϕhat) > 0 || return rs
    config.verbose && print_log("ϕhat: $(rs.ϕhat)")
        # readline()

    ## Statistical summaries
    get_stat_sums!(rs, config)

    # Set estimation region as processed
    rs.proc = true

    # Return
    return rs

end
"""
    `analyze_bam(BAM,FASTA,CONFIG)`

    Function that estimates ϕ for all analysis regions and stores MML, NME, as well as the parameter vector ϕ.

    # Examples
    ```julia-repl
    julia> CpelNano.analyze_bam(bam,fasta,config)
    ```
"""
function analyze_bam(bam::String, fasta::String, config::CpelNanoConfig)::Nothing

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
        if config.informme_mode

            # Partition the chromosome as informME does
            chr_part = get_informME_chr_part(chr_size, config)
        
        else

            # Take into account group when partitioning
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

        end
        
        # Get FASTA record for chr
        fa_rec = get_chr_fa_rec(chr, fasta)
        
        # Process each analysis region in chr
        out_pmap = pmap(x -> pmap_analyze_reg_bam(x, chr, bam, fa_rec, config), chr_part)

        # Write chromosome output
        write_output(out_pmap, config)

    end

    # Return nothing
    return nothing

end