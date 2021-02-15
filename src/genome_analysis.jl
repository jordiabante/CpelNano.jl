###################################################################################################################
# GENOME PROPERTIES
###################################################################################################################
"""

    `get_genome_info(FASTA)`
    
    Get chromosome information.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_genome_info(fasta)
    ```
"""
function get_genome_info(fasta::String)::Tuple{Vector{String},Vector{Int64}}

    # Get info from FASTA file
    fa_reader = open(FASTA.Reader,fasta,index=fasta*".fai")
    chr_sizes = fa_reader.index.lengths
    chr_names = fa_reader.index.names
    close(fa_reader)

    # Return info
    return chr_names,chr_sizes

end
"""

    `get_chr_fa_rec(CHR,FASTA)`
    
    Get CHR chromosome FASTA record.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_chr_fa_rec(chr,fasta)
    ```
"""
function get_chr_fa_rec(chr::String,fasta::String)::FASTA.Record

    # Get info from FASTA file
    fa_reader = open(FASTA.Reader,fasta,index=fasta*".fai")
    fa_rec = fa_reader[chr]
    close(fa_reader)

    # Return info
    return fa_rec

end
"""
    `get_grps_from_cgs(CPG_POS,MIN_GRP_DIST)`

    Function that returns CpG groups.

    # Examples
    ```julia-repl
    julia> CpelNano.get_grps_from_cgs([10,25,28,50],10)
    3-element Array{CpelNano.CpgGrp,1}:
     CpelNano.CpgGrp(5:15, 10, 1:1, true)  
     CpelNano.CpgGrp(20:33, 26, 2:3, false)
     CpelNano.CpgGrp(45:55, 50, 4:4, true)
    ```
"""
function get_grps_from_cgs(cpg_pos::Vector{Int64},min_grp_dist::Int64)::Vector{CpgGrp}

    # Init tmp vars
    cpg_ind = 1:1
    grp_int = cpg_pos[1]:cpg_pos[1]
    length(cpg_pos)==1 && return [CpgGrp(grp_int,cpg_ind)]
    
    # Find groups
    i = 2
    cpg_grps = Vector{CpgGrp}()
    while i<=length(cpg_pos)

        # Check if grp is complete
        if cpg_pos[i]-maximum(grp_int)>min_grp_dist
            # Store CpG group
            push!(cpg_grps,CpgGrp(grp_int,cpg_ind))
            # Create new tmp vars
            cpg_ind = i:i
            grp_int = cpg_pos[i]:cpg_pos[i]
        else
            # Update tmp vars
            cpg_ind = minimum(cpg_ind):i
            grp_int = minimum(grp_int):cpg_pos[i]
        end
        
        # If last one, store
        i==length(cpg_pos) && push!(cpg_grps,CpgGrp(grp_int,cpg_ind))

        # Increase counter
        i += 1

    end

    # Adjust intervals for tails
    @inbounds for g=1:length(cpg_grps)

        # Get coordinates of first and last C in CG context
        grp_st = minimum(cpg_grps[g].grp_int)
        grp_end = maximum(cpg_grps[g].grp_int)

        # Update interval field correcting for tails
        cpg_grps[g].grp_int = (grp_st-5):(grp_end+5)

    end

    # Return groups
    return cpg_grps

end
"""
    `get_ρn(DNA_SEQ,CPG_POS_OFFSET)`

    Function that returns CpG density vector ρn.

    # Examples
    ```julia-repl
    julia> CpelNano.get_ρn(dna_seq,cpg_pos_offset)
    ```
"""
function get_ρn(dna_seq::String,cpg_pos_offset::Vector{Int64};width=1000)::Vector{Float64}

    # Loop over each CpG site
    ρn = Vector{Float64}()
    @inbounds for n=1:length(cpg_pos_offset)
        
        # Get window st and end
        win_st = max(1,cpg_pos_offset[n]-width÷2)
        win_end = min(length(dna_seq),cpg_pos_offset[n]+width÷2)

        # Get num CpG sites
        push!(ρn,length(collect(eachmatch(r"[Cc][Gg]",dna_seq[win_st:win_end]))))

    end

    # Return density
    return ρn./width

end
"""
    `get_ρl(ρn,CPG_GRPS)`

    Function that returns average CpG density per group vector ρl.

    # Examples
    ```julia-repl
    julia> CpelNano.get_ρl(ρl,cpg_grps)
    ```
"""
function get_ρl(ρn::Vector{Float64},cpg_grps::Vector{CpgGrp})::Vector{Float64}

    # Loop over each CpG group
    ρl = Vector{Float64}()
    @inbounds for l=1:length(cpg_grps)

        # Get average 
        aux = sum(ρn[cpg_grps[l].cpg_ind])
        aux /= length(cpg_grps[l].cpg_ind)
        
        # Push
        push!(ρl,aux)

    end

    # Return average group density
    return ρl

end
"""
    `get_dn(CPG_POS)`

    Function that returns distance between contiguous CpG sites.

    # Examples
    ```julia-repl
    julia> CpelNano.get_dn([10,20,30])
    ```
"""
function get_dn(cpg_pos::Vector{Int64})::Vector{Int64}

    # Loop
    dn = Vector{Int64}()
    @inbounds for n=2:length(cpg_pos)
        push!(dn,cpg_pos[n]-cpg_pos[n-1])
    end

    # Return dn
    return dn

end
"""
    `get_dl(CPG_POS,CPG_GRPS)`

    Function that returns distance between contiguous CpG groups.

    # Examples
    ```julia-repl
    julia> CpelNano.get_dl(cpg_pos,cpg_grps)
    ```
"""
function get_dl(cpg_pos::Vector{Int64},cpg_grps::Vector{CpgGrp})::Vector{Int64}

    # Loop
    dl = Vector{Int64}()
    @inbounds for l=2:length(cpg_grps)
        cpg1 = cpg_pos[maximum(cpg_grps[l-1].cpg_ind)]
        cpg2 = cpg_pos[minimum(cpg_grps[l].cpg_ind)]
        push!(dl,cpg2-cpg1)
    end

    # Return dl
    return dl

end
"""
    `get_grp_info!(REG,FA_REC,MIN_GRP_DIST)`

    Function that stores information about each CpG group, separated by a minimum of min_grp_dist, in REG.

    # Examples
    ```julia-repl
    julia> CpelNano.get_grp_info!(reg,fa_rec,min_grp_dist)
    ```
"""
function get_grp_info!(rs::RegStruct,fa_rec::FASTA.Record,min_grp_dist::Int64)::Nothing

    # Get size of chromosome
    chr_size = length(FASTA.sequence(fa_rec))
    
    # Get position CpG sites in analysis region
    dna_seq = FASTA.sequence(String,fa_rec,rs.chrst:rs.chrend)
    cpg_pos = map(x->getfield(x,:offset),eachmatch(r"[Cc][Gg]",dna_seq)) .+ rs.chrst .- 1
        # print_log("$(cpg_pos)")
    
    # Set individual CpG site information
    rs.cpg_pos = cpg_pos
    rs.N = length(cpg_pos)
    rs.N>1 || return nothing

    # Set group information
    cpg_grps = get_grps_from_cgs(cpg_pos,min_grp_dist)
    rs.Nl = map(grp->length(grp.cpg_ind),cpg_grps)
    rs.L = length(cpg_grps)
    rs.cpg_grps = cpg_grps
        # print_log("$(cpg_grps)")

    ## Density
    dna_cont_st = max(1,rs.chrst-500)
    dna_cont_end = min(chr_size,rs.chrend+500)    
    dna_seq = FASTA.sequence(String,fa_rec,dna_cont_st:dna_cont_end)
        # print_log("$(dna_seq)")

    # Store CpG density
    rs.ρn = get_ρn(dna_seq,cpg_pos.-dna_cont_st)
        # print_log("ρn=$(rs.ρn)")

    # Store average group average density
    rs.ρl = get_ρl(rs.ρn,cpg_grps)
        # print_log("ρl=$(rs.ρl)")

    ## Distances

    # Store CpG distances
    rs.dn = get_dn(cpg_pos)
        # print_log("dn=$(rs.dn)")

    # Store group distances
    rs.dl = get_dl(cpg_pos,cpg_grps)
        # print_log("dl=$(rs.dl)")

        # readline()

    # Return nothing
    return nothing

end
"""
    `get_informME_chr_part(CHR_SIZE,CONFIG)`

    Function that partitions chromosome the same way informME does.

    # Examples
    ```julia-repl
    julia> CpelNano.get_informME_chr_part(chr_size,config)
    ```
"""
function get_informME_chr_part(chr_size::Int64,config::CpelNanoConfig)::Vector{UnitRange{Int64}}

    # Loop until entire chromosome is covered
    ok = false
    ar_st = 1
    ar_end = ar_st + config.size_est_reg - 1
    chr_part = Vector{UnitRange{Int64}}()
    while !ok

        # Check end
        ar_end = ar_end>chr_size ? chr_size : ar_end
        ok = ar_end==chr_size ? true : false

        # Push region
        push!(chr_part,ar_st:ar_end)

        # Update window
        ar_st += config.size_est_reg
        ar_end = ar_st + config.size_est_reg - 1

    end

    # Return partition
    return chr_part

end
"""
    `get_chr_part(CHR,CONFIG,FASTA)`

    Function that partitions chromosome into analysis regions. No CpG group is allowed to be split in the boundary.
    There are two methods for this function: (i) entire chromosome is partitioned, and (ii) only targeted regions
    as defined by the `targ_regs` arguments are partitioned. 

    # Examples
    ```julia-repl
    julia> CpelNano.get_chr_part(chr,config,fasta)
    julia> CpelNano.get_chr_part(chr,config,fasta,targ_regs)
    ```
"""
function get_chr_part(chr::String,config::CpelNanoConfig,fasta::String)::Vector{UnitRange{Int64}}

    # Get chromosome sequence
    chr_seq = get_chr_fa_rec(chr,fasta)
    chr_size = length(FASTA.sequence(chr_seq))
    # config.verbose && print_log("chr_seq: $(chr_seq)")

    # Loop over chromosome
    i = 1
    chr_part = Vector{UnitRange{Int64}}()
    while i<chr_size

        # Candidate boundary coordinates
        bndy = i + config.size_est_reg - 1
        # config.verbose && print_log("bndy: $(bndy)")

        # Leave if boundary is beyond size
        if bndy>=chr_size 
            push!(chr_part,i:chr_size)
            break
        end
        
        # Gather DNA sequence from ±1000
        dna_seq = FASTA.sequence(String,chr_seq,(bndy-1000):min(chr_size,bndy+1000))
        # config.verbose && print_log("dna_seq: $(dna_seq[1000:1002])")

        # Obtain absolute position of CpG sites
        cpg_pos = eachmatch(r"[Cc][Gg]",dna_seq)
        cpg_pos = map(x->getfield(x,:offset),cpg_pos) .+ bndy .- 1001
        # config.verbose && print_log("cpg_pos: $(cpg_pos)")

        # Check boundary if there are CpG sites
        if length(cpg_pos)>0
        
            # Obtain groups in window
            cpg_grps = get_grps_from_cgs(cpg_pos,config.min_grp_dist)
            # config.verbose && print_log("cpg_grps: $(cpg_grps)")
            
            # Identify problematic group (if any)
            ind_grp = findfirst(grp->bndy in grp.grp_int,cpg_grps)
            # config.verbose && print_log("ind_grp: $(ind_grp)")

            # Define bndy depending on case
            if !isnothing(ind_grp)
                # Set boundary depending on extend included
                grp_st = minimum(cpg_grps[ind_grp].grp_int)
                grp_end = maximum(cpg_grps[ind_grp].grp_int)
                bndy = 2*bndy-grp_st>grp_end ? grp_end+2 : grp_st-1
            end
            # config.verbose && print_log("Modified bndy: $(bndy)")
            # config.verbose && readline()

        end

        # Store and move on to next boundary
        push!(chr_part,i:bndy)
        i = bndy+1

    end

    # Return partition
    return chr_part

end
function get_chr_part(chr::String,config::CpelNanoConfig,fasta::String,targ_regs::Vector{UnitRange{Int64}})::Vector{UnitRange{Int64}}

    # Get chromosome sequence
    chr_seq = get_chr_fa_rec(chr,fasta)
    chr_size = length(FASTA.sequence(chr_seq))

    # Loop over chromosome
    chr_part = Vector{UnitRange{Int64}}()
    for reg in targ_regs
        
        # Get start and end of target region
        i = minimum(reg)
        reg_end = maximum(reg)

        # Partition target region
        while i<reg_end

            # Candidate boundary coordinates
            bndy = i + config.size_est_reg - 1

            # Leave if boundary is beyond size
            if bndy >= reg_end 
                push!(chr_part,i:reg_end)
                break
            end

            # Gather DNA sequence from ±1000
            dna_seq = FASTA.sequence(String,chr_seq,(bndy-1000):min(chr_size,bndy+1000))
                
            # Obtain absolute position of CpG sites
            cpg_pos = eachmatch(r"[Cc][Gg]",dna_seq)
            cpg_pos = map(x->getfield(x,:offset),cpg_pos) .+ bndy .- 1001
            # config.verbose && print_log("cpg_pos: $(cpg_pos)")

            # Check boundary if there are CpG sites
            if length(cpg_pos)>0
            
                # Obtain groups in window
                cpg_grps = get_grps_from_cgs(cpg_pos,config.min_grp_dist)
                # config.verbose && print_log("cpg_grps: $(cpg_grps)")
                
                # Identify problematic group (if any)
                ind_grp = findfirst(grp->bndy in grp.grp_int,cpg_grps)
                # config.verbose && print_log("ind_grp: $(ind_grp)")

                # Define bndy depending on case
                if !isnothing(ind_grp)
                    # Set boundary depending on extend included
                    grp_st = minimum(cpg_grps[ind_grp].grp_int)
                    grp_end = maximum(cpg_grps[ind_grp].grp_int)
                    bndy = 2*bndy-grp_st>grp_end ? grp_end+2 : grp_st-1
                end
                # config.verbose && print_log("Modified bndy: $(bndy)")
                # config.verbose && readline()
                
            end

            # Store and move on to next boundary
            push!(chr_part,i:bndy)
            i = bndy+1

        end
    end

    # Return partition
    return chr_part

end
###################################################################################################################
# ANALYSIS REGION PROPERTIES
###################################################################################################################
"""
    `get_nls_reg_info!(REG_STRUCT,CONFIG)`

    Function that stores subregion information.

    # Examples
    ```julia-repl
    julia> CpelNano.get_nls_reg_info!(rs,config)
    ```
"""
function get_nls_reg_info!(rs::RegStruct,config::CpelNanoConfig)::Nothing

    # Get size of estimation region
    reg_size = (rs.chrend-rs.chrst)+1

    # Find number of necessary subregions
    k = Int(ceil(reg_size/config.max_size_subreg))

    # Base size of subregions
    base_size = Int(floor(reg_size / k))

    # Set starts
    bnds = fill(base_size,k)

    # Reminder
    rem_size = mod(reg_size,k)

    # Distribute reminder in order
    i = 1
    while rem_size > 0
        bnds[i] += 1
        rem_size -=1
        i = i==k  ? 1 : i+1
    end

    # Add up sizes to obtain offset from start
    bnds = [sum(bnds[1:i]) for i=1:k] 
    
    # Offset all by chromosome start
    bnds .+= rs.chrst .- 1
    
    # Add chromosome start
    pushfirst!(bnds,rs.chrst)
    
    # Construct analysis region intervals
    sub_int = [bnds[i]:bnds[i+1] for i=1:k]

    # Divide CpG sites indices
    cpg_inds = fill(0:0,k)
    for i=1:k
        
        # Find first and last
        fst_ind = findfirst(x->bnds[i]<=x<bnds[i+1],rs.cpg_pos)
        fst_ind = isnothing(fst_ind) ? 0 : fst_ind
        lst_ind = findlast(x->bnds[i]<=x<bnds[i+1],rs.cpg_pos)
        lst_ind = isnothing(lst_ind) ? 0 : lst_ind

        # Store info
        if i<k
            # Excluding the one on right boundary
            cpg_inds[i] = fst_ind:lst_ind
        else
            # Include the CpG on right boundary if exists
            cpg_inds[i] = fst_ind:rs.N
        end

    end
        # print_log("$(rs.cpg_pos)")

    ## Store in struct
    rs.nls_rgs = AnalysisRegions(sub_int,cpg_inds)
        # print_log("$(rs.nls_rgs)")

    # Return nothing
    return nothing

end
"""
    `pmap_subregion_table(REG_INT,FA_REC)`

    Function that returns subregions size and number of CpG sites.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_subregion_table(reg_int,fa_rec)
    ```
"""
function pmap_subregion_table(reg_int::UnitRange{Int64},fa_rec::FASTA.Record,max_size_subreg::Int64)::Array{Int64,2}

    # Get start and end
    chrst = minimum(reg_int)
    chrend = maximum(reg_int)
        # print_log("Range: $(reg_int) ...")

    # Get size of estimation region
    reg_size = chrend-chrst+1

    # Find number of necessary subregions
    k = Int(ceil(reg_size/max_size_subreg))

    # Base size of subregions
    base_size = Int(floor(reg_size / k))

    # Set starts
    bnds = fill(base_size,k)

    # Reminder
    rem_size = mod(reg_size,k)

    # Distribute reminder in order
    i = 1
    while rem_size > 0
        bnds[i] += 1
        rem_size -=1
        i = i==k  ? 1 : i+1
    end

    # Add up sizes to obtain offset from start
    bnds = [sum(bnds[1:i]) for i=1:k] 
    
    # Offset all by chromosome start
    bnds .+= chrst .- 1
    
    # Add chromosome start
    pushfirst!(bnds,chrst)
        # print_log("Bounds: $(bnds)")
    
    # Construct subregion intervals
    sub_int = [bnds[i]:bnds[i+1] for i=1:k]
    sub_len = length.(sub_int) .- 1
        # print_log("Subregion intervals: $(sub_int)")
        # print_log("Subregion lenghts: $(sub_len)")

    # Get position CpG sites
    dna_seq = FASTA.sequence(String,fa_rec,reg_int)
    cpg_pos = map(x->getfield(x,:offset),eachmatch(r"[Cc][Gg]",dna_seq)) .+ minimum(reg_int) .- 1

    # If not enough CpG sites for CpelNano return empty
    length(cpg_pos)>9 || return hcat([],[])

    # Divide CpG sites indices
    cpg_inds = fill(0:0,k)
    for i=1:k
        
        # Find first and last
        fst_ind = findfirst(x->bnds[i]<=x<bnds[i+1],cpg_pos)
        fst_ind = isnothing(fst_ind) ? 0 : fst_ind
        lst_ind = findlast(x->bnds[i]<=x<bnds[i+1],cpg_pos)
        lst_ind = isnothing(lst_ind) ? 0 : lst_ind
        
        # Store info
        if i<k
            # Excluding the one on right boundary
            cpg_inds[i] = fst_ind:lst_ind
        else
            # Include the CpG on right boundary if exists
            cpg_inds[i] = fst_ind:length(cpg_pos)
        end

    end
    n_cpgs = length.(cpg_inds)
        # print_log("CpG sites: $(cpg_inds)")
        # print_log("Num CpG sites: $(n_cpgs)")

    # Return matrix
    return hcat(sub_len,n_cpgs)

end
"""
    `subregion_table(FASTA,CONFIG)`

    Function that stores a tabulated table with columns: 

        | Start of estimation region | End of estimation region | Subregion size | # CpG sites |

    # Examples
    ```julia-repl
    julia> CpelNano.subregion_table(fasta,config)
    ```
"""
function subregion_table(fasta::String,config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names,chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

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
        
        # Get FASTA record
        fa_reader = open(FASTA.Reader,fasta,index=fasta*".fai")
        fa_rec = fa_reader[chr]
        close(fa_reader)
        
        # Process each analysis region in chr
        print_log("Processing chr $(chr) ...")
        out_pmap = pmap(x->pmap_subregion_table(x,fa_rec,config.max_size_subreg),chr_part)

        # Produce a matrix
        out_pmap = vcat(out_pmap...)

        # Write histogram
        print_log("Writing results from chr $(chr) ...")
        open("$(config.out_dir)/subregion_table.txt","a") do io
            writedlm(io,out_pmap)
        end

    end

    # Return nothing
    return nothing

end
###################################################################################################################
# ANALYSIS REGION ANALYSIS
###################################################################################################################
"""
    `pmap_estimation_region_analysis(REG_INT,FA_REC,CONFIG)`

    Function that returns a vector with analysis start and end, number of CpG sites, and number of CpG groups.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_estimation_region_analysis(reg_int,fa_rec,config)
    ```
"""
function pmap_estimation_region_analysis(reg_int::UnitRange{Int64},fa_rec::FASTA.Record,min_grp_dist::Int64)::Vector{Int64}

    # Init
    out_vec = fill(0.0,4)
    out_vec[1] = Float64(minimum(reg_int))
    out_vec[2] = Float64(minimum(reg_int))

    # Get position CpG sites
    dna_seq = FASTA.sequence(String,fa_rec,reg_int)
    cpg_pos = map(x->getfield(x,:offset),eachmatch(r"[Cc][Gg]",dna_seq)) .+ minimum(reg_int) .- 1
    length(cpg_pos)>0 || return out_vec
    
    # Store CpG sites
    out_vec[3] = Float64(length(cpg_pos))
    
    # Get groups 
    out_vec[4] = Float64(length(get_grps_from_cgs(cpg_pos,min_grp_dist)))
    
    # Return vector
    return out_vec

end
"""
    `estimation_region_table(FASTA,CONFIG)`

    Function that stores a tabulated table with columns 

        | Analysis region size | Number of CpG sites | Number of CpG groups |

    # Examples
    ```julia-repl
    julia> CpelNano.estimation_region_table(nano,fasta,config)
    ```
"""
function estimation_region_table(fasta::String,config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names,chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

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
        
        # Get FASTA record
        fa_reader = open(FASTA.Reader,fasta,index=fasta*".fai")
        fa_rec = fa_reader[chr]
        close(fa_reader)
        
        # Process each analysis region in chr
        print_log("Analyzing chr $(chr) ...")
        out_pmap = pmap(x->pmap_estimation_region_analysis(x,fa_rec,config.min_grp_dist),chr_part)

        # Produce a matrix
        out_pmap = transpose(hcat(out_pmap...))

        # Write histogram
        print_log("Writing results from chr $(chr) ...")
        open("$(config.out_dir)/estimation_region_table.txt","a") do io
            writedlm(io,out_pmap)
        end

    end

    # Return nothing
    return nothing

end
###################################################################################################################
# GROUP LENGTH ANALYSIS
###################################################################################################################
"""
    `pmap_grp_size(CHR_ST,FASTA_REC,CONFIG)`

    Function that returns a vector with CG group sizes between coordinates REG_INT.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_grp_size(chrst,fa_rec,config)
    ```
"""
function pmap_grp_size(chrst::Int64,fa_rec::FASTA.Record,min_grp_dist::Int64)::Vector{Int64}

    # Get position CpG sites
    dna_seq = FASTA.sequence(String,fa_rec,chrst:(chrst+100000))
    cpg_pos = map(x->getfield(x,:offset),eachmatch(r"[Cc][Gg]",dna_seq)) .+ chrst .- 1
    
    # Set individual CpG site information
    length(cpg_pos)>0 || return Vector{Int64}()

    # Set group information
    cpg_grps = get_grps_from_cgs(cpg_pos,min_grp_dist)

    # Return vector of sizes
    return map(g->length(g.grp_int),cpg_grps)

end
"""
    `hist_group_length(FASTA,CONFIG)`

    Function that produces a plot with size of CpG groups in reference genome FASTA.

    # Examples
    ```julia-repl
    julia> CpelNano.hist_group_length(nano,fasta,config)
    ```
"""
function hist_group_length(fasta::String,config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names,chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

    # Initialize histogram
    histo = fit(Histogram,[],0.5:1.0:50000.5)

    # Loop over chromosomes
    for i=1:length(chr_names)
        
        # Get chromosome name and size
        chr = chr_names[i]
        chr_size = chr_sizes[i]
        
        # Print info
        print_log("Processing chr $(chr) of length $(chr_size) ...")
        
        # Get analysis region start positions
        print_log("Partitioning chr $(chr) ...")
        chrst = 1:100000:(chr_size-100000)

        # Get FASTA record
        fa_reader = open(FASTA.Reader,fasta,index=fasta*".fai")
        fa_rec = fa_reader[chr]
        close(fa_reader)

        # Process each analysis region in chr
        print_log("Measuring groups in chr $(chr) ...")
        out_pmap = pmap(x->pmap_grp_size(x,fa_rec,config.min_grp_dist),chrst)

        # Produce a single vector
        out_pmap = vcat(out_pmap...)
        
        # Get histogram
        print_log("Updating histogram with chr $(chr) ...")
        histo_chr = fit(Histogram,out_pmap,0.5:1.0:50000.5)
        histo = merge(histo,histo_chr)

    end

    # Write histogram
    open("$(config.out_dir)/len_group_histogram.txt","w") do io
        writedlm(io,histo.weights)
    end

    # Return nothing
    return nothing

end
###################################################################################################################
# GROUP SIZE (# CpGs) ANALYSIS
###################################################################################################################
"""
    `pmap_grp_cpg_num(CHR_ST,FASTA_REC,CONFIG)`

    Function that returns a vector with CG group sizes in terms of numbers of CpG sites between coordinates REG_INT.

    # Examples
    ```julia-repl
    julia> CpelNano.pmap_grp_cpg_num(chrst,fa_rec,config)
    ```
"""
function pmap_grp_cpg_num(chrst::Int64,fa_rec::FASTA.Record,min_grp_dist::Int64)::Vector{Int64}

    # Get position CpG sites
    dna_seq = FASTA.sequence(String,fa_rec,chrst:(chrst+100000))
    cpg_pos = map(x->getfield(x,:offset),eachmatch(r"[Cc][Gg]",dna_seq)) .+ chrst .- 1
    
    # Set individual CpG site information
    length(cpg_pos)>0 || return Vector{Int64}()

    # Set group information
    cpg_grps = get_grps_from_cgs(cpg_pos,min_grp_dist)

    # Return vector of sizes
    return map(g->length(g.cpg_ind),cpg_grps)

end
"""
    `hist_group_cpg_num(FASTA,CONFIG)`

    Function that produces a histogram of number of CpG sites per group in reference genome FASTA.

    # Examples
    ```julia-repl
    julia> CpelNano.hist_group_cpg_num(nano,fasta,config)
    ```
"""
function hist_group_cpg_num(fasta::String,config::CpelNanoConfig)::Nothing

    # Find chromosomes
    chr_names,chr_sizes = get_genome_info(fasta)

    # Create output directory if not existant
    isdir(config.out_dir) || mkdir(config.out_dir)

    # Initialize histogram
    histo = fit(Histogram,[],0.5:1.0:10000.5)

    # Loop over chromosomes
    for i=1:length(chr_names)
        
        # Get chromosome name and size
        chr = chr_names[i]
        chr_size = chr_sizes[i]
        
        # Print info
        print_log("Processing chr $(chr) of length $(chr_size) ...")
        
        # Get analysis region start positions
        print_log("Partitioning chr $(chr) ...")
        chrst = 1:100000:(chr_size-100000)

        # Get FASTA record
        fa_reader = open(FASTA.Reader,fasta,index=fasta*".fai")
        fa_rec = fa_reader[chr]
        close(fa_reader)

        # Process each analysis region in chr
        print_log("Counting # CpG sites in groups in chr $(chr) ...")
        out_pmap = pmap(x->pmap_grp_cpg_num(x,fa_rec,config.min_grp_dist),chrst)

        # Produce a single vector
        out_pmap = vcat(out_pmap...)
        
        # Get histogram
        print_log("Updating histogram with chr $(chr) ...")
        histo_chr = fit(Histogram,out_pmap,0.5:1.0:10000.5)
        histo = merge(histo,histo_chr)

    end

    # Write histogram
    open("$(config.out_dir)/num_cpg_per_grp_histogram.txt","w") do io
        writedlm(io,histo.weights)
    end

    # Return nothing
    return nothing

end