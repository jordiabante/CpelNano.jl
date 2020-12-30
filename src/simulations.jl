##################################################################################################
## Read length
##################################################################################################
"""

    `fastq_len_store(FASTQ)`
    
    Store length of reads in FASTQ file.
    
    # Examples
    ```julia-repl
    julia> CpelNano.fastq_len_store(fastq)
    ```
"""
function fastq_len_store(fastq::String,outdir::String)::Nothing

    # Get prefix
    pref = split(basename(fastq),".fastq")[1]

    # Loop over records in FASTQ file
    lens = Vector{Int64}()
    open(FASTQ.Reader,fastq) do reader
        for record in reader
            # Get lengths
            push!(lens,length(FASTQ.sequence(record)))
        end
    end

    # Write lengths
    print_log("Storing lengths...")
    open("$(outdir)/$(pref)_read_lens.txt","w") do io
        writedlm(io,lens)
    end

    # Return nothing
    return nothing

end
"""
    `gen_read_len(NREADS)`
    
    Returns a vector of read lengths following an exponential distribution obtained performing 
    ML estimation to read length of dataset from:
    
        https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md
    
    # Examples
    ```julia-repl
    julia> CpelNano.gen_read_len(10)
    ```
"""
function gen_read_len(n_reads::Int64)::Vector{Int64}
    
    # Return read lengths
    return Int.(round.(rand(Exponential(8484.84),n_reads)))

end
##################################################################################################
## Generate Nucleotide Sequences
##################################################################################################
"""

    `gen_seq(seq_len,n)`
    
    Generates a random sequence of length `seq_len` with `n` CpG sites.
    
    # Examples
    ```julia-repl
    julia> CpelNano.gen_seq(10,2)
    "CCGGGGCGGA"
    ```
"""
function gen_seq(seq_len::Int64,n::Int64)::String

    # Create random sequence
    seq = randstring("ACGT",seq_len)

    # Get rid of coincidental CpGs
    seq = replace(seq,r"CG" => (c) -> "C" * randstring("AT",1))

    # Determine spacing required
    spacing = floor(Int64,seq_len/n)
    spacing>=2 || return ""

    # Add equally spaced CG's
    seq = split(seq,"")
    ind = floor(Int64,spacing/2):spacing:seq_len
    for i in ind
        seq[i:(i+1)] = ["C","G"]
    end

    # Return sequence with n CpG sites
    return join(seq)

end
"""

    `gen_read(seq,read_len)`
    
    Returns a random read of length `read_len` from sequence `seq`.
    
    # Examples
    ```julia-repl
    julia> CpelNano.gen_read("CCGGGGCGGA",5)
    (5,"GGCGG")
    ```
"""
function gen_read(chr::String,fasta::String,read_len::Int64)::Tuple{Int64,String}
    
    # Get position
    fa_record = CpelNano.get_chr_fa_rec(chr,fasta)
    chrstart = fa_record.sequence[1]
    chrend = fa_record.sequence[length(fa_record.sequence)];
    pos = rand(chrstart:(chrend-read_len+1))

    # Extract read
    seq = convert(String,FASTA.sequence(fa_record,pos:(pos+read_len-1)));
    
    # Return tuple
    return pos,seq

end
"""

    `count_cpg(seq)`
    
    Count the number of CGs in sequence `seq`.
    
    # Examples
    ```julia-repl
    julia> CpelNano.count_cpg("CCGGGGCGGA")
    2
    ```
"""
function count_cpg(seq::String)::Int64
    ind = minimum.(findall(r"CG",seq))
    return length(ind)
end
"""

    `count_mpg(seq)`
    
    Count the number of MGs in sequence `seq`.
    
    # Examples
    ```julia-repl
    julia> CpelNano.count_mpg("CMGGGGCGGA")
    1
    ```
"""
function count_mpg(seq::String)::Int64
    ind = minimum.(findall(r"MG",seq))
    return length(ind)
end
"""

    `meth_seq(seq,methyl_vec)`
    
    Converts a CG in `seq` to MG based on its methylation status from the corresponding index in `methyl_vec`.
    
    # Examples
    ```julia-repl
    julia> CpelNano.meth_seq("CCGGGGCGGA",CpelNano.metro_hast_inst(2,10000,0.5,0.5))
    "CMGGGGMGGA"
    ```
"""
function meth_seq(seq::String,methyl_vec::Vector{Int8})::String
    
    # Get number of CpG sites and spacing
    n = length(methyl_vec)
    spacing = floor(Int64,length(seq)/n)

    # Find indices of C's in CG context
    ind = minimum.(findall(r"CG",seq))
    length(ind)==n || return ""

    # Modify sequence
    seq = split(seq,"")
    seq[ind[methyl_vec.==Int8(1)]] .= "M"
    
    # Return modified sequence
    return join(seq)

end
"""
    `aug_mod_file(reg_file,chr_name,fasta)`
    
    Augments model file. If analysis region is in `reg_file`,uses parameter vector from `reg_file`. Otherwise, 
    randomly samples parameter vector from analysis regions in `reg_file`. This function produces a text file
    with suffix `.aug`.
    
    # Examples
    ```julia-repl
    julia> CpelNano.aug_mod_file(reg_file,chr_name,fasta)
    ```
"""
function aug_mod_file(reg_file::String,chr_name::String,fasta::String)::Nothing

    # Read in model file
    reg_lines = readdlm(reg_file)

    # Get FASTA record for chr
    fa_record = get_chr_fa_rec(chr_name,fasta)

    # Initialize vectors to store parameter vector and partitions from model file
    model_part = Vector{UnitRange{Int64}}()
    ϕhat_1 = Vector{Float64}()
    ϕhat_2 = Vector{Float64}()
    ϕhat_3 = Vector{Float64}()

    # Read in parameter vectors and partitions from model files
    for i in 1:size(reg_lines)[1]
        chrst = reg_lines[i,2]
        chrend = reg_lines[i,3]
        ϕhat = parse.(Float64,split(reg_lines[i,4],','));
        push!(model_part,chrst:chrend)
        push!(ϕhat_1,ϕhat[1])
        push!(ϕhat_2,ϕhat[2])
        push!(ϕhat_3,ϕhat[3])
    end

    # Obtain full set of partitions from whole chromosome
    config = CpelNano.CpelNanoConfig()
    chr_part = CpelNano.get_chr_part(chr_name,config,fasta)

    # Initialize vector to store analysis regions
    model_regs = Vector{CpelNano.RegStruct}()

    # Iterate through every partition in the chromosome
    for i in 1:length(chr_part)
        
        # Check if partition is in model file
        modelpart_index = findfirst(x->x==chr_part[i],model_part)

        # If partition is not in model file,then randomly sample parameter vector from model file
        if isnothing(modelpart_index)
            # Randomly sample parameters
            new_ϕhat_1 = ϕhat_1[rand(1:length(ϕhat_1))]
            new_ϕhat_2 = ϕhat_2[rand(1:length(ϕhat_2))]
            new_ϕhat_3 = ϕhat_2[rand(1:length(ϕhat_3))]
        else
            # Keep the estimated vectors
            new_ϕhat_1 = ϕhat_1[modelpart_index]
            new_ϕhat_2 = ϕhat_2[modelpart_index]
            new_ϕhat_3 = ϕhat_3[modelpart_index]
        end

        # Construct RegStruct object for this analysis region
        reg_data = CpelNano.RegStruct()
        reg_data.chr = chr_name
        reg_data.chrst = minimum(chr_part[i])
        reg_data.chrend = maximum(chr_part[i])
        reg_data.ϕhat = [new_ϕhat_1,new_ϕhat_2,new_ϕhat_3]
        
        # Get group info assuming singletons
        CpelNano.get_grp_info!(reg_data,fa_record,1)

        # Set processed flag as true to store
        reg_data.proc = true

        # Add RegStruct to list of RegStructs for the chromosome
        push!(model_regs,reg_data)

    end

    # Store
    write_output_ϕ("$(reg_file).aug",model_regs)

    # Return nothing
    return nothing

end
"""
    `reads_reg_chr(reg_file,chr_name,fasta,fa_record)`
    
    Returns a list of RegStructs for every partition in chromosome given `reg_file` containing all 
    analysis regions in the chromosome. WARNING: this function assumes there is an existing model
    for each single region.
    
    # Examples
    ```julia-repl
    julia> model_regs = CpelNano.read_reg_chr(reg_file,chr_name,fasta,fa_record)
    ```
"""
function read_reg_chr(reg_file::String,chr_name::String,fasta::String,fa_record::FASTA.Record)::Vector{CpelNano.RegStruct}

    # Read in model file
    reg_lines = readdlm(reg_file)

    # Initialize vectors to store parameter vector and partitions from model file
    model_part = Vector{UnitRange{Int64}}()
    ϕhat_1 = Vector{Float64}()
    ϕhat_2 = Vector{Float64}()
    ϕhat_3 = Vector{Float64}()

    # Read in parameter vectors and partitions from model files
    for i in 1:size(reg_lines)[1]
        chrst = reg_lines[i,2]
        chrend = reg_lines[i,3]
        ϕhat = parse.(Float64,split(reg_lines[i,4],','));
        push!(model_part,chrst:chrend)
        push!(ϕhat_1,ϕhat[1])
        push!(ϕhat_2,ϕhat[2])
        push!(ϕhat_3,ϕhat[3])
    end

    # Initialize vector to store analysis regions
    model_regs = Vector{CpelNano.RegStruct}()

    # Obtain full set of partitions from whole chromosome
    config = CpelNano.CpelNanoConfig();
    chr_part = CpelNano.get_chr_part(chr_name,config,fasta);

    # Iterate through every partition in the chromosome
    for i in 1:length(chr_part)

        # Get model index
        modelpart_index = findfirst(x->x==chr_part[i],model_part)

        # Get corresponding ϕ
        new_ϕhat_1 = ϕhat_1[modelpart_index]
        new_ϕhat_2 = ϕhat_2[modelpart_index]
        new_ϕhat_3 = ϕhat_3[modelpart_index]

        # Construct RegStruct object for this analysis region
        reg_data = CpelNano.RegStruct()
        reg_data.chr = chr_name
        reg_data.chrst = minimum(chr_part[i]) 
        reg_data.chrend = maximum(chr_part[i])
        reg_data.ϕhat = [new_ϕhat_1,new_ϕhat_2,new_ϕhat_3]
        
        # Get group info assuming singletons
        CpelNano.get_grp_info!(reg_data,fa_record,1)
 
        # Add RegStruct to list of RegStructs for the chromosome
        push!(model_regs,reg_data)
 
     end
 
     # Return list of RegStructs in chromosome,in same order as chromosome partitions
     return model_regs
 
 end
"""
    `linearsearch_readstart(startpos,parts)`
    
    Find index of analysis region corresponding to start of read `startpos` given list of chromosome partitions `parts`
    
    # Examples
    ```julia-repl
    julia> CpelNano.linearsearch_readstart(14996959,CpelNano.get_chr_part("chr22",CpelNano.CpelNanoConfig(),"/Users/sandeepk/Downloads/chr22.fa"))
    5000
    ```
"""
function linearsearch_readstart(startpos::Int64,parts::Vector{UnitRange{Int64}})::Int64
    for i in 1:length(parts)
        chrst = parts[i][1]
        chrend = parts[i][length(parts[i])]
        if startpos >= chrst && startpos <= chrend
            return i
        end
    end
    return -1
end
"""
    `linearsearch_readend(startpos,idx,parts)`
    
    Find index of analysis region corresponding to end of read `startpos` given list of chromosome partitions `parts`
    and index of analysis region corresponding to start of read `idx`
    
    # Examples
    ```julia-repl
    julia> CpelNano.linearsearch_readend(20999958,5000,CpelNano.get_chr_part("chr22",CpelNano.CpelNanoConfig(),"/Users/sandeepk/Downloads/chr22.fa"))
    7001
    ```
"""
function linearsearch_readend(endpos::Int64,idx::Int64,parts::Vector{UnitRange{Int64}})::Int64
    for i in idx:length(parts)
        chrst = parts[i][1]
        chrend = parts[i][length(parts[i])]
        if endpos >= chrst && endpos <= chrend
            return i
        end
    end
    return -1
end
"""
    `get_methyl_vector(rd_st,rd_end,startidx,endidx,model_regs)`
    
    Obtain methylation vector for a given read given start of read `rd_st`,end of read `rd_end`,
    index of analysis region corresponding to start of read `startidx`,index of analysis region 
    corresponding to end of read `endidx`,and list of RegStructs in chromosome `model_regs`
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_methyl_vector(rd_st,rd_end,startidx,endidx,model_regs)
    ```
"""
function get_methyl_vector(rd_st::Int64,rd_end::Int64,startidx::Int64,endidx::Int64,model_regs::Vector{RegStruct})::Tuple{Vector{Int8},Vector{Int64}}
    
    # Initialize
    x_out = Vector{Int8}()
    cpg_pos = Vector{Int64}()

    # Concatenate all αs and βs
    @inbounds for i in startidx:endidx

        # Get region struct
        rs = model_regs[i]

        # If no CpG sites, skip
        length(rs.cpg_pos)>0 || continue

        # Get parameter vectors
        α,β = CpelNano.get_αβ_from_ϕ(rs.ϕhat,rs)

        # Generate realization x
        x = length(rs.cpg_pos)==1 ? [Int8(1)] : CpelNano.gen_x_mc(α,β)
        x_out = vcat(x_out,x)

        # Add CpG sites
        cpg_pos = vcat(cpg_pos,rs.cpg_pos)

    end

    # Return nothing
    length(cpg_pos)>0 || return [],[]

    # Keep only overlapping part
    fst_ind = findfirst(x->x>=rd_st,cpg_pos)
    isnothing(fst_ind) && return [],[]
    lst_ind = findlast(x->x<=rd_end,cpg_pos)
    isnothing(lst_ind) && return [],[]
    cpg_pos = cpg_pos[fst_ind:lst_ind]
    x_out = x_out[fst_ind:lst_ind]

    # Return x & cpg positions
    return x_out,cpg_pos

end
##################################################################################################
## Markov chain sampler of fine grain model
##################################################################################################
"""
    `get_log_g_x(n,α,β)`
    
    Compute function log(g_n(x_n)).
    
    # Examples
    ```julia-repl using CpelNano
    julia> CpelNano.get_log_g_x(1,zeros(10),zeros(9))
    (256.0,256.0)
    ```
"""
function get_log_g_x(n::Int64,α::Vector{Float64},β::Vector{Float64})::NTuple{2,Float64}

    # Get size
    N = length(α)

    # If n=N return 0 (log(1))
    n==N && return 0.0,0.0

    # Define u_n for xn = ±1
    unm1 = [-(α[n+1]-β[n])/2.0; +(α[n+1]-β[n])/2.0]
    unp1 = [-(α[n+1]+β[n])/2.0; +(α[n+1]+β[n])/2.0]

    # Define u_N
    uN = [-α[end]/2.0; α[end]/2.0]

    # If longer gap...
    n==(N-1) && return log_vec_vec_mult(unm1,uN),log_vec_vec_mult(unp1,uN)

    # Define Wn
    Wnm1 = [-(α[n+1]-β[n])/2.0+β[n+1]-α[n+2]/2.0 -(α[n+1]-β[n])/2.0-β[n+1]+α[n+2]/2.0;
            +(α[n+1]-β[n])/2.0-β[n+1]-α[n+2]/2.0 +(α[n+1]-β[n])/2.0+β[n+1]+α[n+2]/2.0]
    Wnp1 = [-(α[n+1]+β[n])/2.0+β[n+1]-α[n+2]/2.0 -(α[n+1]+β[n])/2.0-β[n+1]+α[n+2]/2.0;
            +(α[n+1]+β[n])/2.0-β[n+1]-α[n+2]/2.0 +(α[n+1]+β[n])/2.0+β[n+1]+α[n+2]/2.0]

    # If longer gap...
    n==(N-2) && return log_vec_vec_mult(log_vec_mat_mult(unm1,Wnm1),uN),log_vec_vec_mult(log_vec_mat_mult(unp1,Wnp1),uN)

    # Define Wn for n'>n+1
    Wns = Vector{Array{Float64,2}}()
    @inbounds for l=(n+2):(N-1)
        W = [-α[l]/2.0+β[l]-α[l+1]/2.0 -α[l]/2.0-β[l]+α[l+1]/2.0;
             +α[l]/2.0-β[l]-α[l+1]/2.0 +α[l]/2.0+β[l]+α[l+1]/2.0]
        push!(Wns,W)
    end
    
    # Return log g_l(x_l) for x_l∈{-1,+1}
    Wns_prod = mult_log_mats(Wns)
    gl_m1 = log_vec_vec_mult(log_vec_mat_mult(log_vec_mat_mult(unm1,Wnm1),Wns_prod),uN)
    gl_p1 = log_vec_vec_mult(log_vec_mat_mult(log_vec_mat_mult(unp1,Wnp1),Wns_prod),uN)

    # return log unm1'*Wnm1*prod(Wns)*uN,log unp1'*Wnp1*prod(Wns)*uN
    return gl_m1,gl_p1

end
"""
    `gen_x_mc(α,β)`
    
    Generate an realization of X=[X1,…,XN] using α_n and β_n,or X̄=[X̄1,…,X̄L] if using α_l and β_l. 

    # Examples
    ```julia-repl
    julia> CpelNano.gen_x_mc(ones(2),ones(1))
    2-element Array{Int8,1}:
     1
     1
    ```
"""
function gen_x_mc(α::Vector{Float64},β::Vector{Float64})::Vector{Int8}

    # Size of x
    N = length(α)

    # Initialize x
    x = zeros(Int8,N)

    # Sample first CpG site (y1/(y1+y2))
    g = get_log_g_x(1,α,β)
    log_y1 = g[1]-α[1]
    log_y2 = g[2]+α[1]
    c = max(log_y1,log_y2)
    log_den = c + log(exp(log_y1-c)+exp(log_y2-c))
    p = exp(log_y2 - log_den)
    x[1] = rand()<p ? Int8(1) : Int8(-1)
        # print_log("g: $(g)")
        # print_log("p: $(p)")
        # readline()

    # Loop over all CpG sites
    @inbounds for n in 2:(N-1)
        
        # Get g
        g = get_log_g_x(n,α,β)
            # g = get_g_x(n,α,β)

        # Get transition probability
        log_y1 = g[1]-α[n]-β[n-1]*x[n-1]
        log_y2 = g[2]+α[n]+β[n-1]*x[n-1]
        c = max(log_y1,log_y2)
        log_den = c + log(exp(log_y1-c)+exp(log_y2-c))
        p = exp(log_y2 - log_den)
            # p = g[2]*exp(α[n]+β[n-1]*x[n-1])
            # p /= (p + g[1]*exp(-α[n]-β[n-1]*x[n-1]))
            # print_log("p: $(p)")

        # Sample n-th CpG site
        x[n] = rand()<p ? Int8(1) : Int8(-1)

        # Print for debug
        # print_log("g: $(g)")
        # print_log("p: $(p)")
        # readline()

    end

    # Sample last CpG site
    log_y1 = -α[N]-β[N-1]*x[N-1]
    log_y2 = α[N]+β[N-1]*x[N-1]
    c = max(log_y1,log_y2)
    log_den = c + log(exp(log_y1-c)+exp(log_y2-c))
    p = exp(log_y2 - log_den)
    x[N] = rand()<p ? Int8(1) : Int8(-1)
        # p = exp(α[N]+β[N-1]*x[N-1])
        # p /= (p + exp(-α[N]-β[N-1]*x[N-1]))
        # x[N] = rand()<p ? Int8(1) : Int8(-1)
    
    # Return x in final state
    return x

end
###################################################################################################################
# Simulate 1. X̄~q(x̄;ϕ) 2. Y|X̄ ~ 𝒩(μ(X̄),σ(X̄))
###################################################################################################################
"""
    `cpel_samp_ont(m,Nl,αs,βs,pobs,μ_m,σ_m,μ_u,σ_u,𝒜s,ℬs)`
    
    Simulate entire model where observable is distributed according to a Gaussian distribution. The same
    Gaussian distributions are assumed regardless of the k-mer. The mean is given in terms of mA. Here we 
    assume that we are passed 
        
        α_l = (a + b⋅ρ̄_l)⋅N_l ,β_l = c/d_l

    where ρ̄_l and N_l is the average density and the number of CpG sites in the l-th group,respectively.
    
    # Examples
    ```julia-repl
    julia> L = 10; M = 25; αs = fill(0.0,L); βs = fill(0.0,L-1);
    julia> pobs = 1.0; μ_m = 60.0; σ_m = 1.0; μ_u = 70.0; σ_u = 1.0;
    julia> data = CpelNano.cpel_samp_ont(M,αs,βs,pobs,μ_m,σ_m,μ_u,σ_u);
    ```
"""
function cpel_samp_ont(m::Int64,α::Vector{Float64},β::Vector{Float64},pobs::Float64,μ_m::Float64,σ_m::Float64,
    μ_u::Float64,σ_u::Float64)::RegStruct
    
    # Distributions p(y|x̄)
    p_y_x_m = Normal(μ_m,σ_m)
    p_y_x_u = Normal(μ_u,σ_u)

    # Generate calls
    L = length(α)
    calls = Vector{Vector{MethCallCpgGrp}}()
    while length(calls)<m
        # Init call vector
        call = [MethCallCpgGrp() for l=1:L]
        # Get x̄ realization
        xcpel = gen_x_mc(α,β)
        # Fill call vector
        for l in 1:L
            # Check if group is observed
            rand()<pobs || continue
            # Methylation call (convert to bool)
            xbar = xcpel[l]==1
            # Sample normal observable
            y = xbar ? rand(p_y_x_m) : rand(p_y_x_u)
            # Compute log p(y_l|x̄_l)
            log_pyx_u = log(pdf(p_y_x_u,y))
            log_pyx_m = log(pdf(p_y_x_m,y))
            # Set call for x̄_l
            call[l].obs = true
            call[l].log_pyx_u = log_pyx_u
            call[l].log_pyx_m = log_pyx_m
        end
        push!(calls,call)
    end

    # Create simulation structure
    sim_struct = CpelNano.RegStruct(calls)

    # Return data object
    return sim_struct

end
"""
    `cpel_samp_wgbs(m,Nl,αs,βs,pobs,μ_m,σ_m,μ_u,σ_u,𝒜s,ℬs)`
    
    Simulate entire model where observable is distributed according to a Gaussian distribution. The same
    Gaussian distributions are assumed regardless of the k-mer. The mean is given in terms of mA. Here we 
    assume that we are passed 
        
        α_l = (a + b⋅ρ̄_l)⋅N_l ,β_l = c/d_l

    where ρ̄_l and N_l is the average density and the number of CpG sites in the l-th group,respectively.
    
    # Examples
    ```julia-repl
    julia> L = 20; M = 20; αs = fill(0.0,L); βs = fill(0.0,L-1); pobs = 0.5;
    julia> data = CpelNano.cpel_samp_wgbs(M,αs,βs,pobs);
    ```
"""
function cpel_samp_wgbs(m::Int64,α::Vector{Float64},β::Vector{Float64},pobs::Float64)::RegStruct

    # Generate calls
    L = length(α)
    calls = Vector{Vector{MethCallCpgGrp}}()
    while length(calls)<m
        # Init call vector
        call = [MethCallCpgGrp() for l=1:L]
        # Init obs counter
        obs_count = 15
        # Get x̄ realization
        xcpel = gen_x_mc(α,β)
        # Fill call vector
        for l in 1:L
            # Check if group is observed
            if rand()<pobs && obs_count>=15
                obs_count = 0
            end
            # Increase obs counter
            obs_count += 1
            obs_count<10 || continue
            # Methylation call (convert to bool)
            xbar = xcpel[l]==1
            # Compute log p(y_l|x̄_l)
            log_pyx_u = xbar ? log_pyx_wrong_x : log_pyx_right_x
            log_pyx_m = xbar ? log_pyx_right_x : log_pyx_wrong_x
            # Set call for x̄_l
            call[l].obs = true
            call[l].log_pyx_u = log_pyx_u
            call[l].log_pyx_m = log_pyx_m
        end
        push!(calls,call)
    end

    # Create simulation structure
    sim_struct = CpelNano.RegStruct(calls)

    # Return data object
    return sim_struct

end
##################################################################################################
## Metropolis Hastings (NOT USED)
##################################################################################################
"""
    `init_lattice(lat_length)`
    
     Randomly initializes vector of length `lat_length` with 1 and -1 representing methylated
     and unmethylated CGs respectively.
    
    # Examples
    ```julia-repl
    julia> CpelNano.init_lattice(2)
    2-element Array{Int8,1}:
     -1
     -1
    ```
"""
function init_lattice(lat_length::Int64)::Vector{Int8}

    # Return initialized lattice
    return rand([Int8(-1),Int8(1)],lat_length)

end
"""
    `get_ΔU(lattice,i,α,β)`
    
     Computes the change in potential energy U(X) by flipping the methylation state of the `i`th CG 
     in a lattice with vectors `[α1,α2,…,αN]` and `[β1,β2,…,β_{N-1}]`.
    
    # Examples
    ```julia-repl
    julia> CpelNano.get_ΔU([Int8(1),Int8(1)],2,[0.0,0.0],[0.0])
    0.0
    ```
"""
function get_ΔU(lattice::Vector{Int8},i::Int64,α::Vector{Float64},β::Vector{Float64})::Float64
    
    # Modified lattice
    mod_lattice = copy(lattice)
    mod_lattice[i] = -mod_lattice[i]

    # Original contribution of i-th CpG site
    U = -α[i]*lattice[i]
    if i==1 
        U -= β[1]*lattice[i]*lattice[i+1]
    elseif i==length(lattice)
        U -= β[end]*lattice[i]*lattice[i-1]
    else
        U -= β[i-1]*lattice[i-1]*lattice[i] + β[i]*lattice[i]*lattice[i+1]
    end

    # Original contribution of i-th CpG site
    Uflip = -α[i]*mod_lattice[i]
    if i==1 
        Uflip -= β[1]*mod_lattice[i]*mod_lattice[i+1]
    elseif i==length(mod_lattice)
        Uflip -= β[end]*mod_lattice[i]*mod_lattice[i-1]
    else
        Uflip -= β[i-1]*mod_lattice[i-1]*mod_lattice[i] + β[i]*mod_lattice[i]*mod_lattice[i+1]
    end

    # Note: if ΔU>0 → new config is less likely under θ 
    # Return ΔU
    return Uflip-U

end
"""
    `metro_hast_inst(lat_length,n_iters,α,β)`
    
    Run MCMC for `n_iters` iterations on a lattice of length `lat_length` with vectors `[α1,α2,…,αN]` 
    and `[β1,β2,…,β_{N-1}]`.
    
    # Examples
    ```julia-repl
    julia> CpelNano.metro_hast_inst(2,1000,[1.0,1.0],[1.0])
    2-element Array{Int8,1}:
     1
     1
    ```
"""
function metro_hast_inst(lat_length::Int64,n_iters::Int64,α::Vector{Float64},β::Vector{Float64})::Vector{Int8}

    # Initialize lattice
    lattice = init_lattice(lat_length)

    # Iterate til convergence of MC chain
    @inbounds for i in 1:n_iters
        @inbounds for j in 1:lat_length
            # Get energy Hamiltonian differential
            ΔU = get_ΔU(lattice,j,α,β)
            # If flipping leads to lower energy,flip
            if ΔU < 0
                lattice[j] = -lattice[j]
            else
                # If flipping leads to ΔU>0,flip based on weighted likelihood from Boltzmann distribution
                if rand() < exp(-ΔU)
                    lattice[j] = -lattice[j]
                end
            end
        end
    end

    # Return lattice in final state
    return lattice

end
