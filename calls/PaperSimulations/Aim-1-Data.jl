## Dependencies
using Pkg
Pkg.activate("/ibox/afeinbe2/jabante1/CpelNano/CpelNano.jl/")     # <== Modify with your environment path
using CpelNano

## Arguments
out_dir = ARGS[1]       # output directory
out_prefix = ARGS[2]    # output prefix
fa_file = ARGS[3]       # fasta file
mod_file = ARGS[4]      # model file

## Function
function gen_seqs_sims(out_dir,out_prefix,fa_file,mod_file)
    
    # Set variables
    count = 1
    avg_cov = 0.0

    # Load fasta record chromosome 22
    CpelNano.print_log("Loading fasta record ...")
    fa_record = CpelNano.get_chr_fa_rec("chr22",fa_file)
    chrstart = fa_record.sequence[1]
    chrend = fa_record.sequence[length(fa_record.sequence)]

    # Get chromosome 22 partition (use default partition by nanopolish)
    CpelNano.print_log("Partitioning chromosome ...")
    config = CpelNano.CpelNanoConfig()
    chr_part = CpelNano.get_chr_part("chr22",config,fa_file)

    # Get models from chromosome 22
    CpelNano.print_log("Read in models ...")
    model_regs = CpelNano.read_reg_chr(mod_file,"chr22",fa_file,fa_record)

    # Loop until we reach average 25X coverage
    CpelNano.print_log("Starting read generation...")
    while(avg_cov < 25.0)  

        # Generate read_length and read
        read_len = CpelNano.gen_read_len(1)[1]  
        genom_pos,new_read = CpelNano.gen_read("chr22",fa_file,read_len)

        # Filter based on read length
        read_len>100 || continue

        # Discard read if it contains N's at all
        n_match = map(x->getfield(x,:offset),eachmatch(r"N",new_read))
        length(n_match)/read_len<0.05 || continue

        # Get indices of estimation regions
        startidx = CpelNano.linearsearch_readstart(genom_pos,chr_part)
        testidx = CpelNano.linearsearch_readend(genom_pos+read_len-1,startidx,chr_part)

        # Get Ising realization and CpG positions. Skip if no CpG sites
        x_vec,cpg_pos = CpelNano.get_methyl_vector(genom_pos,genom_pos+read_len-1,startidx,testidx,model_regs)
        length(cpg_pos)>0 || continue

        # Output methylation lattice
        io = open(string(out_dir,out_prefix,"_",count,"_x_vec"),"w")
        @inbounds for i in 1:length(x_vec)
            cpg_loc = cpg_pos[i] + genom_pos
            print(io,cpg_loc)
            print(io,'\t')
            println(io,x_vec[i])
        end
        new_methyl_read = CpelNano.meth_seq(new_read,x_vec)
        close(io)

        # Print methylated sequence to a fasta file
        io = open(string(out_dir,out_prefix,"_",count,"_reads.fa"),"w")
        println(io,string(">Read",count))
        println(io,new_methyl_read)
        close(io)

	    # Calculate new average coverage (substracting N's in chr22)
        avg_cov += read_len / (chrend - chrstart + 1 - 11658691)

        # Increase count counter
        count += 1

        # Print progress every 500 reads
	    mod(count,100)==0 &&  CpelNano.print_log("Read ID: $(count). Avg Cov: $(avg_cov)")

    end

    CpelNano.print_log("Finished read generation. Total number reads: $(count)")

    return nothing

end

# Function call
gen_seqs_sims(out_dir,out_prefix,fa_file,mod_file)
