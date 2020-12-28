#################################################################################################################
# AIMS 1-3:
# * Computes E[X] and E[XX] from model file used in aims 1-3
#################################################################################################################
using CpelNano
using LinearAlgebra

#################################################################################################################
# FUNCTIONS
#################################################################################################################
function process_model_file(chr::String,fasta::String,mod_file::String)::Nothing

    # IO
    outdir = dirname(mod_file)
    prefix = split(basename(mod_file),".")[1]
    ex_file = "$(outdir)/$(prefix)_ex.txt"
    exx_file = "$(outdir)/$(prefix)_exx.txt"
    
    # FASTA record
    fa_rec = CpelNano.get_chr_fa_rec(chr,fasta)

    # Go over all lines in file
    open(mod_file) do f
        for line in enumerate(eachline(f))

            # Get line data
            line_vec = split(line[2],"\t")

            # Get chromosome info
            reg_chr = line_vec[1]

            # Init struct
            rs = CpelNano.RegStruct()
            rs.chr = reg_chr
            
            # Get start or end
            rs.chrst = parse(Int64,line_vec[2])
            rs.chrend = parse(Int64,line_vec[3])
            CpelNano.print_log("st: $(rs.chrst) - end: $(rs.chrend)")
            
            # Get genome info
            CpelNano.get_grp_info!(rs,fa_rec,1)
            
            # Get αs & βs
            αs = parse.(Float64,String.(split(line_vec[5],",")))
            βs = parse.(Float64,String.(split(line_vec[6],",")))
            
            # Compute E[X] & E[XX]
            u1 = CpelNano.get_log_u(αs[1])
            uN = CpelNano.get_log_u(αs[end])
            Ws = [CpelNano.get_log_W(αs[n],αs[n+1],βs[n]) for n=1:length(βs)]
            Z = CpelNano.get_log_Z(u1,uN,Ws)
            ex = CpelNano.get_E_X_log(u1,uN,Ws,Z)
            exx = CpelNano.get_E_XX_log(u1,uN,Ws,Z)

            # Write E[X]
            open(ex_file,"a") do io
                for n=1:rs.N
                    write(io,"$(rs.chr)\t$(rs.cpg_pos[n])\t$(rs.cpg_pos[n])\t$(ex[n])\n")
                end
            end

            # Write E[XX]
            open(exx_file,"a") do io
                for n=1:(rs.N-1)
                    write(io,"$(rs.chr)\t$(rs.cpg_pos[n])\t$(rs.cpg_pos[n+1])\t$(exx[n])\n")
                end
            end

        end
    end
    
    # Return nothing
    return nothing

end
#################################################################################################################
# Calls
#################################################################################################################

# Model file aims 1-2
chr = "chr22"
fasta = "/Users/jordiabante/Desktop/chr22.fa"
mod_file = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/GM12878_wgbs_cpelnano_theta_aug.txt"
process_model_file(chr,fasta,mod_file)

# Model file aim 3
chr = "chr22"
fasta = "/Users/jordiabante/Desktop/chr22.fa"
mod_file = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/GM12878/chr22/GM12878_wgbs_cpelnano_theta.txt"
process_model_file(chr,fasta,mod_file)
