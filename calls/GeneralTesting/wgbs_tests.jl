#####################################################################################################
# Read BAM
#####################################################################################################
using CpelNano
using FASTX.FASTA

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/wgbs/"
fasta = "$(dataDir)/reference/ref.fa"
bam = "$(dataDir)/bam/g1_s1.bam"
chr = "chrTest"

# Get FASTA data
fa_rec = CpelNano.get_chr_fa_rec(chr,fasta)

# Get size of chromosome
dna_seq = FASTA.sequence(fa_rec)
chr_size = length(dna_seq)

# Get position CpG sites in analysis region
diff = 151
chr_st = 1
chr_end = 1000-diff
dna_seq = FASTA.sequence(String,fa_rec,chr_st:chr_end)
cpg_pos = map(x->getfield(x,:offset),eachmatch(r"[Cc][Gg]",dna_seq)) .+ chr_st .- 1

# Read bam (get calls)
calls = CpelNano.get_calls_bam(bam,chr,chr_st,chr_end,cpg_pos,false,(0,0,0,0))

#####################################################################################################
# Split BAM file
#####################################################################################################
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/wgbs/"
fasta = "$(dataDir)/reference/ref.fa"
bam = "$(dataDir)/bam/meth_sample.bam"
win_size = 500

# Split BAM
CpelNano.split_bam(bam,fasta,win_size)

#####################################################################################################
# ESTIMATION
#####################################################################################################

## Process chromosome
using CpelNano

# IO
dataDir = "$(dirname(dirname(pathof(CpelNano))))/examples/wgbs/"
fasta = "$(dataDir)/reference/ref.fa"
bam = "$(dataDir)/bam/meth_sample.bam"
outdir = "$(dataDir)/cpelnano/"

# Configuration
min_cov=2.5; max_size_subreg=500; size_est_reg=1000; max_em_init=5; max_em_iters=20; 
config = CpelNano.CpelNanoConfig(min_cov,max_size_subreg,size_est_reg,max_em_init,max_em_iters);
config.out_dir = outdir; config.out_prefix = "test";

# Analyze each region in chr
CpelNano.analyze_bam(bam,fasta,config)

# Output files 
mml_path = "$(dataDir)/cpelnano/test_mml.bedGraph"
nme_path = "$(dataDir)/cpelnano/test_nme.bedGraph"
theta_path = "$(dataDir)/cpelnano/test_theta.bedGraph"

# Print files
# println(String(read(mml_path)))
# println(String(read(nme_path)))
println(String(read(theta_path)))