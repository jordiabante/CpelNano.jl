# Dependencies
push!(LOAD_PATH,"./")
using CpelNano

# IO
dataDir = dirname(dirname(pathof(CpelNano)))
nano = "$(dataDir)/examples/full_example/nanopolish/full_example_noise2.0_methylation_calls.sorted.tsv"
fasta = "$(dataDir)/examples/full_example/reference/hg38_chr17_43023997_43145780.fa"
outdir = "$(dataDir)/examples/full_example/cpelnano/"

# Output files 
mml_path = "$(outdir)/full_example_mml.cpelnano"
nme_path = "$(outdir)/full_example_nme.cpelnano"
theta_path = "$(outdir)/full_example_theta.cpelnano"

# Remove output in outdir
isfile(mml_path) && rm(mml_path)
isfile(nme_path) && rm(nme_path)
isfile(theta_path) && rm(theta_path)

# Configuration
min_cov=2.5; max_size_subreg=500; size_an_reg=3000; max_em_init=1; max_em_iters=20; 
config = CpelNano.CpelNanoConfig(min_cov,max_size_subreg,size_an_reg,max_em_init,max_em_iters);
config.out_dir = outdir; config.out_prefix = "full_example"; config.verbose = false;

# Analyze each region in chr
CpelNano.analyze_nano(nano,fasta,config)

# Print files
# println("*******************************MML*******************************")
# println(String(read(mml_path)))
# println("*****************************************************************")
# readline()
# println("*******************************NME*******************************")
# println(String(read(nme_path)))
# println("*****************************************************************")
# readline()
println("*******************************THETA******************************")
println(String(read(theta_path)))
println("*****************************************************************")
readline()
