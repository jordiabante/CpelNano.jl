# Dependencies
using Distributed
@everywhere push!(LOAD_PATH, "./")
using CpelNano

# IO
dataDir = dirname(dirname(pathof(CpelNano)))
nano = "$(dataDir)/examples/targeted_example/nanopolish/targeted_example_noise2.0_sorted_methylation_calls.sorted.tsv"
fasta = "$(dataDir)/examples/targeted_example/reference/hg38_chr17_42879379_43290399.fa"
bed = "$(dataDir)/examples/targeted_example/regions_of_interest/brca1.bed"
outdir = "$(dataDir)/examples/targeted_example/cpelnano/"

# Output files 
mml_path = "$(outdir)/targeted_example_mml.txt"
nme_path = "$(outdir)/targeted_example_nme.txt"
theta_path = "$(outdir)/targeted_example_theta.txt"

# Remove output in outdir
isfile(mml_path) && rm(mml_path)
isfile(nme_path) && rm(nme_path)
isfile(theta_path) && rm(theta_path)

# Show region of interest
CpelNano.print_log("Regions of interest")
print(String(read(bed)))

# Configuration
min_cov = 10.0; max_size_subreg = 350; size_est_reg = 3000; max_em_init = 5; max_em_iters = 50; 
config = CpelNano.CpelNanoConfig(min_cov, max_size_subreg, size_est_reg, max_em_init, max_em_iters);
config.out_dir = outdir; config.out_prefix = "targeted_example"; config.bed_reg = bed; config.verbose = false;

# Analyze each region in chr
CpelNano.analyze_nano(nano,fasta,config)

# Print files
println("*******************************MML*******************************")
println(String(read(mml_path)))
println("*****************************************************************")
readline()
println("*******************************NME*******************************")
println(String(read(nme_path)))
println("*****************************************************************")
readline()
