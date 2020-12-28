# TEMPLATE STEP 1
# 1. Create an environment for CpelNano in your home directory (only once)
# 2. Install CpelNano in the environment using (only once)
#   add https://github.com/jordiabante/CpelNano.jl.git#TowardsContinuousY
# 3. Copy the content of this file in a file in your home directory
# 4. Replace the environment path for the correct one where the arrow points
# 5. Define the proper paths

# Load dependencies
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("/home-net/home-2/jabante1@jhu.edu/julia_env/CpelNano/") # <=
@everywhere using CpelNano

# Define paths (you can pass arguments and access them through the array ARGS, e.g. ARGS[1])
indir = ""
nano = "$(indir)/"
fasta = "$(indir)/"
outdir = ""
out_prefix = ""

# Configure CpelNano
min_cov = 5.0               # Minimum average depth
max_em_init = 10            # Max number of inits in Cpel EM
max_em_iters = 100          # Max number of iters in Cpel EM

# Configure CpelNano
config = CpelNano.CpelNanoConfig()
config.out_dir = outdir
config.min_cov = min_cov
config.out_prefix = out_prefix
config.max_em_init = max_em_init
config.max_em_iters = max_em_iters

# Analyze nanopore data run
CpelNano.analyze_nano(nano,fasta,config)
