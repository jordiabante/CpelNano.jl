# TEMPLATE STEP 1
# 1. Create an environment for CpelNano in your home directory (only once)
# 2. Install CpelNano in the environment using (only once)
#   add https://github.com/jordiabante/CpelNano.jl.git#NatureGeneticsModel
# 3. Copy the content of this file in a file in your home directory
# 4. Replace the environment path for the correct one where the arrow points
# 5. Define the proper paths

# Load dependencies
using Distributed
@everywhere using Pkg
@everywhere Pkg.activate("$HOME/julia_env/CpelNano/") # <= create environment
@everywhere using CpelNano

# IO
outdir = ""
dataDir = ""
bam = "$(dataDir)/"
fasta = "$(dataDir)/"

# IO
config = CpelNano.CpelNanoConfig();
config.out_prefix = "prefix"
config.max_em_iters = 100
config.out_dir = outdir
config.max_em_init = 10
config.min_cov = 5.0
config.pe = true

# Make true if you want to activate informME mode
config.informme_mode = false

# Analyze BAM file
CpelNano.analyze_bam(bam,fasta,config)
