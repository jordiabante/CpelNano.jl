module CpelNano

###################################################################################################
# DEPENDENCIES
###################################################################################################
using Dates
using Random
using NLsolve
using XAM.BAM
using Calculus
using StatsBase
using BGZFStreams
using Distributed
using FASTX.FASTA
using FASTX.FASTQ
using Combinatorics
using LinearAlgebra
using Distributions
using DelimitedFiles
using MultipleTesting
using GenomicFeatures:eachoverlap

##################################################################################################
## Constants
##################################################################################################

const LOG2 = log(2.0)

# WGBS data
const THRESH_MAPQ = 39                      # MAPQ threshold (-10*log(p)) only true uni-reads
const FLAGS_ALLOWED = [0,16,83,99,147,163]  # Flags allowed in BAM recs
const log_pyx_wrong_x = -50.0
const log_pyx_right_x = 0.0 

# Paramater space boundaries
AMAX = 20.0
BMAX = 150.0
CMAX = 50.0

# Transfer Matrix
const US = UniformScaling(1.0e-1)
const D_PX = [0.0 0.0; 0.0 1.0]
const D_EX = [-1.0 0.0; 0.0 1.0]
const D_EXX = [1.0 -1.0; -1.0 1.0]
const logD1 = log.([1.0 0.0; 0.0 3.0])
const logD2 = log.([3.0 1.0; 1.0 3.0])
const logD3 = log.([2.0 0.0; 0.0 3.0])

# Hypthesis Testing
const LMAX = 1000

###################################################################################################
# SOURCE CODE
###################################################################################################
include("structures.jl")
include("input_output.jl")
include("link_functions.jl")
include("transfer_matrix.jl")
include("em_algorithm.jl")
include("marginal_model.jl")
include("statistical_summaries.jl")
include("sample_averages.jl")
include("simulations.jl")
include("grp_perm_tests.jl")
include("two_samp_perm_tests.jl")
include("genome_analysis.jl")
include("wgbs.jl")
include("nanopore.jl")
include("multiple_hypothesis.jl")

end # module
