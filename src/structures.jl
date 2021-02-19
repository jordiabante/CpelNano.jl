##################################################################################################
## Methylation calls structure & associated methods
##################################################################################################

# Output files
mutable struct OutputFiles
    theta_file::String                  # THETA
    ex_file::String                     # E[X]
    exx_file::String                    # E[XX]
    mml_file::String                    # MML
    nme_file::String                    # NME
    OutputFiles() = new()
    OutputFiles(out_dir, out_prefix) = new(
        "$(out_dir)/$(out_prefix)_theta.txt",
        "$(out_dir)/$(out_prefix)_ex.txt",
        "$(out_dir)/$(out_prefix)_exx.txt",
        "$(out_dir)/$(out_prefix)_mml.txt",
        "$(out_dir)/$(out_prefix)_nme.txt"
    )
end

# Output differential analysis files
mutable struct OutputDiffFiles
    tmml_file::String                       # File containing Tmml results
    tnme_file::String                       # File containing Tnme results
    tcmd_file::String                       # File containing Tcmd results
    OutputDiffFiles() = new()
    OutputDiffFiles(out_dir, out_prefix) = new(
        "$(out_dir)/$(out_prefix)_tmml.txt",
        "$(out_dir)/$(out_prefix)_tnme.txt",
        "$(out_dir)/$(out_prefix)_tcmd.txt"
    )
end

# CpelNano Configuration
mutable struct CpelNanoConfig
    min_cov::Float64                    # Minimum average coverage
    min_grp_dist::Int64                 # Minimum distance between CpG groups
    max_size_subreg::Int64              # Maximum size subregion
    size_est_reg::Int64                 # Average size estimation region
    max_em_init::Int64                  # Maximum number of EM initializations used
    max_em_iters::Int64                 # Maximum number of iterations in each EM instance
    out_dir::String                     # Output directory
    out_prefix::String                  # Output prefix
    mod_nse::Bool                       # Model noise in methylation calls
    informme_mode::Bool                 # Partition genome in the same way informME does
    matched::Bool                       # Matched or unmatched comparison
    filter::Bool                        # Filter hypothesis for more power
    min_pval::Bool                      # Include regions with not enough data for p-val<0.05
    caller::String                      # Methylation caller used
    bed_reg::String                     # Path to bed file with regions of interest
    verbose::Bool                       # Print intermediate results
    trim::NTuple{4,Int64}               # Trimming of reads
    pe::Bool                            # Paired end (used if BS data)
    out_files::OutputFiles              # Name of output files
    out_diff_files::OutputDiffFiles     # Name of output files
    # Two sample testing
    LMAX_TWO_SAMP::Int64                # Maximum number of permutations in 2 sample test
    pval_comp::Bool                     # Compute p-values
    # Init methods
    CpelNanoConfig() = new(10.0,10,350,3000,10,20,"./","cpelnano",
        true,false,false,false,false,"nanopolish","",
        false,(0, 0, 0, 0),false,OutputFiles("./", "cpelnano"),OutputDiffFiles("./", "cpelnano"),
        100,true)
    CpelNanoConfig(min_cov, max_size_subreg, size_est_reg, max_em_init, max_em_iters) = 
        new(min_cov,10,max_size_subreg + 1,size_est_reg,max_em_init,max_em_iters,"","",
        true,false,false,false,false,"nanopolish","",false,(0, 0, 0, 0),false,OutputFiles("", ""),OutputDiffFiles("", ""),
        100,true)
end

##################################################################################################
## Methylation calls structure & associated methods
##################################################################################################

# BAM file
struct AlignTemp
    strand::String                  # Methylation call strand
    R1::BAM.Record                  # First record from left to right
    R2::BAM.Record                  # Second record from left to right
end

mutable struct AllAlignTemps
    paired::Bool                    # Boolean indicating if record has pair
    templates::Array{AlignTemp,1}   # All alignment templates mapping to a region
end

# Structure for matrices
mutable struct LogGs
    pp::Vector{Float64}                                 # Series of log(g1(+x_p))
    pm::Vector{Float64}                                 # Series of log(g1(-x_p))
    qp::Vector{Float64}                                 # Series of log(g2(+x_q))
    qm::Vector{Float64}                                 # Series of log(g2(-x_q))
    # Init Methods
    LogGs() = new([], [], [], [])
end

# Structure for expectations
mutable struct Expectations
    ex::Vector{Float64}                                 # Vector E[X]
    exx::Vector{Float64}                                # Vector E[XX]
    log_g1::Vector{Float64}                             # Vector E[log g1(xp)]
    log_g2::Vector{Float64}                             # Vector E[log g2(xq)]
    # Init Methods
    Expectations() = new([], [], [], [])
end

# Structure for matrices
mutable struct TransferMat
    u1::Vector{Float64}                                 # Vector u1
    uN::Vector{Float64}                                 # Vector uN
    Ws::Vector{Array{Float64,2}}                        # Series of W matrices
    log_u1::Vector{Float64}                             # Vector log(u1)
    log_uN::Vector{Float64}                             # Vector log(uN)
    log_Ws::Vector{Array{Float64,2}}                    # Series of log(W) matrices
    log_gs::LogGs                                       # log(g_i(±⋅))
    # Init Methods
    TransferMat() = new([], [], [], [], [], [], LogGs())
end

# Structure for methylation call at CpG site
mutable struct MethCallCpgGrp
    obs::Bool                                           # Binary vector indicating if observed
    log_pyx_u::Float64                                  # log p(y|x=-1) as computed by pore model
    log_pyx_m::Float64                                  # log p(y|x=+1) as computed by pore model
    # Init methods
    MethCallCpgGrp() = new(false, NaN, NaN)
    MethCallCpgGrp(log_pyx_u, log_pyx_m) = new(true, log_pyx_u, log_pyx_m)
end

# Structure for CpG groups
mutable struct CpgGrp
    # Coordinates
    grp_int::UnitRange{Int64}                           # Genomic interval of CpG group
    cpg_ind::UnitRange{Int64}                           # CpG indices for given group
    # Init methods
    CpgGrp(grp_int, cpg_ind) = new(grp_int, cpg_ind)
end

# Structure for analysis regions (aka subregions)
mutable struct AnalysisRegions
    num::Int64                                          # Number of analysis regions
    chr_int::Vector{UnitRange{Int64}}                   # Genomic intervals
    cpg_ind::Vector{UnitRange{Int64}}                   # Indices of CpG sites
    # Init methods
    AnalysisRegions() = new(0, [], [])
    AnalysisRegions(chr_int, cpg_ind) = new(length(chr_int), chr_int, cpg_ind)
end

# Structure for methylation observation vectors in a region
mutable struct RegStruct
    # Processed
    proc::Bool                                          # Determines if region has been processed
    # Chromosome properties
    chr::String                                         # Chromosome of the region
    chrst::Int64                                        # Start position of region (1-based)
    chrend::Int64                                       # End position of region (1-based)
    N::Int64                                            # Number of CpG sites in region
    L::Int64                                            # Number of CpG groups in region
    Nl::Vector{Float64}                                 # Number of CpG sites per group
    ρn::Vector{Float64}                                 # Vector of density ρ per CpG site
    ρl::Vector{Float64}                                 # Vector of density ρ per CpG group
    dn::Vector{Float64}                                 # Vector of d(n,n+1) per CpG site
    dl::Vector{Float64}                                 # Vector of d(l,l+1) per CpG group
    cpg_pos::Vector{Int64}                              # CpG site positions
    cpg_grps::Vector{CpgGrp}                            # CpG groups intervals
    # Data
    m::Int64                                            # Number of observations
    calls::Vector{Vector{MethCallCpgGrp}}               # Set of meth calls vecs in a region
    # Analysis regions
    nls_rgs::AnalysisRegions                            # Structure containing analysis regions info    
    # Statistical summaries
    ϕhat::Vector{Float64}                               # Estimated parameter vector ϕ
    logZ::Float64                                       # Log partition function evaluated @ ϕ
    mml::Vector{Float64}                                # Mean methylation level per analysis region
    nme::Vector{Float64}                                # Normalized methylation entropy per analysis region
    # Transfer matrix
    tm::TransferMat                                     # Arrays for transfer matrix methods
    # Expectations
    exps::Expectations                                  # Expectations of estimation region
    # Init Methods
    RegStruct() = new(false,
        "",0,0,0,0,[],[],[],[],[],[],[],
        0,[],
        AnalysisRegions(),
        [],NaN,[],[],
        TransferMat(),
        Expectations()
    )
    RegStruct(calls) = new(false,
        "",0,0,0,0,[],[],[],[],[],[],[],
        length(calls),calls,
        AnalysisRegions(),
        [],NaN,[],[],
        TransferMat(),
        Expectations()
    )
end

# Struct methods
get_depth_ith_cpg(i::Int64,calls::Vector{Vector{MethCallCpgGrp}})::Float64 = sum([x[i].obs for x in calls])
get_ave_depth(reg::RegStruct)::Float64 = sum([get_depth_ith_cpg(i, reg.calls) for i = 1:reg.L]) / reg.L
is_grp_obs(i::Int64,calls::Vector{Vector{MethCallCpgGrp}})::Float64 = any([x[i].obs for x in calls] .== true)
perc_gprs_obs(reg::RegStruct)::Float64 = sum([is_grp_obs(i, reg.calls) for i = 1:reg.L]) / reg.L

##################################################################################################
## Hypothesis testing structs
##################################################################################################

# Structure for analysis regions used in testing
mutable struct NlsRegTestStruct
    tmml_test::Vector{NTuple{2,Float64}}                # Pairs (Tmml,Pmml)
    tnme_test::Vector{NTuple{2,Float64}}                # Pairs (Tnme,Pnme)
    tcmd_test::Vector{NTuple{2,Float64}}                # Pairs (Tcmd,Pcmd)
    # Init Methods
    NlsRegTestStruct(num_nls) = new(fill((NaN, NaN), num_nls), fill((NaN, NaN), num_nls), fill((NaN, NaN), num_nls))
end

# Structure for analysis region used in testing
mutable struct RegStatTestStruct
    chr::String                                         # Chromosome of analysis region
    chrst::Int64                                        # Start position of region (1-based)
    chrend::Int64                                       # End position of region (1-based)
    coords::Vector{UnitRange{Int64}}                    # Coordinates of analysis regions
    tests::NlsRegTestStruct                             # Tests results
    # Init method 1
    RegStatTestStruct() = new()
    # Init method 2
    RegStatTestStruct(est_rs::RegStruct) = new(
        est_rs.chr,
        est_rs.chrst,
        est_rs.chrend,
        est_rs.nls_rgs.chr_int,
        NlsRegTestStruct(est_rs.nls_rgs.num)
    )
        
end