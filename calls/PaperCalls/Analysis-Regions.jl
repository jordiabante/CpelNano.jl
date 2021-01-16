#####################################################################################################
# Plots corresponding to analysis region (or subregions) properties
#####################################################################################################
# Deps
using Plots
using Random
using StatsBase
using StatsPlots
using KernelDensity
using LinearAlgebra
using DelimitedFiles
default(titlefont=(14,"arial"),guidefont=(14,"arial"),tickfont=(10,"arial"),legendfont=(12,"arial"))
blnd_frnd_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]

#####################################################################################################
# Functions
#####################################################################################################
function read_in_data(data_dir::String,sub_max_sizes::Vector{Int64})

    # Ranges
    num_cpgs_rng = 0:1:25
    sub_size_rng = 0:20:600
    sts_sub_size = collect(sub_size_rng)
    sts_num_cpgs = collect(num_cpgs_rng)
    cent_sub_size = sts_sub_size .- step(sub_size_rng)/2.0
    cent_num_cpgs = sts_num_cpgs .- step(num_cpgs_rng)/2.0
    
    # Read in data
    size_histos = []
    cpgs_histos = []
    obs_sub_size = []
    for sub_size in sub_max_sizes
        println("$(sub_size)")
        # Read in data
        sub_data = readdlm("$(data_dir)/subregion_$(sub_size)_study/subregion_table.txt")
        # Get histograms
        sub_size_histo = fit(Histogram,sub_data[:,1],sub_size_rng)
        num_cpgs_histo = fit(Histogram,sub_data[:,2],num_cpgs_rng)
        # Normalize them
        sub_size_histo = normalize(sub_size_histo,mode=:probability)
        num_cpgs_histo = normalize(num_cpgs_histo,mode=:probability)
        # Push histos
        push!(size_histos,sub_size_histo)
        push!(cpgs_histos,num_cpgs_histo)
    end
    
    # Return
    return (cent_sub_size,size_histos),(cent_num_cpgs,cpgs_histos)
    
end
#####################################################################################################
# Produce histogram subregion size
#####################################################################################################

# Data dir
base_dir = "/Users/jordiabante/OneDrive - Johns Hopkins"
data_dir = "$(base_dir)/CpelNano/Data/Real-Data/Analysis-Region-Properties-hg38/"

# Produce histos
sub_max_sizes = [100,150,200,250,300,350,400,450,500]
(cent_sub_size,size_histos),(cent_num_cpgs,cpgs_histos) = read_in_data(data_dir,sub_max_sizes)

# Plot subregion sizes
p1s = []
for (i,size) in enumerate(sub_max_sizes)
    p = plot(cent_sub_size,size_histos[i].weights*100,seriestype=:bar,alpha=0.5,label="",title="$(size) bp",
        ylim=(0,100),xlabel="Analysis region size (bp)",ylabel="Percentage (%)")
    push!(p1s,p)
end

# Store plot
plot(p1s...,size=(1000,1000),layout=(3,3))
savefig("$(data_dir)/Analysis-Region-Sizes.pdf")

# Plot number of CpG sites
p2s = []
for (i,size) in enumerate(sub_max_sizes)
    p = plot(cent_num_cpgs,cpgs_histos[i].weights*100,seriestype=:bar,alpha=0.5,label="",xticks=0:2:25,
        ylim=(0,100),xlabel="Number of CpG sites",ylabel="Percentage (%)",title="$(size) bp",tickfont=(10,"arial"))
    push!(p2s,p)
end

# Store plot
plot(p2s...,size=(1000,1000),layout=(3,3))
savefig("$(data_dir)/Analysis-Region-CpG-Histogram.pdf")
