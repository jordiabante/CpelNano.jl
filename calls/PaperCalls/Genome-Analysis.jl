# Deps
using Plots
using Random
using StatsBase
using StatsPlots
using KernelDensity
using DelimitedFiles
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#####################################################################################################
# Produce histogram of group length
#####################################################################################################

# Data dir
data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/Genome-Properties-hg38/"

# Read in data
grp_len_data = readdlm("$(data_dir)/len_group_histogram.txt")[:,1]

# We found that $84.1\%$ (18,770,638) of the CG-groups
grp_len_data[11]
grp_len_data[11] / sum(grp_len_data) * 100

# Whereas the remaining 15.9% (3,549,239) of the CG-groups contain from $12$ to 1318 bases
sum(grp_len_data[12:end])
sum(grp_len_data[12:end]) / sum(grp_len_data) * 100

# Largest group is 1002 nucleotides
findlast(x -> x == 1,grp_len_data .> 0)

# Normalize
grp_len_data = grp_len_data / sum(grp_len_data)
group_size = collect(1:length(grp_len_data))

# Cumululative results
# - 84.06% of them are between 1-11 bp.
# - 97.67% of them are between 1-25 bp.
# - 99.85% of them are between 1-50 bp.
# - 99.99% of them are between 1-100 bp.
map(i -> sum(grp_len_data[1:i] * 100),[11,25,50,100])

#####################################################################################################
# Produce histogram of group number of CGs
#####################################################################################################

# Data dir
data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/Genome-Properties-hg38/"

# Read in data
num_cpg_per_grp_data = readdlm("$(data_dir)/num_cpg_per_grp_histogram.txt")[:,1]

# Group with the most CG sites has 212 CpGs
findlast(x -> x == 1,num_cpg_per_grp_data .> 0)

#####################################################################################################
# Stats on analysis regions (hg38)
#####################################################################################################

# Data dir
data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/Genome-Properties-hg38/"

# Read in data
ar_data = readdlm("$(data_dir)/estimation_region_table_3kb.txt")

# 1. Total number of R regions in the human genome that contain at least 1 CG-group: 919,321
sum(ar_data[:,3] .> 0.0)

# 2. How many of those contain 1 CG-group: 66
sum(ar_data[:,3] .== 1.0)

# 3. The max and min size in # bases of the regions containing 1 CG-group: (4000,4000)
ind = ar_data[:,3] .== 1.0;
minimum(ar_data[ind,1])
maximum(ar_data[ind,1])

# 4. The max and min size in # bases of the regions containing at least 2 CG-groups: (2905,3071)
ind = ar_data[:,3] .>= 2.0;
minimum(ar_data[ind,1])
maximum(ar_data[ind,1])

# 5. the max and min # of CpG sites within the regions that contain 1 CG-group: (1,2)
ind = ar_data[:,3] .== 1.0;
minimum(ar_data[ind,2])
maximum(ar_data[ind,2])

# 6. the max and min # of CpG sites within the regions that contain at least 2 CG-groups
ind = ar_data[:,3] .>= 2.0;
minimum(ar_data[ind,2])
maximum(ar_data[ind,2])

#####################################################################################################
# Produce plots
#####################################################################################################

# Data dir
data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/Genome-Properties-hg38/"

# Read in analysis region data
ar_data = readdlm("$(data_dir)/estimation_region_table_3kb.txt")

# Read in group histogram data
grp_len_data = readdlm("$(data_dir)/len_group_histogram.txt")[:,1]
grp_len_data = grp_len_data / sum(grp_len_data)

# Read in data
num_cpg_per_grp_data = readdlm("$(data_dir)/num_cpg_per_grp_histogram.txt")[:,1]
num_cpg_per_grp_data = num_cpg_per_grp_data / sum(num_cpg_per_grp_data)

## Take random sample
Random.seed!(123);
ran_samp = sample(1:size(ar_data)[1], 10000;replace=false);

## Plot histogram of group size (bp)
p1 = plot(11:40, grp_len_data[11:40] * 100, seriestype=:bar, xlab="CpG Group size (bp)", ylab="Percentage(%)", label="", xlim=(10, 40), ylim=(0, 100));

## Plot histogram of group size (#CGs)
p2 = plot(1:100, num_cpg_per_grp_data[1:100] * 100, seriestype=:bar, xlab="Number of CpGs in groups", ylab="Percentage(%)", label="", xlim=(0, 20), ylim=(0, 100));

## Density plot CpG sites & CpG groups
ind_nonzero = ar_data[:,2] .> 0.0
p3 = plot(ar_data[ind_nonzero,2], seriestype=:density, label="# CpG sites", xlabel="# CpG sites/groups in analysis regions", ylabel="Density");
plot!(p3,ar_data[ind_nonzero,3],seriestype=:density,label="# CpG groups");

## Plot CpG sites vs CpG groups
p4 = plot(ar_data[ran_samp,2], ar_data[ran_samp,3], seriestype=:scatter, xlab="# CpG sites in analysis region", ylab="# CpG groups in analysis region", label="", alpha=0.25);
plot!(p4,1:400,1:400,label="");

# Collage
p = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 750))

savefig(p,"$(data_dir)/Analysis-Region-3kb-hg38.pdf")
