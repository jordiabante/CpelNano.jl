## Deps
using Plots
using StatsBase
using Statistics
using StatsPlots
using LinearAlgebra
using DelimitedFiles

# Default plotting
default(titlefont=(14,"arial"),guidefont=(16,"arial"),tickfont=(12,"arial"))

# Read in data
data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/GM12878/chr22/"
cpel_params = readdlm("$(data_dir)/GM12878_wgbs_cpelnano_theta.txt")[:,4]
cpel_params = map(x->parse.(Float64,split(x,",")),cpel_params)

# Get parameters
a = map(x->x[1],cpel_params)
b = map(x->x[2],cpel_params)
c = map(x->x[3],cpel_params)

a_range = -10.0:0.5:10.0
b_range = -99.5:10.0:99.5
c_range = -20.0:2.0:20.0

a_histo = fit(Histogram,a,a_range)
b_histo = fit(Histogram,b,b_range)
c_histo = fit(Histogram,c,c_range)

a_histo = normalize(a_histo, mode=:probability)
b_histo = normalize(b_histo, mode=:probability)
c_histo = normalize(c_histo, mode=:probability)

sts_a = collect(a_range)
sts_b = collect(b_range)
sts_c = collect(c_range)

cent_a = sts_a.+step(a_range)/2.0
cent_b = sts_b.+step(b_range)/2.0
cent_c = sts_c.+step(c_range)/2.0

p_a = plot(cent_a,a_histo.weights*100,seriestype=:bar,alpha=0.25,xlabel="a",ylabel="Percentage (%)",label="");
p_b = plot(cent_b,b_histo.weights*100,seriestype=:bar,alpha=0.25,xlabel="b",ylabel="Percentage (%)",label="");
p_c = plot(cent_c,c_histo.weights*100,seriestype=:bar,alpha=0.25,xlabel="c",ylabel="Percentage (%)",label="");

p = plot(p_a,p_b,p_c,layout=(3,1),size=(500,800));
savefig(p,"$(data_dir)/Histograms.png")
savefig(p,"$(data_dir)/Histograms.pdf")
