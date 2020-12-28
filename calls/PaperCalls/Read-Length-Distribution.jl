# Deps
using Plots
using Random
using StatsBase
using StatsPlots
using Distributions
using KernelDensity
using LinearAlgebra
using DelimitedFiles
default(titlefont=(14,"arial"),guidefont=(16,"arial"),tickfont=(12,"arial"))

######################################################################################
# Functions
######################################################################################
function get_percs(d::UnivariateDistribution,centers::Vector{Float64})

    # Get intervals
    win_size = centers[2]-centers[1]
    intv_vec = [(0.5+win_size*i):(0.5+win_size*(i+1)) for i=0:(length(centers)-1)]

    # Loop over centers
    percs = zeros(length(intv_vec))
    for i=1:length(percs)
        if i==1
            percs[1] = cdf(d,maximum(intv_vec[1]))
        else
            percs[i] = cdf(d,maximum(intv_vec[i]))-cdf(d,maximum(intv_vec[i-1]))
        end
    end

    # Return percentages
    return percs

end
function plot_dist(lens,reso)

    # Histograms
    len_range = 0.0:reso:50000.0
    sts = collect(len_range)
    centers = sts.+step(len_range)/2.0
    histo = fit(Histogram,lens,len_range)
    histo = normalize(histo, mode=:probability)

    # Plot histogram
    len_plot = plot(sts,histo.weights.*100,seriestype=:bar,bins=100,xlim=(0,maximum(len_range)),
        alpha=0.5,xlabel="Read length (bp)",ylabel="Percentage(%)",label="Observed",size=(700,400));

    # Fit distributions
    # gam_fit = fit_mle(Gamma,lens)
    # gam_logpdf = sum(logpdf.(gam_fit,lens))
    # println("Gamma: $(gam_fit.α), $(gam_fit.θ). AIC: $(4-2*gam_logpdf)")
    exp_fit = fit_mle(Exponential,lens)
    exp_logpdf = sum(logpdf.(exp_fit,lens))
    println("Exp: $(exp_fit.θ). AIC: $(2-2*exp_logpdf)")
    # nb_fit = fit_mle(NegativeBinomial,lens)
    # nb_logpdf = sum(logpdf.(nb_fit,lens))
    # println("NB: $(nb_fit.r),$(nb_fit.p). AIC: $(4-2*nb_logpdf)")
    # gam_perc = get_percs(gam_fit,centers)
    exp_perc = get_percs(exp_fit,centers)
    # nb_perc = get_percs(nb_fit,centers)

    # Plot density
    plot!(len_plot,centers,exp_perc.*100,seriestype=:scatter,lw=2,label="Exp");
    # plot!(len_plot,centers,gam_perc.*100,seriestype=:scatter,lw=2,label="Gamma");
    # plot!(len_plot,centers,nb_perc.*100,seriestype=:scatter,lw=2,label="NB");

    # Return plot
    return len_plot

end
######################################################################################
# Raw lengths
######################################################################################

# IO
data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Read-Length/"

# Read in histogram data
gm_lens = readdlm("$(data_dir)/gm12878_rel_6_read_lens.txt")[:,1]

# Histograms
p1 = plot_dist(gm_lens,500)
savefig(p1,"$(data_dir)/gm12878_read_lens_plot.png")
savefig(p1,"$(data_dir)/gm12878_read_lens_plot.pdf")

######################################################################################
# Simulated lengths
######################################################################################

# IO
data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Read-Length/"

# Read in histogram data
read_len = readdlm("$(data_dir)/simulated_read_lengths.txt")[:,1]

# Histograms
p1 = plot_dist(read_len,500)
savefig(p1,"$(data_dir)/simulated_read_lengths_plot.png")
savefig(p1,"$(data_dir)/simulated_read_lengths_plot.pdf")
