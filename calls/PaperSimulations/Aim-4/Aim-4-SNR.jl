#################################################################################################
# AIM 4: signal-to-noise ratio study of nanopolish model file
#################################################################################################
## Deps
using Pkg
Pkg.activate("./") # <=
using CpelNano
using StatsBase
using StatsPlots
using Distributed
using LinearAlgebra
using Distributions
using DelimitedFiles

## Constants

# IO
const aim_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-4/" # <=

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#################################################################################################
# Calls
#################################################################################################

# IO
dataDir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/Nanopolish-Pore-Model/"
nano = "$(dataDir)/r9.4_450bps.cpg.6mer.template.model"

# Get nanopolish model dictionary
nano_dic = CpelNano.get_nano_pore_mod(nano)

# Plot histogram of Δμ
Δμs = [nano_dic[key][1] - nano_dic[key][3] for key in keys(nano_dic)]
Δμ_range = -5.0:0.2:5.0
Δμ_histo = fit(Histogram, Δμs, Δμ_range)
Δμ_histo = normalize(Δμ_histo, mode=:probability)
sts_Δμ = collect(Δμ_range)
cent_Δμ = sts_Δμ .+ step(Δμ_range) / 2.0
p_Δμ = plot(cent_Δμ, Δμ_histo.weights * 100, seriestype=:bar, alpha=0.25, xlabel="\\Delta\\mu", ylabel="Percentage (%)", label="");

# Plot histogram of σ̄
σ̄s = [0.5 * (nano_dic[key][2] + nano_dic[key][4]) for key in keys(nano_dic)]
σ̄_range = 0.0:0.1:5.0
σ̄_histo = fit(Histogram, σ̄s, σ̄_range)
σ̄_histo = normalize(σ̄_histo, mode=:probability)
sts_σ̄ = collect(σ̄_range)
cent_σ̄ = sts_σ̄ .+ step(σ̄_range) / 2.0
p_σ̄ = plot(cent_σ̄, σ̄_histo.weights * 100, seriestype=:bar, alpha=0.25, xlabel="Average \\sigma", ylabel="Percentage (%)", label="");

# Plot histogram of SNR
snrs = abs.(Δμs ./ σ̄s)
snr_range = 0.0:0.1:5.0
snr_histo = fit(Histogram, snrs, snr_range)
snr_histo = normalize(snr_histo, mode=:probability)
sts_snr = collect(snr_range)
cent_snr = sts_snr .+ step(snr_range) / 2.0
p_snr = plot(cent_snr, snr_histo.weights * 100, seriestype=:bar, alpha=0.25, xlabel="SNR", ylabel="Percentage (%)", label="");

# Plot all of them and save
plot(p_Δμ,p_σ̄,p_snr,layout=(3, 1),size=(600, 800))
savefig("$(aim_dir)/Histogram-Nanopolish-r9-model-file.pdf")

# Median SNR
snr_med = median(snrs)

# Δμ for simulations
Δμ_sims_1 = snr_med * 3.0 / sqrt(4000 / 450)
Δμ_sims_2 = median(Δμs)
