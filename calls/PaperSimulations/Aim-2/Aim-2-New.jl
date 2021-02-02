#####################################################################################################
## Comparison noise-modeling vs. no noise modeling. Pore model taken from:
##   https://github.com/jts/nanopolish/blob/master/etc/r9-models/r9.4_450bps.cpg.6mer.template.model
#####################################################################################################
using Random
using StatsPlots
using Distributed
using Distributions
# @everywhere using Pkg
# @everywhere Pkg.activate("./")
@everywhere using CpelNano

#####################################################################################################
# Functions
#####################################################################################################
function cos_sim(α1::Vector{Float64}, β1::Vector{Float64}, α2::Vector{Float64}, β2::Vector{Float64})::Float64
    
    # Create aux vars
    aux_1 = vcat(α1, β1)
    aux_2 = vcat(α2, β2)

    # Return cosine similarity
    return aux_1' * aux_2 / (sqrt(sum(aux_1.^2)) * sqrt(sum(aux_2.^2)))

end
function mse(α1::Vector{Float64}, β1::Vector{Float64}, α2::Vector{Float64}, β2::Vector{Float64})::Float64
    
    # Create aux vars
    aux_1 = vcat(α1, β1)
    aux_2 = vcat(α2, β2)

    # Return cosine similarity
    return sum((aux_1 .- aux_2).^2) / length(aux_1)

end
function plot_dists(μ_m::Float64, μ_u::Float64, σs::Vector{Float64})

    # Loop over σs
    ps = []
    for (i, σ) in enumerate(σs)
        # Sigmas
        σ_m = σ_u = σ
        
        # Plot distributions
        if i == 1
            p = plot(Normal(μ_m, σ_m), lw=3, label="methylated", title="\\sigma=$(σ)", xlabel="Current (pA)")
            plot!(p, Normal(μ_u, σ_u), lw=3, label="unmethylated")
            plot!(p, Uniform(μ_m - σ_m, μ_m + σ_m), label="")
            plot!(p, Uniform(μ_u - σ_u, μ_u + σ_u), label="")
        else
            p = plot(Normal(μ_m, σ_m), lw=3, label="", title="\\sigma=$(σ)", xlabel="Current (pA)")
            plot!(p, Normal(μ_u, σ_u), lw=3, label="")
            plot!(p, Uniform(μ_m - σ_m, μ_m + σ_m), label="")
            plot!(p, Uniform(μ_u - σ_u, μ_u + σ_u), label="")
        end

        # Push plot
        push!(ps, p)

    end

    # Return plot array
    return ps

end
function sim(μ_m::Float64, μ_u::Float64, σs::Vector{Float64}, seed::Int64, reps::Int64)
    
    # Seed
    Random.seed!(seed)

    # Parameter vector
    a = 1.0; b = -5.0; c = 2.0; ϕ = [a,b,c]; L = 20;
    
    # Configuration
    config = CpelNano.CpelNanoConfig(); config.max_em_iters = 100; 
    config.max_em_init = 20; config.verbose = false;
    
    # Other parameters
    x = CpelNano.RegStruct(); x.L = L; x.Nl = fill(1.0, x.L); 
    x.ρl = rand(0.001:0.001:0.5, x.L); x.dl = fill(10.0, x.L - 1);
    α, β = CpelNano.get_αβ_from_ϕ(ϕ, x);
    
    # We take 6-mers AACGAA vs AAMGAA
    M = 10; pobs = 0.8;
    
    # Loop over σs
    nse_mat = fill(NaN, (length(σs), reps));
    no_nse_mat = fill(NaN, (length(σs), reps));
    for (i, σ) in enumerate(σs)
        
        CpelNano.print_log("Working on $(σ)")

        # Sigmas
        σ_m = σ_u = σ
    
        # Perform reps times
        Threads.@threads for j = 1:reps

            CpelNano.print_log("Working on iteration j=$(j) (thread $(Threads.threadid()) of out $(Threads.nthreads()))")
            
            # Simulate data
            rs = CpelNano.cpel_samp_ont(M, α, β, pobs, μ_m, σ_m, μ_u, σ_u); 
            rs.L = x.L; rs.Nl = x.Nl; rs.ρl = x.ρl; rs.dl = x.dl; 
            rs_bin = deepcopy(rs); CpelNano.bin_calls!(rs_bin);
    
            # Run EM algorithm
            CpelNano.get_ϕhat!(rs, config)
            # CpelNano.print_log("$(rs.ϕhat)")
            CpelNano.get_ϕhat!(rs_bin, config)
            # CpelNano.print_log("$(rs_bin.ϕhat)")
    
            # Compute cos similarity
            if length(rs.ϕhat) > 0
                nse_α, nse_β = CpelNano.get_αβ_from_ϕ(rs.ϕhat, rs)
            end
            if length(rs_bin.ϕhat) > 0 
                no_nse_α, no_nse_β = CpelNano.get_αβ_from_ϕ(rs_bin.ϕhat, rs_bin)
            end

            # Filter out problematic runs
            if mse(α, β, nse_α, nse_β) > 10.0
                CpelNano.print_log("MSE issue: $(rs.ϕhat)")
                continue
            end

            # Record
            if length(rs.ϕhat) > 0
                nse_mat[i,j] = mse(α, β, nse_α, nse_β)
            end
            if length(rs_bin.ϕhat) > 0
                no_nse_mat[i,j] = mse(α, β, no_nse_α, no_nse_β)
            end


        end
    
    
    end
    
    # Return
    return nse_mat, no_nse_mat

end
function bxplt(nse_mat, no_nse_mat, σs, reps)

    # Cat matrices
    nse_mat_vcat = vcat(nse_mat'...)
    nse_sig_vcat = vcat([fill("$(σ)", reps) for σ in σs]...)
    nse_grp_vcat = fill("M1", length(nse_sig_vcat))
    no_nse_mat_vcat = vcat(no_nse_mat'...)
    no_nse_sig_vcat = vcat([fill("$(σ)", reps) for σ in σs]...)
    no_nse_grp_vcat = fill("M2", length(no_nse_sig_vcat))
    data_vcat = vcat(nse_mat_vcat, no_nse_mat_vcat)
    sig_vcat = vcat(nse_sig_vcat, nse_sig_vcat)
    grp_vcat = vcat(nse_grp_vcat, no_nse_grp_vcat)

    # Remove NaNs
    non_nan_ind = .!isnan.(data_vcat)
    data_vcat = data_vcat[non_nan_ind]
    sig_vcat = sig_vcat[non_nan_ind]
    grp_vcat = grp_vcat[non_nan_ind]

    # Grouped boxplots
    p = groupedboxplot(sig_vcat,data_vcat,group=grp_vcat,outliers=false,
        xlabel="\\sigma",ylabel="MSE",size=(700, 600))

    # Return plot
    return p

end
#####################################################################################################
# Calls
#####################################################################################################

# AAAAACGAAAAA
# AAAAAC => log p(e1|AAAAAC) / AAAAAM => log p(e1|AAAAAM)
# AAAACG => log p(e2|AAAACG) / AAAAMG => log p(e2|AAAAMG)
# ...

# Parameters from pore model (AACGAA - AAMGAA)
μ_m = 101.93
μ_u = 103.13

# Plot current level distribution
# σs = [0.2,0.4,0.6,0.8]
σs = collect(0.1:0.1:0.8)
ps = plot_dists(μ_m, μ_u, σs)
plot(ps...,layout=(4, 2),size=(800, 1000))

# Run simulations
seed = 123
reps = 1000
nse_mat, no_nse_mat = sim(μ_m, μ_u, σs, seed, reps)

# Boxplot
p = bxplt(nse_mat, no_nse_mat, σs, reps);
plot(p)

# sum(nse_mat[1,:] .<=  no_nse_mat[1,:])
# sum(nse_mat[2,:] .<=  no_nse_mat[2,:])
# sum(nse_mat[3,:] .<=  no_nse_mat[3,:])
# sum(nse_mat[4,:] .<=  no_nse_mat[4,:])
