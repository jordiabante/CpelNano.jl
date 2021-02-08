#################################################################################################
# AIM 3:
#################################################################################################
## Deps
using StatsPlots
using Distributions
using DelimitedFiles

## Constants
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/GM12878/"
const aim_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-3/"

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#################################################################################################
# Functions
#################################################################################################
function find_comm_regions(sts1::Vector{Int64}, sts2::Vector{Int64})::NTuple{2,Vector{Bool}}

    # Find common regions
    int_sts = intersect(Set(sts1), Set(sts2))
    sts1_ind = map(i -> sts1[i] in int_sts, 1:length(sts1))
    sts2_ind = map(i -> sts2[i] in int_sts, 1:length(sts2))

    # Return
    return sts1_ind, sts2_ind

end
function cos_sim(mat1::Array{Float64,2}, mat2::Array{Float64,2})::Vector{Float64}
    
    # Intermediate quantities
    num = sum(mat1 .* mat2, dims=2)
    den = sqrt.(sum(mat1.^2, dims=2)) .* sqrt.(sum(mat2.^2, dims=2))

    # Return cosine similarity
    return vec(num ./ den)

end
function read_in_model_file(path::String)::Tuple{Vector{Int64},Array{Float64,2}}
    
    # Get data
    stream = GzipDecompressorStream(open(path))
    in_data = readdlm(stream)
    close(stream)

    # Read in estimated regions and parameters
    cpel_reg = convert.(Int64, in_data[:,2])
    cpel_params = in_data[:,4]
    cpel_params = map(x -> parse.(Float64, split(x, ",")), cpel_params)

    # Get parameters
    a = map(x -> x[1], cpel_params)
    b = map(x -> x[2], cpel_params)
    c = map(x -> x[3], cpel_params)
    
    # Return vectors
    return cpel_reg, [a b c]

end
function read_in_ex_file(path::String)::Tuple{Vector{Int64},Vector{Float64}}
    
    # Get data
    stream = GzipDecompressorStream(open(path))
    in_data = readdlm(stream)
    close(stream)

    # Read in estimated regions and parameters
    cpel_reg = convert.(Int64, in_data[:,2])
    cpel_ex = convert.(Float64, in_data[:,4])
    
    # Return vectors
    return cpel_reg, cpel_ex

end
function read_in_exx_file(path::String)::Tuple{Vector{Int64},Vector{Float64}}
    
    # Get data
    stream = GzipDecompressorStream(open(path))
    in_data = readdlm(stream)
    close(stream)

    # Read in estimated regions and parameters
    cpel_cg1 = convert.(Int64, in_data[:,2])
    cpel_exx = convert.(Float64, in_data[:,4])
    
    # Return vectors
    return cpel_cg1, cpel_exx

end
function comp_eyy(ex_cg::Vector{Int64}, ex::Vector{Float64}, exx_cg1::Vector{Int64}, exx::Vector{Float64})::Vector{Float64}
    
    # Compute E[YnYn+1]
    eyy = fill(NaN, length(exx))
    @inbounds for n = 1:length(eyy)
        # Find corresponding CG sites
        cg1_ind = findfirst(isequal(exx_cg1[n]), ex_cg)
        isnothing(cg1_ind) && continue
        # E[YY] = 0.25 *(E[XnXn+1]+E[Xn]+E[Xn+1]+1)
        eyy[n] = 0.25 * (exx[n] + ex[cg1_ind] + ex[cg1_ind + 1] + 1.0)
    end

    # Return 
    return min.(1.0, max.(0.0, eyy))

end
function plt_params_scatter()

    # Files
    gt_mod_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug.txt.gz"
    nano_file = "$(aim_dir)/cpelnano/gm12878_chr22_cpelnano_theta.txt.gz"
    
    # Params
    wgbs_sts, wgbs_params = read_in_model_file(gt_mod_file)
    nano_sts, nano_params = read_in_model_file(nano_file)

    # Find common regions
    wgbs_sts_ind, nano_sts_ind = find_comm_regions(wgbs_sts, nano_sts)
    wgbs_params = wgbs_params[wgbs_sts_ind,:]
    nano_params = nano_params[nano_sts_ind,:]

    # Plot scatter
    xlab = "Fine-grain (WGBS)"
    ylab = "Coarse-grain (Nanopore)"
    p_a = plot(wgbs_params[:,1],nano_params[:,1],seriestype=:scatter,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(-10.0, 10.0),ylim=(-10.0, 10.0),label="",title="a");
    p_b = plot(wgbs_params[:,2],nano_params[:,2],seriestype=:scatter,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(-100.0, 100.0),ylim=(-100.0, 100.0),label="",title="b");
    p_c = plot(wgbs_params[:,3],nano_params[:,3],seriestype=:scatter,xlabel=xlab,ylabel=ylab,
        markersize=0.2,markeralpha=0.1,xlim=(-20.0, 20.0),ylim=(-20.0, 20.0),label="",title="c");
    plot(p_a, p_b, p_c, layout=(3, 1), size=(500, 800))
    savefig("$(aim_dir)/Scatter-Params-Aim-3.pdf")
    savefig("$(aim_dir)/Scatter-Params-Aim-3.png")

    # Return nothing
    return nothing

end
function plt_exp_scatter()

    # Ground-truth quantities
    gt_exx_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_exx.txt.gz"
    gt_ex_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_ex.txt.gz"
    gt_ex_sts, gt_ex = read_in_ex_file(gt_ex_file)
    gt_exx_sts, gt_exx = read_in_exx_file(gt_exx_file)
    gt_ey = 0.5 .* (gt_ex .+ 1.0)
    gt_eyy = comp_eyy(gt_ex_sts, gt_ex, gt_exx_sts, gt_exx)

    # Estimated quantities
    ex_file = "$(aim_dir)/cpelnano/gm12878_chr22_cpelnano_ex.txt.gz"
    exx_file = "$(aim_dir)/cpelnano/gm12878_chr22_cpelnano_exx.txt.gz"
    ex_sts, ex = read_in_ex_file(ex_file)
    exx_sts, exx = read_in_exx_file(exx_file)
    ey = 0.5 .* (ex .+ 1.0)
    eyy = comp_eyy(ex_sts, ex, exx_sts, exx)

    # Find common regions E[X]
    wgbs_sts_ex_ind, ex_sts_ind = find_comm_regions(gt_ex_sts, ex_sts)
    gt_ey = gt_ey[wgbs_sts_ex_ind]
    ey = ey[ex_sts_ind]
    
    # Find common regions E[XX]
    wgbs_sts_exx_ind, exx_sts_ind = find_comm_regions(gt_exx_sts, exx_sts)
    gt_eyy = gt_eyy[wgbs_sts_exx_ind]
    eyy = eyy[exx_sts_ind]

    # Plot scatter
    xlab = "Fine-grain (WGBS)"
    ylab = "Coarse-grain (Nanopore)"
    p_ex = plot(gt_ey,ey,seriestype=:scatter,xlabel=xlab,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(0.0, 1.0),ylim=(0.0, 1.0),label="",title="E[X]");
    p_exx = plot(gt_eyy,eyy,seriestype=:scatter,xlabel=xlab,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(0.0, 1.0),ylim=(0.0, 1.0),label="",title="E[XX]");
    plot(p_ex, p_exx, layout=(2, 1), size=(500, 600))
    savefig("$(aim_dir)/Scatter-Exp-Aim-3.pdf")
    savefig("$(aim_dir)/Scatter-Exp-Aim-2.png")

    # Return nothing
    return nothing

end
#################################################################################################
# Calls
#################################################################################################

# Scatter parameters
println("Generating: Scatter parameters")
plt_params_scatter()

# Scatter Expectations
println("Generating: Scatter expectations")
plt_exp_scatter()
