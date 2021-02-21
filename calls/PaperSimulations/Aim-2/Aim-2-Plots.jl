#################################################################################################
# AIM 2
#################################################################################################
## Deps
using Plots
using CodecZlib
using StatsPlots
using Distributions
using DelimitedFiles
using Plots.PlotMeasures

using Pkg
Pkg.activate("./")
using CpelNano

## Constants
const cov_levels = [5.0,10.0,15.0,20.0,25.0]
const cov_labels = ["$(Int(c))x" for c in cov_levels]
const noise_levels = [2.0,2.5,3.0,3.5]
const noise_labels = ["$(s)" for s in noise_levels]
const aim_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-2/"
const gt_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/ground_truth/"
const blind_friend_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#################################################################################################
# Functions
#################################################################################################
function read_in_mod_file(path::String)::Tuple{Vector{Int64},Array{Float64,2}}
    
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
function read_in_mod_file_dic(path::String)::Dict{String,NTuple{2,Vector{Float64}}}
    
    # Get data
    stream = GzipDecompressorStream(open(path))
    in_data = readdlm(stream)
    close(stream)

    # Remove region with no α,β
    kp_ind = [in_data[i,5] != "" && in_data[i,6] != "" for i = 1:size(in_data)[1]]
    in_data = in_data[kp_ind,:]

    # Read in
    cpel_reg = ["$(in_data[i,1])-$(in_data[i,2])-$(in_data[i,3])" for i = 1:size(in_data)[1]]
    αs = string.(in_data[:,5])
    βs = string.(in_data[:,6])
    αs = map(x -> parse.(Float64, split(x, ",")), αs)
    βs = map(x -> parse.(Float64, split(x, ",")), βs)

    # Keep indices with more than two groups
    kp_ind = length.(βs) .> 2
    cpel_reg = cpel_reg[kp_ind]
    αs = αs[kp_ind]
    βs = βs[kp_ind]

    # Create dictionary
    out_dic = Dict{String,NTuple{2,Vector{Float64}}}()
    for (i, reg) in enumerate(cpel_reg)
        out_dic[reg] = (αs[i], βs[i])
    end

    # Return dictionary
    return out_dic

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
function read_in_diff_file(path::String)::Tuple{Vector{Int64},Vector{Float64}}
    
    # Get data
    stream = GzipDecompressorStream(open(path))
    in_data = readdlm(stream)
    close(stream)

    # Read in estimated regions and parameters
    cpel_reg = convert.(Int64, in_data[:,2])
    cpel_stat = convert.(Float64, in_data[:,4])
    
    # Return vectors
    return cpel_reg, cpel_stat

end
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
function mae(vec1::Vector{Float64}, vec2::Vector{Float64})::Vector{Float64}
    
    # Return mean absolute error
    return vec(abs.(vec1 - vec2))

end
function get_EX_EXX_logZ(αs::Vector{Float64}, βs::Vector{Float64})::Tuple{Vector{Float64},Vector{Float64},Float64}

    # Get matrices
    logu1 = CpelNano.get_log_u(αs[1])
    loguL = CpelNano.get_log_u(αs[end])
    logWs = [CpelNano.get_log_W(αs[l], αs[l + 1], βs[l]) for l = 1:(length(αs) - 1)]

    # Get log Z
    logZ = CpelNano.get_log_Z(logu1, loguL, logWs)

    # Get expectations
    ex = CpelNano.get_E_X_log(logu1, loguL, logWs, logZ)
    exx = CpelNano.get_E_XX_log(logu1, loguL, logWs, logZ)

    # Return
    return ex, exx, logZ

end
function comp_cmd(m1::NTuple{2,Vector{Float64}}, m2::NTuple{2,Vector{Float64}})::Float64

    # Create mixture
    αmix = 0.5 * (m1[1] + m2[1])
    βmix = 0.5 * (m1[2] + m2[2])

    # Compute intermediate quantities
    ex1, exx1, logZ1 = get_EX_EXX_logZ(m1[1], m1[2])
    ex2, exx2, logZ2 = get_EX_EXX_logZ(m2[1], m2[2])
    ex_mix, exx_mix, logZ_mix = get_EX_EXX_logZ(αmix, βmix)

    # Compute entropies
    h1 = logZ1 - m1[1]' * ex1 - m1[2]' * exx1
    h2 = logZ2 - m2[1]' * ex2 - m2[2]' * exx2
    c1 = logZ_mix - αmix' * ex1 - βmix' * exx1
    c2 = logZ_mix - αmix' * ex2 - βmix' * exx2

    # Compute CMD
    cmd = 1.0 - (h1 + h2) / (c1 + c2)

    # Return
    return min(1.0, max(0.0, cmd))
    
end
function comp_kl(p::NTuple{2,Vector{Float64}}, q::NTuple{2,Vector{Float64}})::Float64

    # Compute intermediate quantities
    ex_p, exx_p, logZp = get_EX_EXX_logZ(p[1], p[2])
    ex_q, exx_q, logZq = get_EX_EXX_logZ(q[1], q[2])

    # Compute entropy of true distribution p(x)
    hp = logZp - p[1]' * ex_p - p[2]' * exx_p

    # Compute cross-entropy of estimate q(x) wrt true p(x)
    c = logZq - q[1]' * ex_p - q[2]' * exx_p

    # Return H(P,Q)-H(P)
    return c - hp
    
end
function comp_cmd_set_keys(comm_keys::Vector{String}, dic1::Dict{String,NTuple{2,Vector{Float64}}}, dic2::Dict{String,NTuple{2,Vector{Float64}}})::Vector{Float64}

    # Loop over keys
    cmds = map(key -> comp_cmd(dic1[key], dic2[key]), comm_keys)

    # Return
    return cmds
    
end
function comp_kl_set_keys(comm_keys::Vector{String}, dic1::Dict{String,NTuple{2,Vector{Float64}}}, dic2::Dict{String,NTuple{2,Vector{Float64}}})::Vector{Float64}

    # Loop over keys
    kls = map(key -> comp_kl(dic1[key], dic2[key]), comm_keys)

    # Return
    return kls
    
end
function gen_kl_plt_array()
    
    println("Generating: KL plots")

    # Read in ground-truth
    gt_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug.txt.gz"
    gt_αβ_dic = read_in_mod_file_dic(gt_file)

    # Create plot array
    plt_arr = Vector{Any}()
    for (i, sig) in enumerate(noise_levels)
        println("Noise level: $(sig)")   
        # marg_box_mat = fill(NaN, (length(keys(gt_αβ_dic)), length(cov_labels)))
        nonoise_box_mat = fill(NaN, (length(keys(gt_αβ_dic)), length(cov_labels)))
        cpelnano_box_mat = fill(NaN, (length(keys(gt_αβ_dic)), length(cov_labels)))
        for (j, cove) in enumerate(cov_labels)
            println("Coverage level: $(cove)")
            # Get filenames
            # marg_file = "$(aim_dir)/sample_averages/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_theta.txt.gz"
            nonoise_file = "$(aim_dir)/cpelnano_nonoise/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_theta.txt.gz"
            cpelnano_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_theta.txt.gz"
            # Read in files
            # marg_αβ_dic = read_in_mod_file_dic(marg_file)
            nonoise_αβ_dic = read_in_mod_file_dic(nonoise_file)
            cpelnano_αβ_dic = read_in_mod_file_dic(cpelnano_file)
            # Record performance for CpelNano
            comm_keys = collect(intersect(Set(keys(gt_αβ_dic)), Set(keys(cpelnano_αβ_dic))))
            key_ind = [x in comm_keys for x in keys(gt_αβ_dic)]
            cpelnano_box_mat[key_ind,j] .= comp_kl_set_keys(comm_keys, gt_αβ_dic, cpelnano_αβ_dic)
            # Record performance for CpelNano no-noise
            comm_keys = collect(intersect(Set(keys(gt_αβ_dic)), Set(keys(nonoise_αβ_dic))))
            key_ind = [x in comm_keys for x in keys(gt_αβ_dic)]
            nonoise_box_mat[key_ind,j] .= comp_kl_set_keys(comm_keys, gt_αβ_dic, nonoise_αβ_dic)
            # Record performance for marginal
            # comm_keys = collect(intersect(Set(keys(gt_αβ_dic)), Set(keys(marg_αβ_dic))))
            # key_ind = [x in comm_keys for x in keys(gt_αβ_dic)]
            # marg_box_mat[key_ind,j] .= comp_kl_set_keys(comm_keys, gt_αβ_dic, marg_αβ_dic)
        end

        # Arrange CpelNano data
        cpelnano_data = vcat(cpelnano_box_mat...)
        cpelnano_labels = fill("M1", length(cpelnano_data))
        cpelnano_covs = vcat([fill(cove, size(cpelnano_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(cpelnano_data)
        cpelnano_data = cpelnano_data[keepind]
        cpelnano_labels = cpelnano_labels[keepind]
        cpelnano_covs = cpelnano_covs[keepind]

        # Arrange CpelNano no-noise data
        nonoise_data = vcat(nonoise_box_mat...)
        nonoise_labels = fill("M2", length(nonoise_data))
        nonoise_covs = vcat([fill(cove, size(nonoise_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(nonoise_data)
        nonoise_data = nonoise_data[keepind]
        nonoise_labels = nonoise_labels[keepind]
        nonoise_covs = nonoise_covs[keepind]

        # Arrange marginal data
        # marg_data = vcat(marg_box_mat...)
        # marg_labels = fill("CMD (Emp.)", length(marg_data))
        # marg_covs = vcat([fill(cove, size(marg_box_mat)[1]) for cove in cov_labels]...)
        # keepind = .!isnan.(marg_data)
        # marg_data = marg_data[keepind]
        # marg_labels = marg_labels[keepind]
        # marg_covs = marg_covs[keepind]

        # Push boxplot to array
        plt_data = vcat(cpelnano_data, nonoise_data) # , marg_data)
        plt_labels = vcat(cpelnano_labels, nonoise_labels) # , marg_labels)
        plt_covs = vcat(cpelnano_covs, nonoise_covs) # , marg_covs)
        ylim = (0.0, 25)
        if i == 1
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=ylim,
                title="\\sigma=$(noise_labels[i])",xlab="Coverage",ylab="KL",outliers=false)
        else
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=ylim,
                title="\\sigma=$(noise_labels[i])",xlab="Coverage",ylab="KL",label="",outliers=false)
        end
        push!(plt_arr, plt)
        
    end

    # Return array
    return plt_arr

end
function gen_param_plt_array()

    println("Generating: params")

    # Plot base
    xlab = "Coverage"
    ylab = "Cosine similarity"
    palet = fill("#999999", length(cov_labels))
    plt_base = plot(seriestype=:boxplot, color_palette=palet, xticks=(1:5, cov_labels), 
        ylim=(-1, 1), xlab=xlab, ylab=ylab, labels="", alpha=0.5)

    ## Read in ground-truth
    gt_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug.txt.gz"
    gt_sts, gt_params = read_in_mod_file(gt_file)

    # Create plot array
    plt_arr = Vector{Any}()
    for (i, sig) in enumerate(noise_levels)
        println("Noise level: $(sig)")
        plt_aux = deepcopy(plt_base)
        box_mat = fill(NaN, (length(gt_sts), length(cov_labels)))
        for (j, cove) in enumerate(cov_labels)
            # Get filename of sample
            model_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_theta.txt.gz"
            # Read in model file
            nano_sts, nano_params = read_in_mod_file(model_file)
            # Find intersect
            gt_sts_ind, nano_sts_ind = find_comm_regions(gt_sts, nano_sts)
            # Record metric in matrix
            box_mat[gt_sts_ind,j] .= cos_sim(gt_params[gt_sts_ind,:], nano_params[nano_sts_ind,:])
            # Add to plot
            boxplot!(plt_aux, box_mat[gt_sts_ind,j], labels="", title=noise_labels[i], outliers=false)
        end
        # Push boxplot to array
        push!(plt_arr, plt_aux)
    end

    # Return array
    return plt_arr

end
function gen_stat_plt_array(stat, sig)

    cpelnano_data = []
    nonoise_data = []
    cpelnano_cove = []
    nonoise_cove = []
    for (j, cove) in enumerate(cov_labels)
        # Get filenames
        cpelnano_file = "$(aim_dir)/differential/two_sample_comp_gt_vs_cpelnano_$(sig)_$(cov_labels[j])_$(stat).txt.gz"
        nonoise_file = "$(aim_dir)/differential/two_sample_comp_gt_vs_nonoise_$(sig)_$(cov_labels[j])_$(stat).txt.gz"
        # Read in files
        cpelnano_regs, cpelnano_stat = read_in_diff_file(cpelnano_file)
        nonoise_regs, nonoise_stat = read_in_diff_file(nonoise_file)
        # Record performance metric for E[X] in matrix
        append!(cpelnano_data, abs.(cpelnano_stat))
        append!(cpelnano_cove, fill(cove, length(cpelnano_stat)))
        # Record performance metric for E[X] (no noise modeling) in matrix
        append!(nonoise_data, abs.(nonoise_stat))
        append!(nonoise_cove, fill(cove, length(nonoise_stat)))
    end

    # Labels
    cpelnano_labels = fill("M1", length(cpelnano_data))
    nonoise_labels = fill("M2", length(nonoise_data))
    
    # Plot
    plt_data = vcat(cpelnano_data, nonoise_data)
    plt_labels = vcat(cpelnano_labels, nonoise_labels)
    plt_covs = vcat(cpelnano_cove, nonoise_cove)
    ylab = stat == "Tcmd" ? "CMD" : "Absolute Error"
    if stat == "Tcmd"
        yticks = 0.00:0.05:0.15
        ylim = (0.0, 0.15)
    elseif stat == "Tnme"
        yticks = 0.0:0.1:0.50
        ylim = (0.0, 0.50)
    else
        yticks = 0.0:0.05:0.25
        ylim = (0.0, 0.25)
    end
    plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=ylim,yticks=yticks,
        title="$(stat)",xlab="Coverage",ylab=ylab,outliers=false,top_margin=10px, bottom_margin=20px, 
        left_margin=40px, right_margin=20px)
    
    # Return array
    return plt

end
function gen_ey_plt_array()

    println("Generating: E[Y]=E[0.5*(X+1)]")

    ## Read in ground-truth
    gt_ex_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_ex.txt.gz"
    gt_ex_sts, gt_ex = read_in_ex_file(gt_ex_file)
    
    # Create plot array
    plt_arr = Vector{Any}()
    for (i, sig) in enumerate(noise_levels)
        println("Noise level: $(sig)")   
        ex_box_mat = fill(NaN, (length(gt_ex_sts), length(cov_labels)))
        smp_ex_box_mat = fill(NaN, (length(gt_ex_sts), length(cov_labels)))
        nonoise_ex_box_mat = fill(NaN, (length(gt_ex_sts), length(cov_labels)))
        for (j, cove) in enumerate(cov_labels)
            # Get filenames
            ex_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_ex.txt.gz"
            smp_ex_file = "$(aim_dir)/sample_averages/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_ex.txt.gz"
            nonoise_ex_file = "$(aim_dir)/cpelnano_nonoise/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_ex.txt.gz"
            # Read in files
            ex_sts, ex = read_in_ex_file(ex_file)
            smp_ex_sts, smp_exs = read_in_ex_file(smp_ex_file)
            nonoise_ex_sts, nonoise_ex = read_in_ex_file(nonoise_ex_file)
            # Record performance metric for E[X] in matrix
            gt_sts_ind, ex_sts_ind = find_comm_regions(gt_ex_sts, ex_sts)
            ex_box_mat[gt_sts_ind,j] .= mae(gt_ex[gt_sts_ind], ex[ex_sts_ind]) ./ 2.0
            # Record performance metric for E[X] (no noise modeling) in matrix
            gt_sts_ind, nonoise_ex_sts_ind = find_comm_regions(gt_ex_sts, nonoise_ex_sts)
            nonoise_ex_box_mat[gt_sts_ind,j] .= mae(gt_ex[gt_sts_ind], nonoise_ex[nonoise_ex_sts_ind]) ./ 2.0
            # Record performance metric for X̄ in matrix
            gt_sts_ind, smp_ex_sts_ind = find_comm_regions(gt_ex_sts, smp_ex_sts)
            smp_ex_box_mat[gt_sts_ind,j] .= mae(gt_ex[gt_sts_ind], smp_exs[smp_ex_sts_ind]) ./ 2.0
        end

        # Arrange E[X] data
        ex_data = vcat(ex_box_mat...)
        ex_labels = fill("E[X] (M1)", length(ex_data))
        ex_covs = vcat([fill(cove, size(ex_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(ex_data)
        ex_data = ex_data[keepind]
        ex_labels = ex_labels[keepind]
        ex_covs = ex_covs[keepind]

        # Arrange no noise modeling E[X] data
        nonoise_ex_data = vcat(nonoise_ex_box_mat...)
        nonoise_ex_labels = fill("E[X] (M2)", length(nonoise_ex_data))
        nonoise_ex_covs = vcat([fill(cove, size(nonoise_ex_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(nonoise_ex_data)
        nonoise_ex_data = nonoise_ex_data[keepind]
        nonoise_ex_labels = nonoise_ex_labels[keepind]
        nonoise_ex_covs = nonoise_ex_covs[keepind]

        # Arrange X̄ data
        smp_ex_data = vcat(smp_ex_box_mat...)
        smp_ex_labels = fill("Xbar", length(smp_ex_data))
        smp_ex_covs = vcat([fill(cove, size(smp_ex_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(smp_ex_data)
        smp_ex_data = smp_ex_data[keepind]
        smp_ex_labels = smp_ex_labels[keepind]
        smp_ex_covs = smp_ex_covs[keepind]

        # Push boxplot to array
        plt_data = vcat(ex_data, nonoise_ex_data, smp_ex_data)
        plt_labels = vcat(ex_labels, nonoise_ex_labels, smp_ex_labels)
        plt_covs = vcat(ex_covs, nonoise_ex_covs, smp_ex_covs)
        if i == 1
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=(0.0, 0.5),
                title="\\sigma=$(noise_labels[i])",xlab="Coverage",ylab="Absolute error",outliers=false)
        else
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=(0.0, 0.5),
                title="\\sigma=$(noise_labels[i])",xlab="Coverage",ylab="Absolute error",label="",outliers=false)
        end
    push!(plt_arr, plt)

    end

    # Return array
    return plt_arr

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
function gen_eyy_plt_array()

    println("Generating: E[YY']=E[0.5(X+1)0.5(X'+1)]")

    ## Read in ground-truth
    gt_ex_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_ex.txt.gz"
    gt_exx_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_exx.txt.gz"
    gt_ex_cg, gt_ex = read_in_ex_file(gt_ex_file)
    gt_exx_cg1, gt_exx = read_in_exx_file(gt_exx_file)
    gt_eyy = comp_eyy(gt_ex_cg, gt_ex, gt_exx_cg1, gt_exx)
    
    # Create plot array
    plt_arr = Vector{Any}()
    for (i, sig) in enumerate(noise_levels)
        println("Noise level: $(sig)")
        eyy_box_mat = fill(NaN, (length(gt_exx_cg1), length(cov_labels)))
        smp_eyy_box_mat = fill(NaN, (length(gt_exx_cg1), length(cov_labels)))
        nonoise_eyy_box_mat = fill(NaN, (length(gt_exx_cg1), length(cov_labels)))
        for (j, cove) in enumerate(cov_labels)
            println("Coverage: $(cov_labels[j])")
            # Get filenames
            ex_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_ex.txt.gz"
            exx_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_exx.txt.gz"
            smp_ex_file = "$(aim_dir)/sample_averages/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_ex.txt.gz"
            smp_exx_file = "$(aim_dir)/sample_averages/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_exx.txt.gz"
            nonoise_ex_file = "$(aim_dir)/cpelnano_nonoise/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_ex.txt.gz"
            nonoise_exx_file = "$(aim_dir)/cpelnano_nonoise/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_exx.txt.gz"
            # Read in files
            ex_cg1, ex = read_in_ex_file(ex_file)
            smp_ex_cg1, smp_ex = read_in_ex_file(smp_ex_file)
            nonoise_ex_cg1, nonoise_ex = read_in_ex_file(nonoise_ex_file)
            exx_cg1, exx = read_in_exx_file(exx_file)
            smp_exx_cg1, smp_exx = read_in_exx_file(smp_exx_file)
            nonoise_exx_cg1, nonoise_exx = read_in_exx_file(nonoise_exx_file)
            # Compute E[YY] where Y=0.5*(X+1)
            eyy = comp_eyy(ex_cg1, ex, exx_cg1, exx)
            smp_eyy = comp_eyy(smp_ex_cg1, smp_ex, smp_exx_cg1, smp_exx)
            nonoise_eyy = comp_eyy(nonoise_ex_cg1, nonoise_ex, nonoise_exx_cg1, nonoise_exx)
            # Record metric for E[YY] in matrix
            gt_sts_ind, exx_sts_ind = find_comm_regions(gt_exx_cg1, exx_cg1)
            eyy_box_mat[gt_sts_ind,j] .= mae(gt_eyy[gt_sts_ind], eyy[exx_sts_ind])
            # Record metric for E[YY] (no noise) in matrix
            gt_sts_ind, nonoise_exx_sts_ind = find_comm_regions(gt_exx_cg1, nonoise_exx_cg1)
            nonoise_eyy_box_mat[gt_sts_ind,j] .= mae(gt_eyy[gt_sts_ind], nonoise_eyy[nonoise_exx_sts_ind])
            # Record metric for X̄ in matrix
            gt_sts_ind, smp_exx_sts_ind = find_comm_regions(gt_exx_cg1, smp_exx_cg1)
            smp_eyy_box_mat[gt_sts_ind,j] .= mae(gt_eyy[gt_sts_ind], smp_eyy[smp_exx_sts_ind])
        end

        # Arrange E[YY] data
        eyy_data = vcat(eyy_box_mat...)
        eyy_labels = fill("E[XX] (M1)", length(eyy_data))
        eyy_covs = vcat([fill(cove, size(eyy_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(eyy_data)
        eyy_data = eyy_data[keepind]
        eyy_labels = eyy_labels[keepind]
        eyy_covs = eyy_covs[keepind]

        # Arrange E[YY] (no noise) data
        nonoise_eyy_data = vcat(nonoise_eyy_box_mat...)
        nonoise_eyy_labels = fill("E[XX] (M2)", length(nonoise_eyy_data))
        nonoise_eyy_covs = vcat([fill(cove, size(nonoise_eyy_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(nonoise_eyy_data)
        nonoise_eyy_data = nonoise_eyy_data[keepind]
        nonoise_eyy_labels = nonoise_eyy_labels[keepind]
        nonoise_eyy_covs = nonoise_eyy_covs[keepind]

        # Arrange X̄X̄ data
        smp_eyy_data = vcat(smp_eyy_box_mat...)
        smp_eyy_labels = fill("XXbar", length(smp_eyy_data))
        smp_eyy_covs = vcat([fill(cove, size(smp_eyy_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(smp_eyy_data)
        smp_eyy_data = smp_eyy_data[keepind]
        smp_eyy_labels = smp_eyy_labels[keepind]
        smp_eyy_covs = smp_eyy_covs[keepind]

        # Push boxplot to array
        plt_data = vcat(eyy_data, nonoise_eyy_data, smp_eyy_data)
        plt_covs = vcat(eyy_covs, nonoise_eyy_covs, smp_eyy_covs)
        plt_labels = vcat(eyy_labels, nonoise_eyy_labels, smp_eyy_labels)
        if i == 1
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=(0.0, 0.5),
                title="\\sigma=$(noise_labels[i])",xlab="Coverage",ylab="Absolute error",outliers=false)
        else
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=(0.0, 0.5),
                title="\\sigma=$(noise_labels[i])",xlab="Coverage",ylab="Absolute error",label="",outliers=false)
        end
        push!(plt_arr, plt)
        
    end

    # Return array
    return plt_arr

end
function plt_params_scatter(cove, sig)

    # Find index of coverage label
    cov_ind = findfirst(isequal(cove), cov_levels)

    # Files
    gt_mod_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug.txt.gz"
    nano_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[cov_ind])_theta.txt.gz"
    
    # Params
    wgbs_sts, wgbs_params = read_in_mod_file(gt_mod_file)
    nano_sts, nano_params = read_in_mod_file(nano_file)

    # Find common regions
    wgbs_sts_ind, nano_sts_ind = find_comm_regions(wgbs_sts, nano_sts)
    wgbs_params = wgbs_params[wgbs_sts_ind,:]
    nano_params = nano_params[nano_sts_ind,:]

    # Plot scatter
    xlab = "Fine-grain"
    ylab = "Coarse-grain"
    p_a = plot(wgbs_params[:,1],nano_params[:,1],seriestype=:scatter,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(-10.0, 10.0),ylim=(-10.0, 10.0),label="",title="a");
    p_b = plot(wgbs_params[:,2],nano_params[:,2],seriestype=:scatter,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(-100.0, 100.0),ylim=(-100.0, 100.0),label="",title="b");
    p_c = plot(wgbs_params[:,3],nano_params[:,3],seriestype=:scatter,xlabel=xlab,ylabel=ylab,
        markersize=0.2,markeralpha=0.1,xlim=(-20.0, 20.0),ylim=(-20.0, 20.0),label="",title="c");
    plot(p_a, p_b, p_c, layout=(3, 1), size=(500, 800))
    savefig("$(aim_dir)/Scatter-Params-cov-$(cove)-s-$(sig)-Aim-2.pdf")
    savefig("$(aim_dir)/Scatter-Params-cov-$(cove)-s-$(sig)-Aim-2.png")

    # Return nothing
    return nothing

    end
function get_htmp_mat(vec1, vec2, axis_rng)

    outmat = fill(0.0, (length(axis_rng) - 1, length(axis_rng) - 1))
    for i in (length(axis_rng) - 1):-1:1
            for j in 1:(length(axis_rng) - 1)
            ind_1 = [axis_rng[i] <= x < axis_rng[i + 1] for x in vec1]
            ind_2 = [axis_rng[j] <= x < axis_rng[j + 1] for x in vec2]
            outmat[i,j] = sum(ind_1 .& ind_2)
        end
    end
    
    outmat .+= 1
    outmat ./= sum(outmat)
outmat = -log10.(outmat)
    
    return outmat

end
function plt_params_heatmap(cove, sig)

    # Find index of coverage label
    cov_ind = findfirst(isequal(cove), cov_levels)

    # Files
    gt_mod_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug.txt.gz"
    nano_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[cov_ind])_theta.txt.gz"
    
    # Params
    wgbs_sts, wgbs_params = read_in_mod_file(gt_mod_file)
    nano_sts, nano_params = read_in_mod_file(nano_file)

    # Find common regions
    wgbs_sts_ind, nano_sts_ind = find_comm_regions(wgbs_sts, nano_sts)
    wgbs_params = wgbs_params[wgbs_sts_ind,:]
    nano_params = nano_params[nano_sts_ind,:]

    # Ranges
    a_axis_rng = collect(-10.0:1:10.0)
    b_axis_rng = collect(-100.0:10:100.0)
    c_axis_rng = collect(-20.0:2:20.0)

    # Construct matrix
    a_mat = get_htmp_mat(nano_params[:,1], wgbs_params[:,1], a_axis_rng)
    b_mat = get_htmp_mat(nano_params[:,2], wgbs_params[:,2], b_axis_rng)
    c_mat = get_htmp_mat(nano_params[:,3], wgbs_params[:,3], c_axis_rng)

    # Plot scatter
    xlab = "Fine-grain"
    ylab = "Coarse-grain"
    a_axis_labs = [isodd(i) ? "$(a_axis_rng[i])" : "" for i = 1:(length(a_axis_rng))]
    a_axis_ticks = (1:length(a_axis_labs), a_axis_labs)
    b_axis_labs = [isodd(i) ? "$(b_axis_rng[i])" : "" for i = 1:(length(b_axis_rng))]
    b_axis_ticks = (1:length(b_axis_labs), b_axis_labs)
    c_axis_labs = [isodd(i) ? "$(c_axis_rng[i])" : "" for i = 1:(length(c_axis_rng))]
    c_axis_ticks = (1:length(c_axis_labs), c_axis_labs)
    p_a = heatmap(a_mat,xlabel=xlab, ylabel=ylab, xticks=a_axis_ticks, yticks=a_axis_ticks,
        label="", c=cgrad([:white, :black]), clims=(1, 4), xrotation=45,
        title="\\alpha (\\rho=$(round(cor(nano_params[:,1], wgbs_params[:,1]);digits=2)))");
    p_b = heatmap(b_mat,xlabel=xlab, ylabel=ylab, xticks=b_axis_ticks, yticks=b_axis_ticks,
        label="", c=cgrad([:white, :black]), clims=(1, 4), xrotation=45, 
        title="\\beta (\\rho=$(round(cor(nano_params[:,2], wgbs_params[:,2]);digits=2)))");
    p_c = heatmap(c_mat,xlabel=xlab, ylabel=ylab, xticks=c_axis_ticks, yticks=c_axis_ticks,
        label="", c=cgrad([:white, :black]), clims=(1, 4), xrotation=45, 
        title="\\gamma (\\rho=$(round(cor(nano_params[:,3], wgbs_params[:,3]);digits=2)))");
    plot(p_a, p_b, p_c, layout=(3, 1), size=(800, 1500),top_margin=10px, bottom_margin=20px, 
        left_margin=40px, right_margin=20px)
    savefig("$(aim_dir)/Heatmap-Params-cov-$(cove)-s-$(sig)-Aim-2.pdf")

    # Return nothing
    return nothing

end
function plt_exp_heatmap(cove, sig)

    # Ground-truth quantities
    gt_exx_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_exx.txt.gz"
    gt_ex_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_ex.txt.gz"
    gt_ex_sts, gt_ex = read_in_ex_file(gt_ex_file)
    gt_exx_sts, gt_exx = read_in_exx_file(gt_exx_file)
    gt_ey = 0.5 .* (gt_ex .+ 1.0)
    gt_eyy = comp_eyy(gt_ex_sts, gt_ex, gt_exx_sts, gt_exx)

    # Find index of coverage label
    cov_ind = findfirst(isequal(cove), cov_levels)

    # Estimated quantities
    ex_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[cov_ind])_ex.txt.gz"
    exx_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[cov_ind])_exx.txt.gz"
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

    # Construct matrix
    res_grid = 0.05
    axis_rng = collect(0:res_grid:1.0)
    mat_y = fill(0.0, (length(axis_rng) - 1, length(axis_rng) - 1))
    mat_yy = fill(0.0, (length(axis_rng) - 1, length(axis_rng) - 1))
    for i in (length(axis_rng) - 1):-1:1
            for j in 1:(length(axis_rng) - 1)
            ind_1 = [axis_rng[i] <= x < axis_rng[i + 1] for x in ey]
            ind_2 = [axis_rng[j] <= x < axis_rng[j + 1] for x in gt_ey]
            mat_y[i,j] = sum(ind_1 .& ind_2)
            ind_1 = [axis_rng[i] <= x < axis_rng[i + 1] for x in eyy]
            ind_2 = [axis_rng[j] <= x < axis_rng[j + 1] for x in gt_eyy]
    mat_yy[i,j] = sum(ind_1 .& ind_2)
    end
    end
    
    # Plot
    mat_y .+= 1
    mat_yy .+= 1
    mat_y ./= sum(mat_y)
    mat_yy ./= sum(mat_yy)
    mat_y_log10 = -log10.(mat_y)
    mat_yy_log10 = -log10.(mat_yy)
xlab = "Fine-grain (WGBS)"
    ylab = "Coarse-grain (Nanopore)"
    axis_labs = [isodd(i) ? "$(axis_rng[i])" : "" for i = 1:(length(axis_rng))]
    axis_ticks = (1:length(axis_labs), axis_labs)
    p_ex = heatmap(mat_y_log10, xlabel=xlab, ylabel=ylab, xticks=axis_ticks, yticks=axis_ticks,
        label="", c=cgrad([:white, :black]), clims=(1, 5), xrotation=45,
        title="E[X] (\\rho=$(round(cor(gt_ey, ey);digits=2)))");
    p_exx = heatmap(mat_yy_log10, xlabel=xlab, ylabel=ylab, xticks=axis_ticks, yticks=axis_ticks,
        label="", c=cgrad([:white, :black]), clims=(1, 5), xrotation=45, 
        title="E[XX] (\\rho=$(round(cor(gt_eyy, eyy);digits=2)))");
    plot(p_ex, p_exx, layout=(1, 2), size=(1100, 500), top_margin=10px, bottom_margin=20px, 
        left_margin=20px, right_margin=20px)
    savefig("$(nano_dir)/Heatmap-Exp-cov-$(cove)-s-$(sig)-Aim-2-Log-Res-$(res_grid).pdf")

    # Return nothing
    return nothing

end
#################################################################################################
# Calls
#################################################################################################

## Parameter estimates

# Generate plot array
plt_arr = gen_param_plt_array()

# Make plot
plot(plt_arr...,layout=(2, 2),size=(1000, 1000))
savefig("$(aim_dir)/Boxplots-Params-Aim-2.pdf")

# # Generate plot array
# plt_arr = gen_kl_plt_array()

# # Make plot
# plot(plt_arr...,layout=(2, 2),size=(1000, 1000))
# savefig("$(aim_dir)/Boxplots-KL-Params-Aim-2.pdf")

## E[X] vs X̄

# # Generate plot array
# plt_arr = gen_ey_plt_array()

# # Make plot
# plot(plt_arr...,layout=(2, 2),size=(1000, 1000))
# savefig("$(aim_dir)/Boxplots-EX-Aim-2.pdf")

## E[XX] vs X̄

# # Generate plot array
# plt_arr = gen_eyy_plt_array()

# # Make plot
# plot(plt_arr...,layout=(2, 2),size=(1000, 1000))
# savefig("$(aim_dir)/Boxplots-EXX-Aim-2.pdf")

## Heatmap plots for specific noise and cov

# Heatmap parameters
println("Generating: Scatter parameters")
plt_params_heatmap(10.0,3.0)
plt_params_heatmap(20.0,3.0)

# Heatmap Expectations
println("Generating: Scatter expectations")
plt_exp_heatmap(10.0,2.0)
plt_exp_heatmap(20.0,2.0)
plt_exp_heatmap(10.0,3.0)
plt_exp_heatmap(20.0,3.0)
plt_exp_heatmap(10.0,3.5)
plt_exp_heatmap(20.0,3.5)

## Analysis regions

Tmml_plt = gen_stat_plt_array("Tmml", 3.0)
Tnme_plt = gen_stat_plt_array("Tnme", 3.0)
Tcmd_plt = gen_stat_plt_array("Tcmd", 3.0)
plot(Tmml_plt,Tnme_plt,Tcmd_plt,layout=(3, 1),size=(1000, 2000))
savefig("$(aim_dir)/Boxplots-Analysis-Regions-Aim-2.pdf")
