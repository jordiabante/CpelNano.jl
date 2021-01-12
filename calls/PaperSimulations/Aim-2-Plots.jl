#################################################################################################
# AIM 2
#################################################################################################
## Deps
using Plots
using StatsPlots
using Distributions
using DelimitedFiles

## Constants
const cov_levels = [5.0,10.0,15.0,20.0,25.0]
const noise_levels = [0.5,1.0,1.5,2.0,2.5,3.0]
const cov_labels = ["5x","10x","15x","20x","25x"]
const noise_labels = ["0.5","1.0","1.5","2.0","2.5","3.0"]
const aim_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-2/"
const gt_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/ground_truth/"
const blind_friend_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]

## Default attributes
default(titlefont=(14,"arial"),guidefont=(16,"arial"),tickfont=(12,"arial"))

#################################################################################################
# Functions
#################################################################################################
function read_in_model_file(path::String)::Tuple{Vector{Int64},Array{Float64,2}}
    
    # Read in estimated regions and parameters
    cpel_reg = convert.(Int64,readdlm(path)[:,2])
    cpel_params = readdlm(path)[:,4]
    cpel_params = map(x->parse.(Float64,split(x,",")),cpel_params)

    # Get parameters
    a = map(x->x[1],cpel_params)
    b = map(x->x[2],cpel_params)
    c = map(x->x[3],cpel_params)
    
    # Return vectors
    return cpel_reg,[a b c]

end
function read_in_exp_file(path::String)::Tuple{Vector{Int64},Vector{Float64}}
    
    # Read in estimated regions and parameters
    cpel_reg = convert.(Int64,readdlm(path)[:,2])
    cpel_exp = convert.(Float64,readdlm(path)[:,4])
    
    # Return vectors
    return cpel_reg,cpel_exp

end
function find_comm_regions(sts1::Vector{Int64},sts2::Vector{Int64})::NTuple{2,Vector{Bool}}

    # Find common regions
    int_sts = intersect(Set(sts1),Set(sts2))
    sts1_ind = map(i-> sts1[i] in int_sts,1:length(sts1))
    sts2_ind = map(i-> sts2[i] in int_sts,1:length(sts2))

    # Return
    return sts1_ind,sts2_ind

end
function cos_sim(mat1::Array{Float64,2},mat2::Array{Float64,2})::Vector{Float64}
    
    # Intermediate quantities
    num = sum(mat1 .* mat2,dims=2)
    den = sqrt.(sum(mat1.^2,dims=2)) .* sqrt.(sum(mat2.^2,dims=2))

    # Return cosine similarity
    return vec(num./den)

end
function mse(vec1::Vector{Float64},vec2::Vector{Float64})::Vector{Float64}
    
    # Return mean-squared error
    return vec((vec1-vec2).^2)

end
function mae(vec1::Vector{Float64},vec2::Vector{Float64})::Vector{Float64}
    
    # Return mean absolute error
    return vec(abs.(vec1-vec2))

end
function gen_param_plt_array()

    println("Generating: params")

    # Plot base
    xlab = "Coverage"
    ylab="Cosine similarity"
    palet = fill("#999999",length(cov_labels))
    plt_base = plot(seriestype=:boxplot,color_palette=palet,xticks=(1:5,cov_labels),ylim=(-1,1),xlab=xlab,ylab=ylab,labels="",alpha=0.5)

    ## Read in ground-truth
    gt_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug.txt"
    gt_sts,gt_params = read_in_model_file(gt_file)

    # Create plot array
    plt_arr = Vector{Any}()
    for (i,sig) in enumerate(noise_levels)
        println("Noise level: $(sig)")
        plt_aux = deepcopy(plt_base)
        box_mat = fill(NaN,(length(gt_sts),length(cov_labels)))
        for (j,cove) in enumerate(cov_labels)
            # Get filename of sample
            model_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_theta.txt"
            # Read in model file
            nano_sts,nano_params = read_in_model_file(model_file)
            # Find intersect
            gt_sts_ind,nano_sts_ind = find_comm_regions(gt_sts,nano_sts)
            # Record metric in matrix
            box_mat[gt_sts_ind,j] .= cos_sim(gt_params[gt_sts_ind,:],nano_params[nano_sts_ind,:])
            # Add to plot
            boxplot!(plt_aux,box_mat[gt_sts_ind,j],labels="",title=noise_labels[i],outliers=false)
        end
        # Push boxplot to array
        push!(plt_arr,plt_aux)
    end

    # Return array
    return plt_arr

end
function gen_ex_plt_array()

    println("Generating: E[X]")

    ## Read in ground-truth
    gt_ex_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_ex.txt"
    gt_ex_sts,gt_ex = read_in_exp_file(gt_ex_file)

    # Create plot array
    plt_arr = Vector{Any}()
    for (i,sig) in enumerate(noise_levels)
        println("Noise level: $(sig)")   
        ex_box_mat = fill(NaN,(length(gt_ex_sts),length(cov_labels)))
        smp_ex_box_mat = fill(NaN,(length(gt_ex_sts),length(cov_labels)))
        for (j,cove) in enumerate(cov_labels)
            # Get filenames
            ex_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_ex.txt"
            smp_ex_file = "$(aim_dir)/sample_averages/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_ex.txt"
            # Read in files
            ex_sts,ex = read_in_exp_file(ex_file)
            smp_ex_sts,smp_exs = read_in_exp_file(smp_ex_file)
            # Record performance metric for E[X] in matrix
            gt_sts_ind,ex_sts_ind = find_comm_regions(gt_ex_sts,ex_sts)
            ex_box_mat[gt_sts_ind,j] .= mae(gt_ex[gt_sts_ind],ex[ex_sts_ind]) ./ 2.0
            # Record performance metric for X̄ in matrix
            gt_sts_ind,smp_ex_sts_ind = find_comm_regions(gt_ex_sts,smp_ex_sts)
            smp_ex_box_mat[gt_sts_ind,j] .= mae(gt_ex[gt_sts_ind],smp_exs[smp_ex_sts_ind]) ./ 2.0
        end

        # Arrange E[X] data
        ex_data = vcat(ex_box_mat...)
        ex_labels = fill("E[X]",length(ex_data))
        ex_covs = vcat([fill(cove,size(ex_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(ex_data)
        ex_data = ex_data[keepind]
        ex_labels = ex_labels[keepind]
        ex_covs = ex_covs[keepind]

        # Arrange X̄ data
        smp_ex_data = vcat(smp_ex_box_mat...)
        smp_ex_labels = fill("Xbar",length(smp_ex_data))
        smp_ex_covs = vcat([fill(cove,size(smp_ex_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(smp_ex_data)
        smp_ex_data = smp_ex_data[keepind]
        smp_ex_labels = smp_ex_labels[keepind]
        smp_ex_covs = smp_ex_covs[keepind]

        # Push boxplot to array
        plt_data = vcat(ex_data,smp_ex_data)
        plt_labels = vcat(ex_labels,smp_ex_labels)
        plt_covs = vcat(ex_covs,smp_ex_covs)
        if i==1
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=(0.0,0.5),
                title=noise_labels[i],xlab="Coverage",ylab="Absolute error",outliers=false)
        else
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=(0.0,0.5),
                title=noise_labels[i],xlab="Coverage",ylab="Absolute error",label="",outliers=false)
        end
        push!(plt_arr,plt)
        
    end

    # Return array
    return plt_arr

end
function gen_exx_plt_array()

    println("Generating: E[XX]")

    ## Read in ground-truth
    gt_exx_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_exx.txt"
    gt_exx_sts,gt_exx = read_in_exp_file(gt_exx_file)

    # Create plot array
    plt_arr = Vector{Any}()
    for (i,sig) in enumerate(noise_levels)
        println("Noise level: $(sig)")   
        exx_box_mat = fill(NaN,(length(gt_exx_sts),length(cov_labels)))
        smp_exx_box_mat = fill(NaN,(length(gt_exx_sts),length(cov_labels)))
        for (j,cove) in enumerate(cov_labels)
            # Get filenames
            exx_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_exx.txt"
            smp_exx_file = "$(aim_dir)/sample_averages/sigma_$(sig)/gm12878_chr22_$(cov_labels[j])_exx.txt"
            # Read in files
            exx_sts,exx = read_in_exp_file(exx_file)
            smp_exx_sts,smp_exxs = read_in_exp_file(smp_exx_file)
            # Record metric for E[X] in matrix
            gt_sts_ind,exx_sts_ind = find_comm_regions(gt_exx_sts,exx_sts)
            exx_box_mat[gt_sts_ind,j] .= mae(gt_exx[gt_sts_ind],exx[exx_sts_ind]) ./ 2.0
            # Record metric for X̄ in matrix
            gt_sts_ind,smp_exx_sts_ind = find_comm_regions(gt_exx_sts,smp_exx_sts)
            smp_exx_box_mat[gt_sts_ind,j] .= mae(gt_exx[gt_sts_ind],smp_exxs[smp_exx_sts_ind]) ./ 2.0
        end

        # Arrange E[XX] data
        exx_data = vcat(exx_box_mat...)
        exx_labels = fill("E[XX]",length(exx_data))
        exx_covs = vcat([fill(cove,size(exx_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(exx_data)
        exx_data = exx_data[keepind]
        exx_labels = exx_labels[keepind]
        exx_covs = exx_covs[keepind]

        # Arrange X̄X̄ data
        smp_exx_data = vcat(smp_exx_box_mat...)
        smp_exx_labels = fill("XXbar",length(smp_exx_data))
        smp_exx_covs = vcat([fill(cove,size(smp_exx_box_mat)[1]) for cove in cov_labels]...)
        keepind = .!isnan.(smp_exx_data)
        smp_exx_data = smp_exx_data[keepind]
        smp_exx_labels = smp_exx_labels[keepind]
        smp_exx_covs = smp_exx_covs[keepind]

        # Push boxplot to array
        plt_data = vcat(exx_data,smp_exx_data)
        plt_labels = vcat(exx_labels,smp_exx_labels)
        plt_covs = vcat(exx_covs,smp_exx_covs)
        if i==1
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=(0.0,0.5),
                title=noise_labels[i],xlab="Coverage",ylab="Absolute error",outliers=false)
        else
            plt = groupedboxplot(plt_covs,plt_data,group=plt_labels,bar_width=0.5,ylim=(0.0,0.5),
                title=noise_labels[i],xlab="Coverage",ylab="Absolute error",label="",outliers=false)
        end
        push!(plt_arr,plt)
        
    end

    # Return array
    return plt_arr

end
function plt_params_scatter(cove,sig)

    # Find index of coverage label
    cov_ind = findfirst(isequal(cove),cov_levels)

    # Files
    gt_mod_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug.txt"
    nano_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[cov_ind])_theta.txt"
    
    # Params
    wgbs_sts,wgbs_params = read_in_model_file(gt_mod_file)
    nano_sts,nano_params = read_in_model_file(nano_file)

    # Find common regions
    wgbs_sts_ind,nano_sts_ind = find_comm_regions(wgbs_sts,nano_sts)
    wgbs_params = wgbs_params[wgbs_sts_ind,:]
    nano_params = nano_params[nano_sts_ind,:]

    # Plot scatter
    xlab = "Fine-grain"
    ylab = "Coarse-grain"
    p_a = plot(wgbs_params[:,1],nano_params[:,1],seriestype=:scatter,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(-10.0,10.0),ylim=(-10.0,10.0),label="",title="a");
    p_b = plot(wgbs_params[:,2],nano_params[:,2],seriestype=:scatter,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(-100.0,100.0),ylim=(-100.0,100.0),label="",title="b");
    p_c = plot(wgbs_params[:,3],nano_params[:,3],seriestype=:scatter,xlabel=xlab,ylabel=ylab,
        markersize=0.2,markeralpha=0.1,xlim=(-20.0,20.0),ylim=(-20.0,20.0),label="",title="c");
    plot(p_a,p_b,p_c,layout=(3,1),size=(500,800))
    savefig("$(aim_dir)/Scatter-Params-cov-$(cove)-s-$(sig)-Aim-2.pdf")
    savefig("$(aim_dir)/Scatter-Params-cov-$(cove)-s-$(sig)-Aim-2.png")

    # Return nothing
    return nothing

end
function plt_exp_scatter(cove,sig)

    # Ground-truth files
    gt_exx_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_exx.txt"
    gt_ex_file = "$(gt_dir)/GM12878_wgbs_cpelnano_theta_aug_ex.txt"
    gt_ex_sts,gt_ex = read_in_exp_file(gt_ex_file)
    gt_exx_sts,gt_exx = read_in_exp_file(gt_exx_file)
    
    # Find index of coverage label
    cov_ind = findfirst(isequal(cove),cov_levels)

    # Estimated
    ex_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[cov_ind])_ex.txt"
    exx_file = "$(aim_dir)/cpelnano/sigma_$(sig)/gm12878_chr22_$(cov_labels[cov_ind])_exx.txt"
    ex_sts,ex = read_in_exp_file(ex_file)
    exx_sts,exx = read_in_exp_file(exx_file)

    # Find common regions E[X]
    wgbs_sts_ex_ind,ex_sts_ind = find_comm_regions(gt_ex_sts,ex_sts)
    gt_ex = gt_ex[wgbs_sts_ex_ind]
    ex = ex[ex_sts_ind]
    
    # Find common regions E[XX]
    wgbs_sts_exx_ind,exx_sts_ind = find_comm_regions(gt_exx_sts,exx_sts)
    gt_exx = gt_exx[wgbs_sts_exx_ind]
    exx = exx[exx_sts_ind]

    # Transform to {0,1}
    ex = 0.5 .* (ex .+ 1.0)
    exx = 0.5 .* (exx .+ 1.0)
    gt_ex = 0.5 .* (gt_ex .+ 1.0)
    gt_exx = 0.5 .* (gt_exx .+ 1.0)

    # Plot scatter
    xlab = "Fine-grain"
    ylab = "Coarse-grain"
    p_ex = plot(gt_ex,ex,seriestype=:scatter,xlabel=xlab,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(0.0,1.0),ylim=(0.0,1.0),label="",title="E[X]");
    p_exx = plot(gt_exx,exx,seriestype=:scatter,xlabel=xlab,ylabel=ylab,markersize=0.2,
        markeralpha=0.1,xlim=(0.0,1.0),ylim=(0.0,1.0),label="",title="E[XX]");
    plot(p_ex,p_exx,layout=(2,1),size=(500,600))
    savefig("$(aim_dir)/Scatter-Exp-cov-$(cove)-s-$(sig)-Aim-2.pdf")
    savefig("$(aim_dir)/Scatter-Exp-cov-$(cove)-s-$(sig)-Aim-2.png")

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
plot(plt_arr...,layout=(3,2),size=(1000,1200))
savefig("$(aim_dir)/Boxplots-Params-Aim-2.pdf")
savefig("$(aim_dir)/Boxplots-Params-Aim-2.png")

## E[X] vs X̄

# Generate plot array
plt_arr = gen_ex_plt_array()

# Make plot
plot(plt_arr...,layout=(3,2),size=(1000,1200))
savefig("$(aim_dir)/Boxplots-EX-Aim-2.pdf")
savefig("$(aim_dir)/Boxplots-EX-Aim-2.png")

## E[XX] vs X̄

# Generate plot array
plt_arr = gen_exx_plt_array()

# Make plot
plot(plt_arr...,layout=(3,2),size=(1000,1200))
savefig("$(aim_dir)/Boxplots-EXX-Aim-2.pdf")
savefig("$(aim_dir)/Boxplots-EXX-Aim-2.png")

## Scatter plots for specific noise and cov

# Scatter parameters
println("Generating: Scatter parameters")
plt_params_scatter(10.0,1.0)
plt_params_scatter(10.0,2.0)

# Scatter Expectations
println("Generating: Scatter expectations")
plt_exp_scatter(10.0,1.0)
plt_exp_scatter(10.0,2.0)
