#################################################################################################
# AIM 3:
# * Simulations were generated using 
# * DeepSignal: r9.4_180mv_450bps_6mer (R9.4 pore model)
# * DeepSimulator: R9.4 pore model
#################################################################################################
## Deps
using StatsPlots
using Distributions
using DelimitedFiles

## Constants
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Real-Data/GM12878/chr22/"
const aim_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-3/"
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
function find_comm_regions(sts1::Vector{Int64},sts2::Vector{Int64})::NTuple{2,Vector{Bool}}

    # Find common regions
    int_sts = intersect(sts1,sts2)
    sts1_ind = map(i-> sts1[i] in int_sts,1:length(sts1))
    sts2_ind = map(i-> sts2[i] in int_sts,1:length(sts2))

    # Return
    return sts1_ind,sts2_ind

end
function mse(mat1::Array{Float64,2},mat2::Array{Float64,2})::Vector{Float64}
    
    # Return mean-squared error
    return vec(sum((mat1-mat2).^2,dims=2)./3)

end
function rmse(mat1::Array{Float64,2},mat2::Array{Float64,2})::Vector{Float64}
    
    # Return root mean-squared error
    return sqrt(mse(mat1,mat2))

end
function cos_sim(mat1::Array{Float64,2},mat2::Array{Float64,2})::Vector{Float64}
    
    # Intermediate quantities
    num = sum(mat1 .* mat2,dims=2)
    den = sqrt.(sum(mat1.^2,dims=2)) .* sqrt.(sum(mat2.^2,dims=2))

    # Return cosine similarity
    return vec(num./den)

end
#################################################################################################
# Calls
#################################################################################################

## Read in

# Files
wgbs_path = "$(data_dir)/GM12878_wgbs_cpelnano_theta.txt"
nano_path = "$(data_dir)/GM12878_wgbs_cpelnano_theta.txt"
# nano_path = "$(data_dir)/chr22/GM12878_nanopore_cpelnano_theta.txt"

# Params
wgbs_sts,wgbs_params = read_in_model_file(wgbs_path)
nano_sts,nano_params = read_in_model_file(nano_path)

# Find common regions
wgbs_sts_ind,nano_sts_ind = find_comm_regions(wgbs_sts,nano_sts)
wgbs_params = wgbs_params[wgbs_sts_ind,:]
nano_params = nano_params[nano_sts_ind,:]

# Introduce noise (testing) <================================ REMOVE!!
nano_params .+= rand(Normal(0,0.5),size(nano_params)[1],size(nano_params)[2])

# Correlations
cor(wgbs_params,nano_params)

## Plots
sts = wgbs_sts[wgbs_sts_ind]

# Plot MSE
xlab = "Genomic coordinate (bp)"
ylab = "Mean-squared error (MSE)"
params_mse = mse(wgbs_params,nano_params)
plot(sts,params_mse,seriestype=:scatter,xlabel=xlab,ylabel=ylab,markersize=0.2,markeralpha=0.2,label="")

# Plot RMSE
xlab = "Genomic coordinate (bp)"
ylab = "Root mean-squared error (RMSE)"
params_rmse = sqrt.(params_mse)
plot(sts,params_rmse,seriestype=:scatter,xlabel=xlab,ylabel=ylab,markersize=0.2,markeralpha=0.2,label="")

# Plot cos similarity
xlab = "Genomic coordinate (bp)"
ylab = "Cosine similarity"
params_cos_sim = cos_sim(wgbs_params,nano_params)
plot(sts,params_cos_sim,seriestype=:scatter,xlabel=xlab,ylabel=ylab,markersize=0.2,markeralpha=0.2,
    ylim=(-1.0,1.0),label="")

# Plot scatter
xlab = "Fine-grain"
ylab = "Coarse-grain"
p_a = plot(wgbs_params[:,1],nano_params[:,1],seriestype=:scatter,ylabel=ylab,markersize=0.2,
    markeralpha=0.2,xlim=(-10.0,10.0),ylim=(-10.0,10.0),label="",title="a");
p_b = plot(wgbs_params[:,2],nano_params[:,2],seriestype=:scatter,ylabel=ylab,markersize=0.2,
    markeralpha=0.2,xlim=(-100.0,100.0),ylim=(-100.0,100.0),label="",title="b");
p_c = plot(wgbs_params[:,3],nano_params[:,3],seriestype=:scatter,xlabel=xlab,ylabel=ylab,markersize=0.2,
    markeralpha=0.2,xlim=(-20.0,20.0),ylim=(-20.0,20.0),label="",title="c");
plot(p_a,p_b,p_c,layout=(3,1),size=(500,800))
savefig("$(aim_dir)/Scatter-Aim-3.pdf")
savefig("$(aim_dir)/Scatter-Aim-3.png")
