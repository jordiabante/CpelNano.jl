#################################################################################################
# AIM 1: ROC
#################################################################################################
## Deps
# using Distributed
using StatsPlots
using Distributions
using DelimitedFiles
using Plots.PlotMeasures

## Constants
const noise_levels = [0.5,1.0,1.5,2.0,2.5,3.0]
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-1"

# Thresholds
const mega_thrs = collect(-3.0:0.03:0.0)
const deepsig_thrs = collect(0.0:0.01:1.0)
const nano_thrs = vcat(collect(-400.0:10.0:-20.0),collect(-20.0:0.1:20.0),collect(20.0:10.0:400.0))

# Colors
const blnd_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]
const call_thrs = Dict("nanopolish" => nano_thrs, "deepsignal" => deepsig_thrs,"megalodon" => mega_thrs)
const cllr_color_code = Dict("nanopolish" => blnd_col[1],"deepsignal" => blnd_col[end], "megalodon" => blnd_col[2])
const s_color_cde = Dict(0.5=>blnd_col[1],1.0=>blnd_col[2],1.5=>blnd_col[3],2.0=>blnd_col[4],2.5=>blnd_col[5],3.0=>blnd_col[6])

## Default attributes
default(titlefont=(14,"arial"),guidefont=(16,"arial"),tickfont=(12,"arial"))

# Function to threshold calls
function thresh_calls(true_x,conf_x,thrs)

    # Init output vectors
    out_true_x = []
    out_conf_x = []

    # Compare confidence score to threshold to determine call
    @inbounds for i=1:length(true_x)

        if conf_x[i]=="n/a"
            # Fails to detect
            push!(out_conf_x,-1)
            push!(out_true_x,true_x[i])
        elseif conf_x[i] < thrs
            # Negative
            push!(out_conf_x, -1)
            push!(out_true_x, true_x[i])
        elseif conf_x[i] >= thrs
            # Positive
            push!(out_conf_x, 1)
            push!(out_true_x, true_x[i])
        else
            println("Something went wrong when thresholding calls ...")
        end

    end

    return out_true_x,out_conf_x

end

# Function to compute sensitivity (recall, or true positive rate)
function comp_sens(truth,pred)

    # Get positives
    pos_ind = pred .== 1

    # Get true positives
    if length(truth[pos_ind]) != 0
        true_pos = sum(truth[pos_ind].==pred[pos_ind])
    else
        return 0.0
    end

    # Return sensitivity
    return true_pos/sum(truth.==1)

end

# Function to compute specificity (selectivity or true negative rate)
function comp_spec(truth,pred)

    # Get negatives
    neg_ind = pred .== -1

    # Get prediction & truth
    if length(truth[neg_ind]) != 0
        true_neg = sum(truth[neg_ind].==pred[neg_ind])
    else
        return 0.0
    end
    
    # Return specificity
    return true_neg/sum(truth.==-1)

end

function get_tpr_fpr(caller,s)
    
    # Read in data
    in_data = readdlm("$(data_dir)/$(caller)/gm12878_chr22_sigma_$(s)_$(caller)_tuples_wth_missed_sample.tsv")
        
    # Get true calls and confidence scores for predictions
    true_x = in_data[:,1]
    conf_x = in_data[:,2]
    
    # Threshold calls
    tpr_vec = []
    fpr_vec = []
    for t in call_thrs[caller]
        true_x_thresh,thresh_x = thresh_calls(true_x,conf_x,t)
        push!(tpr_vec,comp_sens(true_x_thresh,thresh_x))
        push!(fpr_vec,1.0-comp_spec(true_x_thresh,thresh_x))
    end

    # Return tuple
    return tpr_vec,fpr_vec

end

function get_auc(fpr_vec,tpr_vec)
    
    # Integrate ROC curve
    auc = 0.0
    @inbounds for i in 1:(length(fpr_vec)-1)
        x = fpr_vec[i]-fpr_vec[i+1]
        y = tpr_vec[i]
        auc += x * y
    end

    # Return AUC
    return auc

end

#################################################################################################
# ROC & AUC
#################################################################################################

# Init plots
p1 = plot(ylabel="TPR",xlabel="FPR",title="Nanopolish ROC Curve");
p2 = plot(ylabel="TPR",xlabel="FPR",title="DeepSignal ROC Curve");
p3 = plot(ylabel="TPR",xlabel="FPR",title="Megalodon ROC Curve");
roc_plots = Dict("nanopolish" => p1,"deepsignal" => p2, "megalodon" => p3)

# Init AUROC array
auc_nano = []
auc_deep = []
auc_mega = []
aucs = Dict("nanopolish" => auc_nano,"deepsignal" => auc_deep, "megalodon" => auc_mega)

# Calculate ROC and AUROC for each caller
for caller in keys(cllr_color_code)

    # Print caller
    println("Working on $(caller)")

    # Loop over noise levels
    for s in noise_levels

        # Print noise level
        println("Working on sigma=$(s)")

        # Get TPR & FPR
        tpr_vec,fpr_vec = get_tpr_fpr(caller,s)

        # Plot curve
        col = s_color_cde[s]
        plot!(roc_plots[caller],fpr_vec,tpr_vec,seriestype=:line,alpha=0.5,color=col,label=string("noise ",s))

        # Calculate AUC
        push!(aucs[caller],get_auc(fpr_vec,tpr_vec))

    end

end

# Plot ROC ROC_Curves
plot(p1,p2,p3,layout=(3,1),size=(600,1000),top_margin=10px,bottom_margin=10px,left_margin=20px,right_margin=20px)
savefig("$(data_dir)/ROC_Curves.pdf")

# Plot AUROC
p4 = plot(xlabel="Noise level",ylabel="AUC",title="Performance callers X_{n}",ylim=(0,1));
for caller in keys(cllr_color_code)
    col = cllr_color_code[caller]
    plot!(p4,noise_levels,aucs[caller],seriestype=:scatter,markershape=:circle,color=col,label=caller)
    plot!(p4,noise_levels,aucs[caller],seriestype=:line,alpha=0.5,color=col,label="")
end
plot(p4,size=(700,600))
savefig("$(data_dir)/AUROC.pdf")
