#################################################################################################
# AIM 1: PR
#################################################################################################
## Deps
using Distributed
using CodecZlib
using StatsPlots
using Distributions
using DelimitedFiles
using Plots.PlotMeasures

## Constants
const noise_levels = [0.5,1.0,1.5,2.0,2.5,3.0]
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-1"

# Thresholds
const mega_thrs = collect(-3.0:0.01:0.0)
const deepsig_thrs = collect(0.0:0.01:1.0)
const nano_thrs = vcat(collect(-520.0:100.0:-20.0),collect(-20.0:0.25:20.0),collect(20.0:100.0:520.0))

# Colors
const blnd_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]
const cllr_color_code = Dict("Nanopolish"=>"#D55E00","DeepSignal"=>"#009E73","Megalodon"=>"#56B4E9")
const call_thrs = Dict("Nanopolish" => nano_thrs, "DeepSignal" => deepsig_thrs,"Megalodon" => mega_thrs)
const s_color_cde = Dict(0.5=>blnd_col[1],1.0=>blnd_col[2],1.5=>blnd_col[3],2.0=>blnd_col[4],2.5=>blnd_col[5],3.0=>blnd_col[6])

## Default attributes
default(titlefont=(14,"arial"),guidefont=(16,"arial"),tickfont=(12,"arial"))

# Function to threshold calls
@everywhere function thresh_calls(conf_x,thrs)

    # Init output vectors
    out_conf_x = []

    # Compare confidence score to threshold to determine call
    @inbounds for i=1:length(conf_x)

        if conf_x[i]=="n/a"
            # Fails to detect
            push!(out_conf_x,-1)
        elseif conf_x[i] < thrs
            # Negative
            push!(out_conf_x,-1)
        elseif conf_x[i] >= thrs
            # Positive
            push!(out_conf_x, 1)
        else
            println("Something went wrong when thresholding calls ...")
        end

    end

    return out_conf_x

end

# Function to compute sensitivity (recall, true positive rate)
@everywhere function comp_sens(truth,pred)

    # Get positives
    pos_ind = pred .== 1

    # Get true positives
    true_pos = sum(pos_ind)>0 ? sum(truth[pos_ind].==pred[pos_ind]) : 0.0

    # Return sensitivity
    return true_pos/sum(truth.==1)

end

# Function to compute precision
@everywhere function comp_prec(truth,pred)

    # Get positives
    pos_ind = pred .== 1

    # Get true positives
    true_pos = sum(pos_ind)>0 ? sum(truth[pos_ind].==pred[pos_ind]) : 0.0

    # Return precision
    return true_pos/sum(pos_ind)

end

@everywhere function pmap_thresh(true_x,conf_x,t)

    # Threshold calls
    thresh_x = thresh_calls(conf_x,t)

    # Return
    return comp_prec(true_x,thresh_x),comp_sens(true_x,thresh_x)

end

function get_prec_recall(caller,s)
    
    # Read in data
    stream = GzipDecompressorStream(open("$(data_dir)/$(caller)/gm12878_chr22_sigma_$(s)_$(caller)_tuples_wth_missed.tsv.gz"))
    in_data = readdlm(stream)
    close(stream)
        
    # Get true calls and confidence scores for predictions
    true_x = in_data[:,1]
    conf_x = in_data[:,2]

    # Get tpr and fpr
    pmap_out = pmap(t->pmap_thresh(true_x,conf_x,t),call_thrs[caller])

    # Disentangle data
    prec_vec = [x[1] for x in pmap_out]
    recall_vec = [x[2] for x in pmap_out]

    # Clean tuple
    kp_in = isfinite.(recall_vec) .& isfinite.(prec_vec)

    # Return tuple
    return recall_vec[kp_in],prec_vec[kp_in]

end

function get_auc(recall_vec,prec_vec)
    
    # Integrate ROC curve
    auc = 0.0
    @inbounds for i in 1:(length(prec_vec)-1)
        x = recall_vec[i]-recall_vec[i+1]
        y = prec_vec[i]
        auc += x * y
    end

    # Return AUC
    return auc

end

#################################################################################################
# ROC & AUC
#################################################################################################

# Init plots
p1 = plot(xlabel="sensitivity (recall)",ylabel="precision",title="Nanopolish",ylim=(0.5,1.0));
p2 = plot(xlabel="sensitivity (recall)",ylabel="precision",title="DeepSignal",ylim=(0.5,1.0));
p3 = plot(xlabel="sensitivity (recall)",ylabel="precision",title="Megalodon",ylim=(0.5,1.0));
roc_plots = Dict("Nanopolish" => p1,"DeepSignal" => p2, "Megalodon" => p3)

# Init AUROC array
auc_nano = []
auc_deep = []
auc_mega = []
aucs = Dict("Nanopolish" => auc_nano,"DeepSignal" => auc_deep, "Megalodon" => auc_mega)

# Calculate ROC and AUROC for each caller
for caller in ["Megalodon","DeepSignal","Nanopolish"]

    # Print caller
    println("Working on $(caller)")

    # Loop over noise levels
    for s in noise_levels

        # Print noise level
        println("Working on sigma=$(s)")

        # Get precision & recall
        recall_vec,prec_vec = get_prec_recall(caller,s)

        # Plot curve
        col = s_color_cde[s]
        plot!(roc_plots[caller],recall_vec,prec_vec,seriestype=:line,alpha=0.5,color=col,label="sd=$(s)")

        # Calculate AUC
        push!(aucs[caller],get_auc(recall_vec,prec_vec))

    end

end

# Plot ROC ROC_Curves
plot(p1,p2,p3,layout=(3,1),size=(600,1000),top_margin=10px,bottom_margin=10px,left_margin=20px,right_margin=20px,legend=:bottomleft)
savefig("$(data_dir)/PRC-Curves.pdf")

# Plot AUROC
p4 = plot(xlabel="signal noise level (sd)",ylabel="AUC-PRC",ylim=(0,1));
for caller in keys(cllr_color_code)
    col = cllr_color_code[caller]
    plot!(p4,noise_levels,aucs[caller],seriestype=:scatter,markershape=:circle,color=col,label=caller)
    plot!(p4,noise_levels,aucs[caller],seriestype=:line,alpha=0.5,color=col,label="")
end
plot(p4,size=(700,600),legend=:bottomleft)
savefig("$(data_dir)/AUC-PRC.pdf")
