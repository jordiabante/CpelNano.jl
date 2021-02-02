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
const noise_levels = ["$(s)" for s in collect(2.0:0.5:3.5)]
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-1"

# Thresholds
const nano_thrs = vcat(collect(-500.0:50.0:-20.0), collect(-20.0:0.1:20.0), collect(20.0:50.0:500.0))

# Colors
const call_thrs = Dict("Nanopolish" => nano_thrs)
const cllr_color_code = Dict("Nanopolish" => "#D55E00")
const blnd_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]
const s_color_cde = Dict("0.5" => blnd_col[1], "1.0" => blnd_col[2], "1.5" => blnd_col[3], "2.0" => blnd_col[4], 
                         "2.5" => blnd_col[5], "3.0" => blnd_col[6], "3.5" => blnd_col[7], "4.0" => blnd_col[8])

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

# Function to threshold calls
@everywhere function thresh_calls(conf_x, thrs)

    # Init output vectors
    out_conf_x = []

    # Compare confidence score to threshold to determine call
    @inbounds for i = 1:length(conf_x)

        if conf_x[i] == "n/a"
            # Fails to detect
            push!(out_conf_x, 0)
        elseif conf_x[i] < thrs
            # Negative
            push!(out_conf_x, -1)
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
@everywhere function comp_sens(truth, pred)

    # Get positives
    pos_ind = pred .== 1

    # Get true positives
    true_pos = sum(pos_ind) > 0 ? sum(truth[pos_ind] .== pred[pos_ind]) : 0.0

    # Return sensitivity
    return true_pos / sum(truth .== 1)

end

# Function to compute precision
@everywhere function comp_prec(truth, pred)

    # Get positives
    pos_ind = pred .== 1

    # Get true positives
    true_pos = sum(pos_ind) > 0 ? sum(truth[pos_ind] .== pred[pos_ind]) : 0.0

    # Return precision
    return true_pos / sum(pos_ind)

end

@everywhere function pmap_thresh(true_x, conf_x, t)

    # Threshold calls
    thresh_x = thresh_calls(conf_x, t)

    # Return
    return comp_prec(true_x, thresh_x), comp_sens(true_x, thresh_x)

end

function get_prec_recall(caller, s)
    
    # Read in data
    stream = GzipDecompressorStream(open("$(data_dir)/$(caller)/gm12878_chr22_sigma_$(s)_$(caller)_tuples_sample.tsv.gz"))
    in_data = readdlm(stream)
    close(stream)
        
    # Remove NaN
    kp_ind = .! map(isequal("n/a"), in_data[:,2])
    in_data = in_data[kp_ind,:]

    # Get true calls and confidence scores for predictions
    true_x = in_data[:,1]
    conf_x = in_data[:,2]

    # Get tpr and fpr
    pmap_out = pmap(t -> pmap_thresh(true_x, conf_x, t), call_thrs[caller])

    # Disentangle data
    prec_vec = [round(x[1], digits=3) for x in pmap_out]
    recall_vec = [round(x[2], digits=3) for x in pmap_out]

    # Remove NaN
    kp_ind = isfinite.(recall_vec) .& isfinite.(prec_vec)
    recall_vec = recall_vec[kp_ind]
    prec_vec = prec_vec[kp_ind]

    # Clean data
    kp_ind = []
    sort_ind = sortperm(recall_vec)
    recall_vec = recall_vec[sort_ind]
    prec_vec = prec_vec[sort_ind]
    for i = 2:length(prec_vec)
        if prec_vec[i] < prec_vec[i - 1]
            push!(kp_ind, i)
        end
    end
    recall_vec = recall_vec[kp_ind]
    prec_vec = prec_vec[kp_ind]

    # Return tuple
    return recall_vec, prec_vec

end

function get_auc(recall_vec, prec_vec)
    
    # Integrate ROC curve
    auc = 0.0
    @inbounds for i in 1:(length(prec_vec) - 1)
        x = recall_vec[i + 1] - recall_vec[i]
        y = prec_vec[i]
        auc += x * y
    end

    # Return AUC
    return auc

end

#################################################################################################
# PRC & AUC
#################################################################################################

# Init plots
p = plot(xlabel="sensitivity (recall)", ylabel="precision", title="PRC - Nanopolish", xlim=(0.0, 1.0), ylim=(0.5, 1.0));

# Init AUROC array
auc_nano = []

# Loop over noise levels
for s in noise_levels

    # Print noise level
    println("Working on sigma=$(s)")

    # Get precision & recall
    recall_vec, prec_vec = get_prec_recall("Nanopolish", s)

    # Plot curve
    col = s_color_cde[s]
    plot!(p, recall_vec, prec_vec, seriestype=:line, alpha=0.5, color=col, label="sd=$(s)")

    # Calculate AUC
    push!(auc_nano, get_auc(recall_vec, prec_vec))

end

# Plot ROC ROC_Curves
plot(p, size=(700, 600))
savefig("$(data_dir)/PRC-Curves.pdf")

# Plot AUROC
col = cllr_color_code["Nanopolish"]
p4 = plot(xlabel="signal noise level (sd)", ylabel="AUC-PRC", ylim=(0.5, 1.0));
plot!(p4, noise_levels, auc_nano, seriestype=:scatter, markershape=:circle, color=col, label="")
plot!(p4, noise_levels, auc_nano, seriestype=:line, alpha=0.5, color=col, label="")
plot(p4, size=(700, 600), legend=:bottomleft)
savefig("$(data_dir)/AUC-PRC.pdf")
