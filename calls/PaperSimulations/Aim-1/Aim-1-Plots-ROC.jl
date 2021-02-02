#################################################################################################
# AIM 1: ROC
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
const nano_thrs = vcat(collect(-500.0:10.0:-20.0), collect(-20.0:0.1:20.0), collect(20.0:10.0:500.0))

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

# Function to compute sensitivity (recall, or true positive rate)
@everywhere function comp_sens(truth, pred)

    # Get positives
    pos_ind = pred .== 1

    # Get true positives
    if length(truth[pos_ind]) != 0
        true_pos = sum(truth[pos_ind] .== pred[pos_ind])
    else
        return 0.0
    end

    # Return sensitivity
    return true_pos / sum(truth .== 1)

end

# Function to compute specificity (selectivity or true negative rate)
@everywhere function comp_spec(truth, pred)

    # Get negatives
    neg_ind = pred .== -1

    # Get prediction & truth
    if length(truth[neg_ind]) != 0
        true_neg = sum(truth[neg_ind] .== pred[neg_ind])
    else
        return 0.0
    end
    
    # Return specificity
    return true_neg / sum(truth .== -1)

end

@everywhere function pmap_thresh(true_x, conf_x, t)

    # Threshold calls
    thresh_x = thresh_calls(conf_x, t)

    # Return
    return comp_sens(true_x, thresh_x), 1.0 - comp_spec(true_x, thresh_x)

end

function get_tpr_fpr(caller, s)
    
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
    tpr_vec = [x[1] for x in pmap_out]
    fpr_vec = [x[2] for x in pmap_out]

    # Return tuple
    return tpr_vec, fpr_vec

end

function get_auc(fpr_vec, tpr_vec)
    
    # Integrate ROC curve
    auc = 0.0
    @inbounds for i in 1:(length(fpr_vec) - 1)
        x = fpr_vec[i] - fpr_vec[i + 1]
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
p = plot(ylabel="sensitivity", xlabel="1-specificity", title="ROC - Nanopolish", xlim=(0, 1));

# Init AUROC array
auc_nano = []

# Loop over noise levels
for s in noise_levels

    # Print noise level
    println("Working on sigma=$(s)")

    # Get TPR & FPR
    tpr_vec, fpr_vec = get_tpr_fpr("Nanopolish", s)

    # Plot curve
    col = s_color_cde[s]
    plot!(p, fpr_vec, tpr_vec, seriestype=:line, alpha=0.5, color=col, label="sd=$(s)")

    # Calculate AUC
    push!(auc_nano, get_auc(fpr_vec, tpr_vec))

end

# Plot ROC ROC_Curves
plot(p,size=(700, 600))
savefig("$(data_dir)/ROC-Curves.pdf")

# Plot AUROC
col = cllr_color_code["Nanopolish"]
p4 = plot(xlabel="signal noise level (sd)", ylabel="AUC-ROC", ylim=(0.5, 1.0));
plot!(p4, noise_levels, auc_nano, seriestype=:scatter, markershape=:circle, color=col, label="")
plot!(p4, noise_levels, auc_nano, seriestype=:line, alpha=0.5, color=col, label="")
plot(p4,size=(700, 600))
savefig("$(data_dir)/AUC-ROC.pdf")
