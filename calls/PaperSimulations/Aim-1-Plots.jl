#################################################################################################
# Simulations were generated using 
# DeepSignal: r9.4_180mv_450bps_6mer (R9.4 pore model)
# DeepSimulator: R9.4 pore model
#################################################################################################
## Deps
using Distributed
@everywhere using StatsPlots
@everywhere using Distributions
@everywhere using DelimitedFiles

## Constants
#const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Caller-Error/"
const data_dir = "/Users/sandeepk/Downloads/Caller-Error"
const noise_levels = [0.5,1.0,1.5,2.0,2.5,3.0]
const blind_friend_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]
const color_code = Dict("nanopolish" => blind_friend_col[1],"deepsignal" => blind_friend_col[end], "megalodon" => blind_friend_col[2])
const nanopolish_thresholds = collect(-500.0:1:500.0)
const deepsignal_thresholds = collect(0.0:0.001:1.0)
const caller_thresholds = Dict("nanopolish" => nanopolish_thresholds, "deepsignal" => deepsignal_thresholds,"megalodon" => deepsignal_thresholds)
const color_code_noise = Dict(0.5 => blind_friend_col[1],1.0 => blind_friend_col[2], 1.5 => blind_friend_col[3],2.0 => blind_friend_col[4], 2.5 => blind_friend_col[5], 3.0 => blind_friend_col[6])

## Default attributes
default(titlefont=(14,"arial"),guidefont=(16,"arial"),tickfont=(12,"arial"))

# Function to compute accuracy
function trans_data(in_data)
    
    # Transform {0,1}→{-1,1}
    in_data[:,3] .= (2.0.*in_data[:,3].-1)

    # Return transformed data
    return in_data

end

# Function to threshold calls
function thresh_calls(true_x, conf_x, thresh)

    out_true_x = []
    out_conf_x = []

    # Compare confidence score to threshold to determine call
    for i in range(1, length=(length(conf_x)))
        if conf_x[i] >= thresh
            push!(out_conf_x, 1)
            push!(out_true_x, true_x[i])
        elseif conf_x[i] < thresh
            push!(out_conf_x, -1)
            push!(out_true_x, true_x[i])
        end
    end

    return out_true_x, out_conf_x

end

# Function to compute accuracy
function comp_accu(truth,pred) 

    # Return accuracy
    return sum(truth.==pred)/length(truth)

end

# Function to compute precision
function comp_prec(truth,pred)

    # Get positives
    pos_ind = pred .== 1

    # Get true positives
    true_pos = sum(truth[pos_ind].==pred[pos_ind])

    # Return precision
    return true_pos/sum(pos_ind)

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

function pmap_noise_ex(caller,s)
    
    # Read in 
    in_data = readdlm("$(data_dir)/$(caller)/noise_$(s)_$(caller)_tuples.txt")

    # Get data 
    true_x = in_data[:,1]
    pred_x = in_data[:,3]

    # Return tuple
    return comp_accu(true_x,pred_x),comp_prec(true_x,pred_x),comp_sens(true_x,pred_x),comp_spec(true_x,pred_x)

end

function pmap_noise_exx(caller,s)
    
    # Read in 
    in_data = readdlm("$(data_dir)/$(caller)/noise_$(s)_$(caller)_tuples.txt")

    # Get covariances 
    true_xx = in_data[1:(end-1),1] .* in_data[2:end,1]
    pred_xx = in_data[1:(end-1),3] .* in_data[2:end,3]

    # Return tuple
    return comp_accu(true_xx,pred_xx),comp_prec(true_xx,pred_xx),comp_sens(true_xx,pred_xx),comp_spec(true_xx,pred_xx)

end

#################################################################################################
# Performance of methylation callers with X_{n}
#################################################################################################

# Init plots
p1 = plot(ylabel="Accuracy (%)",title="Performance callers X_{n}",ylim=(0,100));
p2 = plot(ylabel="Precision (%)",ylim=(0,100));
p3 = plot(ylabel="Sensitivity (%)",ylim=(0,100));
p4 = plot(xlabel="Gaussian Noise at signal-level (\\sigma)",ylabel="Specificity (%)",ylim=(0,100));

for caller in keys(color_code)

    # Get error in call
    pmap_out = pmap(s->pmap_noise_ex(caller,s) ,noise_levels)
    
    # Unravel pmap out
    accu_vec = [x[1] for x in pmap_out]
    prec_vec = [x[2] for x in pmap_out]
    sens_vec = [x[3] for x in pmap_out]
    spec_vec = [x[4] for x in pmap_out]

    # Update plot
    col = color_code[caller]
    plot!(p1,noise_levels,accu_vec*100,seriestype=:scatter,markershape=:circle,color=col,label=caller)
    plot!(p1,noise_levels,accu_vec*100,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p2,noise_levels,prec_vec*100,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p2,noise_levels,prec_vec*100,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p3,noise_levels,sens_vec*100,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p3,noise_levels,sens_vec*100,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p4,noise_levels,spec_vec*100,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p4,noise_levels,spec_vec*100,seriestype=:line,alpha=0.5,color=col,label="")

end

# The sd that corresponds to the real base-calling accuracy (87%-90%) is σ=2.0-2.5
p_x = plot(p1,p2,p3,p4,layout=(4,1),size=(600,900))
savefig("/Users/jordiabante/Desktop/Benchmark-Callers-EX.pdf")

#################################################################################################
# Performance of methylation callers with X_{n}⋅X_{n+1}
#################################################################################################

# Init plots
p1 = plot(ylabel="Accuracy (%)",title="Performance callers X_{n} X_{n+1}",ylim=(0,100));
p2 = plot(ylabel="Precision (%)",ylim=(0,100));
p3 = plot(ylabel="Sensitivity (%)",ylim=(0,100));
p4 = plot(xlabel="Gaussian Noise at signal-level (\\sigma)",ylabel="Specificity (%)",ylim=(0,100));

# Loop over different callers
for caller in keys(color_code)

    # Get error in call
    pmap_out = pmap(s->pmap_noise_exx(caller,s) ,noise_levels)
    
    # Unravel pmap out
    accu_vec = [x[1] for x in pmap_out]
    prec_vec = [x[2] for x in pmap_out]
    sens_vec = [x[3] for x in pmap_out]
    spec_vec = [x[4] for x in pmap_out]
    
    # Update plot
    col = color_code[caller]
    plot!(p1,noise_levels,accu_vec*100,seriestype=:scatter,markershape=:circle,color=col,label=caller)
    plot!(p1,noise_levels,accu_vec*100,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p2,noise_levels,prec_vec*100,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p2,noise_levels,prec_vec*100,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p3,noise_levels,sens_vec*100,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p3,noise_levels,sens_vec*100,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p4,noise_levels,spec_vec*100,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p4,noise_levels,spec_vec*100,seriestype=:line,alpha=0.5,color=col,label="")

end

# The sd that corresponds to the real base-calling accuracy (87%-90%) is σ=2.0-2.5
p_xx = plot(p1,p2,p3,p4,layout=(4,1),size=(600,900))
savefig("/Users/jordiabante/Desktop/Benchmark-Callers-EXX.pdf")

#################################################################################################
# ROC & AUC with X_{n}
#################################################################################################

# Init plots
p1 = plot(ylabel="TPR",xlabel="FPR",title="Nanopolish ROC Curve");
p2 = plot(ylabel="TPR",xlabel="FPR",title="DeepSignal ROC Curve");
p3 = plot(ylabel="TPR",xlabel="FPR",title="Megalodon ROC Curve");
roc_plots = Dict("nanopolish" => p1,"deepsignal" => p2, "megalodon" => p3)

# Init AUROC array
auroc_nanopolish = []
auroc_deepsignal = []
auroc_megalodon = []
aurocs = Dict("nanopolish" => auroc_nanopolish,"deepsignal" => auroc_deepsignal, "megalodon" => auroc_megalodon)

# Calculate ROC and AUROC for each caller
for caller in keys(color_code)
    for s in noise_levels
        tpr_vec = []
        fpr_vec = []
        # Read in 
        in_data = readdlm("$(data_dir)/$(caller)/noise_$(s)_$(caller)_tuples.txt")
        
        # Transform {0,1}→{-1,1}
        in_data = trans_data(in_data)
        
        # Get true calls and confidence scores for predictions
        true_x = in_data[:,1]
        conf_x = in_data[:,2]

        # Get covariances 
        true_xx = in_data[1:(end-1),1] .* in_data[2:end,1]
        pred_xx = in_data[1:(end-1),3] .* in_data[2:end,3]

        # Threshold calls
        for t in caller_thresholds[caller]
            true_x_thresh, thresh_x = thresh_calls(true_x, conf_x, t, caller)
            push!(tpr_vec,comp_sens(true_x_thresh,thresh_x))
            push!(fpr_vec,1-comp_spec(true_x_thresh,thresh_x))

        end

        # Calculate AUROC
        auroc = 0
        for i in 1:(length(fpr_vec)-1)
            x = fpr_vec[i+1] - fpr_vec[i]
            y = tpr_vec[i]
            auroc = auroc + x * y
        end
        push!(aurocs[caller], auroc)

        col = color_code_noise[s]
        plot!(roc_plots[caller],fpr_vec,tpr_vec, seriestype=:line,alpha=0.5,color=col,label=string("noise ", s))

    end
end

# Plot ROC ROC_Curves
plot(p1, p2, p3, layout = (3,1)) #arg layout
savefig("/Users/sandeepk/Downloads/Caller-Error/ROC_Curves.pdf")

# Plot AUROC
p4 = plot(ylabel="AUROC)",title="Performance callers X_{n}",ylim=(0,1));
for caller in keys(color_code)
    col = color_code[caller]
    plot!(p4,noise_levels,-1*aurocs[caller],seriestype=:scatter,markershape=:circle,color=col,label=caller)
    plot!(p4,noise_levels,-1*aurocs[caller],seriestype=:line,alpha=0.5,color=col,label="")
end
plot(p4)
savefig("/Users/sandeepk/Downloads/Caller-Error/AUROC.pdf")

#################################################################################################
# ROC & AUC with X_{n}X_{n+1}
#################################################################################################
