#################################################################################################
# AIM 1: X
#################################################################################################
## Deps
using CodecZlib
using StatsPlots
using Distributions
using DelimitedFiles
using Plots.PlotMeasures

## Constants
const noise_levels = ["$(s)" for s in collect(2.0:0.5:3.5)]
const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/Simulations/Aim-1"

# Colors
const blnd_col = ["#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7"]
const cllr_color_code = Dict("Nanopolish" => "#D55E00", "DeepSignal" => "#009E73", "Megalodon" => "#56B4E9")
const s_color_cde = Dict(0.5 => blnd_col[1], 1.0 => blnd_col[2], 1.5 => blnd_col[3], 2.0 => blnd_col[4], 
                         2.5 => blnd_col[5], 3.0 => blnd_col[6], 3.5 => blnd_col[7], 4.0 => blnd_col[8])

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

# Function to compute accuracy
function comp_accu(truth, pred) 

    # Return accuracy
    return sum(truth .== pred) / length(truth)

end

# Function to compute precision
function comp_prec(truth, pred)

    # Get positives
    pos_ind = pred .== 1

    # Get true positives
    true_pos = sum(truth[pos_ind] .== pred[pos_ind])

    # Return precision
    return true_pos / sum(pos_ind)

end

# Function to compute sensitivity (recall, or true positive rate)
function comp_sens(truth, pred)

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
function comp_spec(truth, pred)

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

function map_noise_ex(caller, s)
    
    # Print noise level
    println("Working on sigma=$(s)")

    # Read in 
    stream = GzipDecompressorStream(open("$(data_dir)/$(caller)/gm12878_chr22_sigma_$(s)_$(caller)_tuples_sample.tsv.gz"))
    in_data = readdlm(stream)
    close(stream)

    # Remove NaN
    kp_ind = .! map(isequal("n/a"), in_data[:,2])
    in_data = in_data[kp_ind,:]

    # Get data 
    true_x = in_data[:,1]
    pred_x = in_data[:,3]

    # Return tuple
    return comp_accu(true_x, pred_x), comp_prec(true_x, pred_x), comp_sens(true_x, pred_x), comp_spec(true_x, pred_x)

end

#################################################################################################
# Performance of methylation callers with X_{n}
#################################################################################################

# Init plots
p1 = plot(xlabel="signal noise level (sd)", ylabel="accuracy", ylim=(0.5, 1.0));
p2 = plot(xlabel="signal noise level (sd)", ylabel="precision", ylim=(0.5, 1.0));
p3 = plot(xlabel="signal noise level (sd)", ylabel="sensitivity", ylim=(0.5, 1.0));
p4 = plot(xlabel="signal noise level (sd)", ylabel="specificity", ylim=(0.5, 1.0));

for caller in ["Nanopolish"]

    # Print caller
    println("Working on $(caller)")
    
    # Get error in call
    map_out = map(s -> map_noise_ex(caller, s), noise_levels)
    
    # Unravel map out
    accu_vec = [x[1] for x in map_out]
    prec_vec = [x[2] for x in map_out]
    sens_vec = [x[3] for x in map_out]
    spec_vec = [x[4] for x in map_out]
    println("Acc.= $(accu_vec)\nPrec.= $(prec_vec)\nSens.=$(sens_vec)\nSpec.=$(spec_vec)")

    # Update plot
    col = cllr_color_code[caller]
    plot!(p1, noise_levels, accu_vec, seriestype=:scatter, markershape=:circle, color=col, label="")
    plot!(p1, noise_levels, accu_vec, seriestype=:line, alpha=0.5, color=col, label="")
    plot!(p2, noise_levels, prec_vec, seriestype=:scatter, markershape=:circle, color=col, label="")
    plot!(p2, noise_levels, prec_vec, seriestype=:line, alpha=0.5, color=col, label="")
    plot!(p3, noise_levels, sens_vec, seriestype=:scatter, markershape=:circle, color=col, label="")
    plot!(p3, noise_levels, sens_vec, seriestype=:line, alpha=0.5, color=col, label="")
    plot!(p4, noise_levels, spec_vec, seriestype=:scatter, markershape=:circle, color=col, label="")
    plot!(p4, noise_levels, spec_vec, seriestype=:line, alpha=0.5, color=col, label="")

end

# Generate plot & store
plot(p1,p2,p3,p4,layout=(4, 1),size=(600, 1000),top_margin=10px,bottom_margin=10px,left_margin=20px,right_margin=20px,legend=(0.2, 0.3))
savefig("$(data_dir)/Benchmark-Callers-X.pdf")
