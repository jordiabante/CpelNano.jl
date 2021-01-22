#################################################################################################
# AIM 1: XX
#################################################################################################
## Deps
using CodecZlib
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
const cllr_color_code = Dict("Nanopolish"=>"#D55E00","DeepSignal"=>"#009E73","Megalodon"=>"#56B4E9")
const call_thrs = Dict("Nanopolish" => nano_thrs, "DeepSignal" => deepsig_thrs,"Megalodon" => mega_thrs)
const s_color_cde = Dict(0.5=>blnd_col[1],1.0=>blnd_col[2],1.5=>blnd_col[3],2.0=>blnd_col[4],2.5=>blnd_col[5],3.0=>blnd_col[6])

## Default attributes
default(titlefont=(14,"arial"),guidefont=(16,"arial"),tickfont=(12,"arial"))

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

function map_noise_exx(caller,s)
    
    # Print noise level
    println("Working on sigma=$(s)")
    
    # Read in 
    stream = GzipDecompressorStream(open("$(data_dir)/$(caller)/gm12878_chr22_sigma_$(s)_$(caller)_tuples_wth_missed.tsv.gz"))
    in_data = readdlm(stream)
    close(stream)

    # Get covariances 
    true_xx = in_data[1:(end-1),1] .* in_data[2:end,1]
    pred_xx = in_data[1:(end-1),3] .* in_data[2:end,3]

    # Return tuple
    return comp_accu(true_xx,pred_xx),comp_prec(true_xx,pred_xx),comp_sens(true_xx,pred_xx),comp_spec(true_xx,pred_xx)

end

#################################################################################################
# Performance of methylation callers with X_{n}X_{n+1}
#################################################################################################

# Init plots
p1 = plot(xlabel="signal noise level (sd)",ylabel="accuracy",ylim=(0,1));
p2 = plot(xlabel="signal noise level (sd)",ylabel="precision",ylim=(0,1));
p3 = plot(xlabel="signal noise level (sd)",ylabel="sensitivity",ylim=(0,1));
p4 = plot(xlabel="signal noise level (sd)",ylabel="specificity",ylim=(0,1));

for caller in ["Megalodon","DeepSignal","Nanopolish"]

    # Print caller
    println("Working on $(caller)")

    # Get error in call
    map_out = map(s->map_noise_exx(caller,s),noise_levels)
    
    # Unravel map out
    accu_vec = [x[1] for x in map_out]
    prec_vec = [x[2] for x in map_out]
    sens_vec = [x[3] for x in map_out]
    spec_vec = [x[4] for x in map_out]

    # Update plot
    col = cllr_color_code[caller]
    plot!(p1,noise_levels,accu_vec,seriestype=:scatter,markershape=:circle,color=col,label=caller)
    plot!(p1,noise_levels,accu_vec,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p2,noise_levels,prec_vec,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p2,noise_levels,prec_vec,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p3,noise_levels,sens_vec,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p3,noise_levels,sens_vec,seriestype=:line,alpha=0.5,color=col,label="")
    plot!(p4,noise_levels,spec_vec,seriestype=:scatter,markershape=:circle,color=col,label="")
    plot!(p4,noise_levels,spec_vec,seriestype=:line,alpha=0.5,color=col,label="")

end

# Generate plot & store
plot(p1,p2,p3,p4,layout=(4,1),size=(600,1000),top_margin=10px,bottom_margin=10px,left_margin=20px,right_margin=20px,legend=(0.2,0.3))
savefig("$(data_dir)/Benchmark-Callers-XX-jordi.pdf")
