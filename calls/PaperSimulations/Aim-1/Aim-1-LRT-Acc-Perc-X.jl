#################################################################################################
# AIM 1: nanopolish's LRT vs Accuracy
#################################################################################################
## Deps
using Dierckx
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
const s_color_cde = Dict("0.5" => blnd_col[1], "1.0" => blnd_col[2], "1.5" => blnd_col[3], "2.0" => blnd_col[4], 
                         "2.5" => blnd_col[5], "3.0" => blnd_col[6], "3.5" => blnd_col[7], "4.0" => blnd_col[8])

## Default attributes
default(titlefont=(14, "arial"),guidefont=(16, "arial"),tickfont=(12, "arial"))

#################################################################################################
# Functions
#################################################################################################

# Function to compute accuracy
function comp_accu(truth, pred) 

    # Return accuracy
    return sum(truth .== pred) / length(truth)

end

function get_perc_calls(s, lr_thrs)

    # Read in 
    stream = GzipDecompressorStream(open("$(data_dir)/nanopolish/gm12878_chr22_sigma_$(s)_nanopolish_tuples_sample.tsv.gz"))
    in_data = readdlm(stream)
    close(stream)
    
    # Remove NaN
    kp_ind = .! map(isequal("n/a"), in_data[:,2])
    in_data = in_data[kp_ind,:]

    # Get data
    perc_calls = []
    total_calls = size(in_data)[1]
    for lr in lr_thrs
        # Get LRs
        obs_lrs = in_data[:,2]
        lr_ind =  [abs(x) >= lr for x in obs_lrs]
        if sum(lr_ind) > 0
            # Push
            push!(perc_calls, sum(lr_ind) / total_calls)
        else
            # Push
            push!(perc_calls, 0)
        end
        println("$(lr)")
    end

    # Return
    return perc_calls

end

function get_acc_calls(s, lr_thrs)

    # Read in 
    stream = GzipDecompressorStream(open("$(data_dir)/nanopolish/gm12878_chr22_sigma_$(s)_nanopolish_tuples_sample.tsv.gz"))
    in_data = readdlm(stream)
    close(stream)
    
    # Remove NaN
    kp_ind = .! map(isequal("n/a"), in_data[:,2])
    in_data = in_data[kp_ind,:]

    # Get data
    accs = []
    for lr in lr_thrs
        # Get LRs
        obs_lrs = in_data[:,2]
        lr_ind =  [abs(x) >= lr for x in obs_lrs]
        if sum(lr_ind) > 0
            # Get calls
            true_x = in_data[lr_ind,1]
            pred_x = in_data[lr_ind,3]
            # Push
            push!(accs, comp_accu(true_x, pred_x))
        else
            # Push
            push!(accs, 1.0)
        end
        println("$(lr)")
    end

    # Return
    return accs

end

#################################################################################################
# Calls
#################################################################################################
max_lrt_plt = 20.0
lr_thrs = vcat(collect(0.0:0.1:25.0), collect(25.0:100.0:425.0))

# Plot percentage of calls
plt_perc = plot(xlabel="LRT", ylabel="percentage of calls", xlim=(0, max_lrt_plt), ylim=(0, 1.0));
for s in noise_levels
    println("Sigma $(s)")
    
    # Get data
    perc_calls = get_perc_calls(s, lr_thrs)
    
    # Plot spline
    plot!(plt_perc, lr_thrs, perc_calls, label="\\sigma=$(s)", color=s_color_cde["$(s)"])

end

# Add vertical line
vline!(plt_perc,[2.5],linestyle=:dash,color="black",label="LRT=2.5")

# Save
plot(plt_perc)
# savefig("$(data_dir)/Nanopolish-LRT-vs-Percentage-Calls-X.pdf")

# Plot accuracy of calls per noise level
plt_err = plot(xlabel="LRT", ylabel="percentage call error", xlim=(0, max_lrt_plt), ylim=(0, 0.31));
for s in noise_levels
    println("Sigma $(s)")
    
    # Get data
    acc_calls = get_acc_calls(s, lr_thrs)

    # Fit spline
    spl = Spline1D(lr_thrs, 1.0 .- acc_calls, k=3, s=1e-4)
    
    # Plot spline
    plot!(plt_err, lr_thrs, spl(lr_thrs), label="\\sigma=$(s)", color=s_color_cde["$(s)"])

end

# Add vertical line
vline!(plt_err,[2.5],linestyle=:dash,color="black",label="LRT=2.5")

# Save individual plot
plot(plt_err)
# savefig("$(data_dir)/Nanopolish-LRT-vs-Error-Calls-Per-Sigma-X.pdf")

# Save combined plot
plot(plt_err,plt_perc,layout=(2, 1),size=(700, 800))
# savefig("$(data_dir)/Nanopolish-LRT-vs-Error-LRT-vs-Perc-X.pdf")

# Plot accuracy of calls per noise level
plt_err_acc = plot(xlabel="percentage calls", ylabel="percentage call error", xlim=(0, 1), ylim=(0, 0.31), legend=:topleft);
for s in noise_levels
    println("Sigma $(s)")
    
    # Get data
    perc_calls = get_perc_calls(s, lr_thrs)
    acc_calls = get_acc_calls(s, lr_thrs)

    # Fit spline
    spl = Spline1D(1.0 .- perc_calls, 1.0 .- acc_calls, k=4, s=1e-2)
    
    # Plot spline
    # plot!(plt_err_acc, perc_calls, 1.0 .- acc_calls, label="\\sigma=$(s)", color=s_color_cde["$(s)"])
    plot!(plt_err_acc, 1.0 .- perc_calls, spl(perc_calls), label="\\sigma=$(s)", color=s_color_cde["$(s)"])

end

# Save
plot(plt_err_acc)
# savefig("$(data_dir)/Nanopolish-Perc-vs-Error-Calls-X.pdf")

plot(plt_err,plt_perc,plt_err_acc,layout=(3, 1),size=(700, 1000))
savefig("$(data_dir)/Nanopolish-LRT-All-X.pdf")
