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
    n_pairs = length(in_data[:,2]) - 1
    for lr in lr_thrs
        println("$(lr)")
        # Get valid pairs
        val_pairs = [all(abs.(in_data[i:(i + 1),2]) .>= lr) for i = 1:n_pairs]
        # Push
        push!(perc_calls, sum(val_pairs) / n_pairs)
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

    # Get all XX data
    n_pairs = length(in_data[:,2]) - 1
    true_xx = [in_data[i,1] * in_data[i + 1,1] for i = 1:n_pairs]
    pred_xx = [in_data[i,3] * in_data[i + 1,3] for i = 1:n_pairs]

    # Get data
    accs = []
    for lr in lr_thrs
        println("$(lr)")
        # Get valid pairs
        val_pairs = [all(abs.(in_data[i:(i + 1),2]) .>= lr) for i = 1:n_pairs]
        if sum(val_pairs) > 0
            # Keep valid pairs
            true_xx_thrs = true_xx[val_pairs]
            pred_xx_thrs = pred_xx[val_pairs]
            # Push
            push!(accs, comp_accu(true_xx_thrs, pred_xx_thrs))
        else
            # Push
            push!(accs, 1.0)
        end
    end

    # Return
    return accs

end

#################################################################################################
# Calls
#################################################################################################
max_lrt_plt = 20.0
lr_thrs = vcat(collect(0.0:0.1:10.0), collect(10.0:0.5:25.0), collect(25.0:25.0:125.0))

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
# savefig("$(data_dir)/Nanopolish-LRT-vs-Percentage-Calls-XX.pdf")

# Plot accuracy of calls per noise level
plt_err = plot(xlabel="LRT", ylabel="percentage call error", xlim=(0, max_lrt_plt), ylim=(0, 0.31));
for s in noise_levels
    println("Sigma $(s)")
    
    # Get data
    acc_calls = get_acc_calls(s, lr_thrs)

    # Fit spline
    spl = Spline1D(lr_thrs, 1.0 .- acc_calls, k=3, s=1e-2)
    
    # Plot spline
    plot!(plt_err, lr_thrs, spl(lr_thrs), label="\\sigma=$(s)", color=s_color_cde["$(s)"])

end

# Add vertical line
vline!(plt_err,[2.5],linestyle=:dash,color="black",label="LRT=2.5")

# Save individual plot
plot(plt_err)
# savefig("$(data_dir)/Nanopolish-LRT-vs-Error-Calls-Per-Sigma-XX.pdf")

# Plot accuracy of calls per noise level
plt_err_acc = plot(xlabel="percentage calls", ylabel="percentage call error", xlim=(0, 1), ylim=(0, 0.31), legend=:topleft);
for s in noise_levels
    println("Sigma $(s)")
    
    # Get data
    perc_calls = get_perc_calls(s, lr_thrs)
    acc_calls = get_acc_calls(s, lr_thrs)

    # Add (0,0) point
    push!(perc_calls, 0.0)
    push!(acc_calls, 1.0)

    # Fit spline
    spl = Spline1D(1.0 .- perc_calls, 1.0 .- acc_calls, k=4, s=1e-1)
    
    # Plot spline
    plot!(plt_err_acc, 1.0 .- perc_calls, spl(perc_calls), label="\\sigma=$(s)", color=s_color_cde["$(s)"])

end

# Save
plot(plt_err_acc)
# savefig("$(data_dir)/Nanopolish-Perc-vs-Error-Calls-XX.pdf")

plot(plt_err,plt_perc,plt_err_acc,layout=(3, 1),size=(700, 1000))
savefig("$(data_dir)/Nanopolish-LRT-All-XX.pdf")
