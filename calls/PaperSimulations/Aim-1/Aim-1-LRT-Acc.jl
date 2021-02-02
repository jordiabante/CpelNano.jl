#################################################################################################
# AIM 1: nanopolish's LRT vs Accuracy
#################################################################################################
## Deps
using CodecZlib
using StatsPlots
using Distributions
using DelimitedFiles
using Plots.PlotMeasures

## Constants
const noise_levels = ["$(s)" for s in collect(2.0:0.5:3.5)]
const lr_ints = [-1000:0.01:-10,-10:0.01:-7.5,-7.5:0.01:-5,-5:0.01:-2.5,
                 -2.5:0.01:-1.25,-1.25:0.01:-0.5,-0.5:0.01:0.5,0.5:0.01:1.25,1.25:0.01:2.5,
                 2.5:0.01:5,5:0.01:7.5,7.5:0.01:10,10:0.01:1000]
const lr_labs = ["$(minimum(x)):$(maximum(x))" for x in lr_ints]
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

function map_noise_ex(s, lr_ints)

    # Read in 
    stream = GzipDecompressorStream(open("$(data_dir)/nanopolish/gm12878_chr22_sigma_$(s)_nanopolish_tuples_sample.tsv.gz"))
    in_data = readdlm(stream)
    close(stream)
    
    # Remove NaN
    kp_ind = .! map(isequal("n/a"), in_data[:,2])
    in_data = in_data[kp_ind,:]

    # Get data
    accs = []
    n_pts = []
    total_pts = 0
    for lr in lr_ints
        # Get LRs
        obs_lrs = in_data[:,2]
        lr_ind =  [x in lr for x in obs_lrs]
        if sum(lr_ind) > 0
            # Get calls
            true_x = in_data[lr_ind,1]
            pred_x = in_data[lr_ind,3]
            # Update total calls
            total_pts += sum(lr_ind)
            # Push
            push!(n_pts, sum(lr_ind))
            push!(accs, comp_accu(true_x, pred_x))
        else
            # Push
            push!(n_pts, 0)
            push!(accs, 1.0)
        end
        println("$(lr)")
    end

    # Return accuracies
    return accs, n_pts / total_pts

end

#################################################################################################
# Functions
#################################################################################################

# Plot
n_bins_int = div(length(lr_ints), 2)
plt_crv = plot(xlabel="percentage of calls", ylabel="accuracy", xlim=(0, 1), ylim=(0.5, 1.0), legend=:bottomleft);
for s in noise_levels
    println("Sigma $(s)")
    
    # Get data
    accs, n_pts = map_noise_ex(s, lr_ints)
    
    # Plot bars
    p1 = plot(lr_labs,accs,seriestype=:bar,xrotation=45,ylabel="accuracy",title="\\sigma=$(s)",
        ylim=(0.45, 1),label="",top_margin=10px);
    p2 = plot(lr_labs,n_pts,seriestype=:bar,xrotation=45,xlabel="LR",ylabel="percentage of calls",
        ylim=(0.0, 0.51),label="");
    p = plot(p1, p2, layout=(2, 1), size=(700, 700))
    savefig("$(data_dir)/Nanopolish-LR-vs-Accuracy-Sigma-$(s).pdf")
    
    # Plot curve percent vs accuracy
    acc_curve = 0.5 * ( accs[1:(1 + n_bins_int)] + accs[end:-1:(1 + n_bins_int)])
    perc_curve = n_pts[1:n_bins_int] + n_pts[end:-1:(end - n_bins_int + 1)]
    push!(perc_curve, n_pts[n_bins_int + 1])
    perc_curve = [sum(perc_curve[1:i]) for i = 1:length(perc_curve)]
    acc_curve = [sum(acc_curve[1:i] .* perc_curve[1:i]) / sum(perc_curve[1:i]) for i = 1:length(perc_curve)]
    plot!(plt_crv, perc_curve, acc_curve, label="\\sigma=$(s)")

end
plot(plt_crv, size=(700, 700))
savefig("$(data_dir)/Nanopolish-Accuracy-vs-Percentage.pdf")
