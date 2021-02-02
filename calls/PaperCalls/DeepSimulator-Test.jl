using Plots
using Statistics
using DelimitedFiles

const data_dir = "/Users/jordiabante/OneDrive - Johns Hopkins/CpelNano/Data/DeepSimulator/"

###############################################################################################
function read_in_seq(seq, σs)
    
    # Loop over noise levels
    mat_arr = []
    for σ in σs
        # Read in data
        norm_data = readdlm("$(data_dir)/deepsim_$(seq)_$(σ).txt")
        # Normalize data (from DeepSimulator code)
        norm_data =  (norm_data .- 14.0) / 5.7 
        # Push
        push!(mat_arr, norm_data)
    end
    
    # Return array of matrices
    return mat_arr

end
function read_in_seq_true_sig(index)
    
    # Read in data
    u_norm_data = readdlm("$(data_dir)/true_sd/deepsimu_current_fasta_$(index).txt")
    m_norm_data = readdlm("$(data_dir)/true_sd/deepsimu_current_fasta_$(index)m.txt")
    
    # Normalize data (from DeepSimulator code)
    u_norm_data =  (u_norm_data .- 14.0) / 5.7 
    m_norm_data =  (m_norm_data .- 14.0) / 5.7 

    # Return matrices
    return u_norm_data, m_norm_data

end
function print_grep_comm(seq)

    print("grep ")
    for i = 1:7
        print("-e $(seq[i:(i + 5)]) ")
    end

    return nothing

end
###############################################################################################

## AAAAACGAAAAA

# Read in
σs = [0.5,1.0,2.5,3.0]
u_seq = "AAAAACGAAAAA"
m_seq = "AAAAAMGAAAAA" # "AAAAAMGAAAAA"
print_grep_comm(u_seq)
print_grep_comm(m_seq)
u_mat_arr = read_in_seq(u_seq, σs)
m_mat_arr = read_in_seq(m_seq, σs)

# Plot standard deviation
kmer_ind = ["$i" for i = 1:12]
p = plot(xlabel="k-mer index", ylabel="sd", ylim=(0, 5));
sd_mat = fill(NaN, (12, length(σs)))
for (i, σ) in enumerate(σs)
    σ_dat = u_mat_arr[i]
    sd_mat[:,i] .= [std(σ_dat[j,:]) for j = 1:size(σ_dat)[1]]
    plot!(p, kmer_ind, sd_mat[:,i], label="\\sigma=$(σ)")
end
plot(p,size=(700, 600))
savefig("$(data_dir)/Observed_Stdev_$(u_seq).pdf")

# Plot
ps = []
μ_u = [83.4465,76.9282,87.1845,103.13,85.6849,80.6163,87.4578]
μ_m = [84.0604,80.6872,85.5882,101.933,84.1578,87.4578,79.6029]
kmer_ind = ["$i" for i = 1:12]
for (i, σ) in enumerate(σs)
    p = plot(xlab="k-mer index", ylab="current (pA)", title="\\sigma=$(σ)")
    plot!(p, kmer_ind[1:length(μ_u)], μ_u, color="blue", label="unmethylated")
    plot!(p, kmer_ind, u_mat_arr[i][:,1], color="blue", alpha=0.025, label="")
    plot!(p, kmer_ind, u_mat_arr[i], color="blue", alpha=0.025, label="")
    plot!(p, kmer_ind[1:length(μ_m)], μ_m, color="red", label="methylated")
    # plot!(p, kmer_ind, m_mat_arr[i][:,1], color="red", alpha=0.025,)
    # plot!(p, kmer_ind, m_mat_arr[i], color="red", alpha=0.025, label="")
    push!(ps, p)
end

# Plot array
plot(ps...,layout=(2, 2),size=(800, 700))
savefig("$(data_dir)/Current_$(u_seq).pdf")

## True sigma

u_mat_arr, m_mat_arr = read_in_seq_true_sig(1)
u_gt_mean = [83.4465,76.9282,87.1845,103.13,85.6849,80.6163]
m_gt_mean = [84.0604,80.6872,85.5882,101.933,84.1578,79.6029]
u_gt_sd = [1.65993, 1.70501, 1.55453, 2.24722, 2.36795, 1.67748]
m_gt_sd = [1.23283, 1.37188, 1.36722, 2.2369, 2.1373, 1.21716]

# Plot mean (scaled)
kmer_ind = ["$i" for i = 1:6]
mean_u = [mean(u_mat_arr[i,:]) for i = 1:6]
mean_m = [mean(m_mat_arr[i,:]) for i = 1:6]
p1 = plot(xlabel="k-mer index", ylabel="average current (pA)", ylim=(50, 110));
plot!(p1, kmer_ind, mean_u, color="blue",label="sample unmethylated");
plot!(p1, kmer_ind, u_gt_mean, seriestype=:scatter,label="true unmethylated",color="blue");
plot!(p1, kmer_ind, mean_m, color="red",label="sample methylated");
plot!(p1, kmer_ind, m_gt_mean, seriestype=:scatter,label="true methylated",color="red");

# Plot standard deviation (scaled)
kmer_ind = ["$i" for i = 1:6]
sd_u = [std(u_mat_arr[i,:]) for i = 1:6]
sd_m = [std(m_mat_arr[i,:]) for i = 1:6]
p2 = plot(xlabel="k-mer index", ylabel="sd", ylim=(0, 10));
plot!(p2, kmer_ind, sd_u, color="blue",label="");
plot!(p2, kmer_ind, 2.98 * u_gt_sd, seriestype=:scatter,label="",color="blue");
plot!(p2, kmer_ind, sd_m, color="red",label="");
plot!(p2, kmer_ind, 2.98 * m_gt_sd, seriestype=:scatter,label="",color="red");

plot(p1,p2,size=(700, 700),layout=(2, 1))
savefig("$(data_dir)/true_sd/mean_std_index_1.pdf")
