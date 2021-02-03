########################################################################################################
# Generates replicates from two base model files for group comparison testing
########################################################################################################

# Dependencies
using Distributions
using DelimitedFiles

# Args
n_reps_g1 = parse(Int64, ARGS[1])
n_reps_g2 = parse(Int64, ARGS[2])

# Read base models
g1_mod = readdlm("examples/grp_comp_example/mod_files/g1_base_model.txt")
g2_mod = readdlm("examples/grp_comp_example/mod_files/g2_base_model.txt")

## Replicates

# Group 1
for rep = 1:n_reps_g1
    rep_mod = deepcopy(g1_mod)
    for i = 1:size(rep_mod)[1]
        θ = parse.(Float64, split(rep_mod[i,4], ','))
        rep_mod[i,4] = "$(θ[1] + rand(Normal())),$(θ[2] + rand(Normal())),$(θ[3] + rand(Normal()))"
    end
    writedlm("examples/grp_comp_example/mod_files/g1_rep$(rep).txt", rep_mod)
end

# Group 2
for rep = 1:n_reps_g2
    rep_mod = deepcopy(g2_mod)
    for i = 1:size(rep_mod)[1]
        θ = parse.(Float64, split(rep_mod[i,4], ','))
        rep_mod[i,4] = "$(θ[1] + rand(Normal())),$(θ[2] + rand(Normal())),$(θ[3] + rand(Normal()))"
    end
    writedlm("examples/grp_comp_example/mod_files/g2_rep$(rep).txt", rep_mod)
end
