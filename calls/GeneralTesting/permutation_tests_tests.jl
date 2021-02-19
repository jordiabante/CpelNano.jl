######################################################################################
## Tcmd matrix
######################################################################################
using CpelNano

# Parameter vectors
a1 = -10.0; b1 = -100.0; c1 = 20.0; ϕ1 = 0.5 * [a1,b1,c1];
a2 = 10.0; b2 = 100.0; c2 = 20.0; ϕ2 = 0.5 * [a2,b2,c2];

# Get group data
s1 = 4; s2 = 4; nse_sc = 1.0; seed = 123;
ms_g1, ms_g2 = CpelNano.gen_grp_comp_data(s1, s2, ϕ1, ϕ2, nse_sc, seed);

# Compute matrix of CMDs
cmd_tbl = CpelNano.get_unmat_cmd_tbl(ms_g1, ms_g2)

######################################################################################
## All unmatched tests
######################################################################################
using CpelNano

# Parameter vectors
a1 = -10.0; b1 = -100.0; c1 = 20.0; ϕ1 = 0.5 * [a1,b1,c1];
a2 = 10.0; b2 = 100.0; c2 = 20.0; ϕ2 = 0.5 * [a2,b2,c2];

# Get group data
s1 = 4; s2 = 4; nse_sc = 0.01; seed = 123;
ms_g1, ms_g2 = CpelNano.gen_grp_comp_data(s1, s2, ϕ1, ϕ2, nse_sc, seed);

# Compute statistic
test_cont_1 = CpelNano.unmat_est_reg_test(ms_g1, ms_g1)
test_cont_2 = CpelNano.unmat_est_reg_test(ms_g2, ms_g2)
test_treat = CpelNano.unmat_est_reg_test(ms_g1, ms_g2)

######################################################################################
## All matched tests
######################################################################################
using CpelNano

# Parameter vectors
a1 = -10.0; b1 = -100.0; c1 = 20.0; ϕ1 = 0.05 * [a1,b1,c1];
a2 = 10.0; b2 = 100.0; c2 = 20.0; ϕ2 = 0.05 * [a2,b2,c2];

# Get group data
s1 = 4; s2 = 4; nse_sc = 1.0; seed = 123;
ms_g1, ms_g2 = CpelNano.gen_grp_comp_data(s1, s2, ϕ1, ϕ2, nse_sc, seed);

# Perform testing
test_cont_1 = CpelNano.mat_est_reg_test(ms_g1, ms_g1)
test_cont_2 = CpelNano.mat_est_reg_test(ms_g2, ms_g2)
test_treat = CpelNano.mat_est_reg_test(ms_g1, ms_g2)

######################################################################################
## User interface
######################################################################################
using CpelNano

## Unmatched

# Fasta
fasta = "examples/full_example/reference/hg38_chr17_43023997_43145780.fa"

# Model files
mod_fls_g1 = ["examples/grp_comp_example/mod_files/g1_rep$(i).txt" for i = 1:6]
mod_fls_g2 = ["examples/grp_comp_example/mod_files/g2_rep$(i).txt" for i = 1:6]

# Config
config = CpelNano.CpelNanoConfig()
config.out_dir = "examples/grp_comp_example/tests/"
config.out_prefix = "cpelnano_grp_comp_unmat"
config.matched = false
config.verbose = true

# Differential analysis
CpelNano.diff_grp_comp(mod_fls_g1,mod_fls_g2,fasta,config)

## Matched

# Fasta
fasta = "examples/full_example/reference/hg38_chr17_43023997_43145780.fa"

# Model files
mod_fls_g1 = ["examples/grp_comp_example/mod_files/g1_rep$(i).txt" for i = 1:6]
mod_fls_g2 = ["examples/grp_comp_example/mod_files/g2_rep$(i).txt" for i = 1:6]

# Config
config = CpelNano.CpelNanoConfig()
config.out_dir = "examples/grp_comp_example/tests/"
config.out_prefix = "cpelnano_grp_comp_mat"
config.matched = true
config.verbose = true

# Differential analysis
CpelNano.diff_grp_comp(mod_fls_g1,mod_fls_g2,fasta,config)

#################################################################################################
## Two sample test
#################################################################################################
using Distributed
addprocs(8)
@everywhere using Pkg
@everywhere Pkg.activate("./")
@everywhere using CpelNano
using CpelNano

# Fasta
fasta = "examples/full_example/reference/hg38_chr17_43023997_43145780.fa"

# Model files
nano_s1 = "examples/full_example/nanopolish/full_example_noise2.0_methylation_calls.sorted.tsv"
nano_s2 = "examples/full_example/nanopolish/full_example_noise2.0_methylation_calls.sorted.tsv"
mod_path_s1 = "examples/two_samp_example/cpelnano/full_example_mod_nse_true_theta.txt"
mod_path_s2 = "examples/two_samp_example/cpelnano/full_example_mod_nse_true_theta.txt"

# Config
config = CpelNano.CpelNanoConfig()
# config.out_dir = "examples/two_samp_example/tests/"
config.out_dir = "/Users/jordiabante/Desktop/two_samp_example/"
config.out_prefix = "cpelnano_two_samp"
config.pval_comp = false
config.max_em_iters = 25
config.max_em_init = 4
config.verbose = true

# Differential analysis
CpelNano.diff_two_samp_comp(mod_path_s1,mod_path_s2,nano_s1,nano_s2,fasta,config)
