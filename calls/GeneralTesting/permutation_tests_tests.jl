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
## All two-sample tests
######################################################################################
using CpelNano

L = 20; N = L; M = 50; pobs = 1.0; μ_m = 2.0; σ_m = 2.0; μ_u = -2.0; σ_u = 2.0; 
Nl = fill(1.0, L); ρl = vcat(fill(1 / 1000, 8), fill(1 / 10, 4), fill(1 / 1000, 8)...); 
dl = fill(150.0, L - 1); chrst = 1; chrend = 3000; cpg_pos = collect(100:150:chrend);

# Config
config = CpelNano.CpelNanoConfig()
config.max_em_iters = 50
config.max_em_init = 8
config.verbose = true

# Parameter vector
a = 0.0; b = -50.0; c = 0.0; ϕ = [a,b,c]; 

# Other parameters
x = CpelNano.RegStruct(); x.N = N; x.L = L; x.Nl = Nl; x.ρl = ρl; x.dl = dl; 
α, β = CpelNano.get_αβ_from_ϕ(ϕ, x); 

# Simulate data
rs = CpelNano.cpel_samp_ont(M, α, β, pobs, μ_m, σ_m, μ_u, σ_u); 
rs.L = x.L; rs.Nl = x.Nl; rs.ρl = x.ρl; rs.dl = x.dl; rs.N = x.L; rs.ρn = x.ρl; rs.dn = x.dl; 
rs.chrst = chrst; rs.chrend = chrend; rs.cpg_pos = cpg_pos;

# Run EM algorithm
CpelNano.get_ϕhat!(rs,config); rs.ϕhat;
CpelNano.get_stat_sums!(rs, config);
rs_1 = rs;

# Parameter vector
a = 0.0; b = 50.0; c = 0.0; ϕ = [a,b,c]; 

# Other parameters
x = CpelNano.RegStruct(); x.N = N; x.L = L; x.Nl = Nl; x.ρl = ρl; x.dl = dl; 
α, β = CpelNano.get_αβ_from_ϕ(ϕ, x); 

# Simulate data
rs = CpelNano.cpel_samp_ont(M, α, β, pobs, μ_m, σ_m, μ_u, σ_u); 
rs.L = x.L; rs.Nl = x.Nl; rs.ρl = x.ρl; rs.dl = x.dl; rs.N = x.L; rs.ρn = x.ρl; rs.dn = x.dl; 
rs.chrst = chrst; rs.chrend = chrend; rs.cpg_pos = cpg_pos;

# Run EM algorithm
CpelNano.get_ϕhat!(rs,config); rs.ϕhat;
CpelNano.get_stat_sums!(rs, config);
rs_2 = rs;

# Config
config.LMAX_TWO_SAMP = 100
config.max_em_iters = 50
config.max_em_init = 2
config.verbose = false

# Perform testing
test_treat = CpelNano.pmap_diff_two_samp_comp_test(rs_1, rs_2, config)

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
