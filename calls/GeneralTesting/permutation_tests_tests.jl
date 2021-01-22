######################################################################################
## Unmatched test Tcmd
######################################################################################
using CpelNano

# Parameter vectors
a1 = -10.0; b1 = -100.0; c1 = 20.0; ϕ1 = 0.05 * [a1,b1,c1];
a2 = 10.0; b2 = 100.0; c2 = 20.0; ϕ2 = 0.05 * [a2,b2,c2];

# Get group data
s1 = 4; s2 = 4; nse_sc = 10.0; seed = 123;
mods_g1, mods_g2 = CpelNano.gen_grp_comp_data(s1, s2, ϕ1, ϕ2, nse_sc, seed);

# Compute statistic
cmd_cont_1 = CpelNano.comp_unmat_stat_cmd(mods_g1, mods_g1)
cmd_cont_2 = CpelNano.comp_unmat_stat_cmd(mods_g1, mods_g1)
cmd_treat = CpelNano.comp_unmat_stat_cmd(mods_g1, mods_g2)

# Test
CpelNano.unmat_reg_test_tcmd(mods_g1,mods_g2)

######################################################################################
## All unmatched tests
######################################################################################
using CpelNano

# Parameter vectors
a1 = -10.0; b1 = -100.0; c1 = 20.0; ϕ1 = 0.05 * [a1,b1,c1];
a2 = 10.0; b2 = 100.0; c2 = 20.0; ϕ2 = 0.05 * [a2,b2,c2];
a1 = -10.0; b1 = -100.0; c1 = 20.0; ϕ1 = 0.05 * [a1,b1,c1];
a2 = 10.0; b2 = 100.0; c2 = 20.0; ϕ2 = 0.05 * [a2,b2,c2];

# Get group data
s1 = 4; s2 = 4; nse_sc = 10.0; seed = 123;
mods_g1, mods_g2 = CpelNano.gen_grp_comp_data(s1, s2, ϕ1, ϕ2, nse_sc, seed);

# Compute statistic
test_cont_1 = CpelNano.unmat_est_reg_test(mods_g1, mods_g1)
test_cont_2 = CpelNano.unmat_est_reg_test(mods_g2, mods_g2)
test_treat = CpelNano.unmat_est_reg_test(mods_g1, mods_g2)
