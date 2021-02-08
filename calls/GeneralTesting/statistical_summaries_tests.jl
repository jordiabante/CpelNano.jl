######################################################################################
## log(g(x)) & E[log(g(x))]
######################################################################################
using CpelNano

# Configuration
config = CpelNano.CpelNanoConfig()

# Parameter vector
a = 20.0; b = 100.0; c = 20.0; ϕ = [a,b,c];

# Other parameters
rs = CpelNano.RegStruct(); 
rs.cpg_pos = collect(100:20:3000); rs.N = length(rs.cpg_pos); 
rs.chrst = 1; rs.chrend = 3000; rs.ϕhat = ϕ;
rs.ρn = rand(0.001:0.001:0.1, rs.N); rs.dn = fill(20.0, rs.N - 1); 

# Re-scale model
CpelNano.rscle_grp_mod!(rs)

# Get matrices
CpelNano.get_rs_lgtrck_mats!(rs)

# Partition function
CpelNano.get_rs_logZ!(rs)

# Analysis region info
CpelNano.get_nls_reg_info!(rs,config)

# Get them all
CpelNano.get_rs_log_gs!(rs)

# Compute E[log g(X)]
CpelNano.get_rs_exp_log_g1!(rs)
CpelNano.get_rs_exp_log_g2!(rs)
rs.exps.log_g1
rs.exps.log_g2

# Compute E[X]
rs.exps.ex = CpelNano.get_E_X_log(rs.tm.log_u1, rs.tm.log_uN, rs.tm.log_Ws, rs.logZ)
rs.exps.exx = CpelNano.get_E_XX_log(rs.tm.log_u1, rs.tm.log_uN, rs.tm.log_Ws, rs.logZ)

# Compute MML
CpelNano.comp_mml!(rs)
rs.mml

# Compute NME
CpelNano.comp_nme!(rs)
rs.nme

# Compute CMD
CpelNano.comp_cmd(rs,rs)

######################################################################################
## Cross-entropy
######################################################################################
using Random
using CpelNano

# Configuration
config = CpelNano.CpelNanoConfig()

# Genomic properties
Random.seed!(123);
chrst = 1; chrend = 3000;
cpg_pos = collect(25:25:2975);
ρn = rand(0.001:0.001:0.1, length(cpg_pos));
dn = fill(25.0, length(cpg_pos) - 1)

# Parameter vectors
a1 = -10.0; b1 = -100.0; c1 = 20.0; ϕ1 = [a1,b1,c1];
a2 = 10.0; b2 = 100.0; c2 = 20.0; ϕ2 = [a2,b2,c2];

# Scale
ϕ1 = 0.5 .* ϕ1
ϕ2 = 0.5 .* ϕ2

## Region genomic properties

# Init structure 1
rs1 = CpelNano.RegStruct(); 
rs1.cpg_pos = cpg_pos; rs1.N = length(cpg_pos); rs1.chrst = chrst; 
rs1.chrend = chrend; rs1.ϕhat = ϕ1; rs1.ρn = ρn; rs1.dn = dn;

# Init structure 2
rs2 = CpelNano.RegStruct(); 
rs2.cpg_pos = cpg_pos; rs2.N = length(cpg_pos); rs2.chrst = chrst; 
rs2.chrend = chrend; rs2.ϕhat = ϕ2; rs2.ρn = ρn; rs2.dn = dn;

# Re-scale models
CpelNano.rscle_grp_mod!(rs1); CpelNano.rscle_grp_mod!(rs2);

# Get matrices
CpelNano.get_rs_lgtrck_mats!(rs1); CpelNano.get_rs_lgtrck_mats!(rs2);

# Partition function
CpelNano.get_rs_logZ!(rs1); CpelNano.get_rs_logZ!(rs2);

# Analysis region info
CpelNano.get_nls_reg_info!(rs1,config); CpelNano.get_nls_reg_info!(rs2, config);

# Get them all
CpelNano.get_rs_log_gs!(rs1); CpelNano.get_rs_log_gs!(rs2);

# Compute E[log g(X)]
CpelNano.get_rs_exp_log_g1!(rs1); CpelNano.get_rs_exp_log_g2!(rs1);
CpelNano.get_rs_exp_log_g1!(rs2); CpelNano.get_rs_exp_log_g2!(rs2);

# Compute E[X] and E[XX]
CpelNano.get_rs_exps!(rs1); CpelNano.get_rs_exps!(rs2);

# Mixture
mix = CpelNano.RegStruct()
mix.cpg_pos = rs1.cpg_pos; mix.N = rs1.N; 
mix.ρn = ρn; mix.dn = dn; mix.nls_rgs = rs1.nls_rgs;
mix.Nl = rs1.Nl; mix.L = rs1.L; mix.ρl = rs1.ρl; 
mix.dl = rs1.dl; mix.ϕhat = 0.5 .* (rs1.ϕhat + rs2.ϕhat);

# Update rest of information (CHECK WHAT'S NEEDED)
CpelNano.get_rs_lgtrck_mats!(mix)
CpelNano.get_rs_logZ!(mix)
CpelNano.get_rs_exps!(mix)
CpelNano.get_rs_log_gs!(mix)

# Mixture
exps_log_g1 = CpelNano.comp_crs_exps_log_g1(rs1, mix)
exps_log_g2 = CpelNano.comp_crs_exps_log_g2(rs2, mix)

# Cross-entropy w/ θ1 & θ2
cross_1 = CpelNano.comp_crs_ent(rs1, mix)
cross_2 = CpelNano.comp_crs_ent(rs2, mix)

######################################################################################
## CMD
######################################################################################
using Random
using CpelNano

# Configuration
config = CpelNano.CpelNanoConfig()

# Genomic properties
Random.seed!(123);
chrst = 1; chrend = 3000;
cpg_pos = collect(25:25:2975);
ρn = rand(0.001:0.001:0.1, length(cpg_pos));
dn = fill(25.0, length(cpg_pos) - 1)

# Parameter vectors
a1 = -10.0; b1 = -100.0; c1 = 20.0; ϕ1 = [a1,b1,c1];
a2 = 10.0; b2 = 100.0; c2 = 20.0; ϕ2 = [a2,b2,c2];

# Scale
ϕ1 = 0.5 .* ϕ1
ϕ2 = 0.5 .* ϕ2

## Region genomic properties

# Init structure 1
rs1 = CpelNano.RegStruct(); 
rs1.cpg_pos = cpg_pos; rs1.N = length(cpg_pos); rs1.chrst = chrst; 
rs1.chrend = chrend; rs1.ϕhat = ϕ1; rs1.ρn = ρn; rs1.dn = dn;

# Init structure 2
rs2 = CpelNano.RegStruct(); 
rs2.cpg_pos = cpg_pos; rs2.N = length(cpg_pos); rs2.chrst = chrst; 
rs2.chrend = chrend; rs2.ϕhat = ϕ2; rs2.ρn = ρn; rs2.dn = dn;

# Re-scale models
CpelNano.rscle_grp_mod!(rs1); CpelNano.rscle_grp_mod!(rs2);

# Get matrices
CpelNano.get_rs_lgtrck_mats!(rs1); CpelNano.get_rs_lgtrck_mats!(rs2);

# Partition function
CpelNano.get_rs_logZ!(rs1); CpelNano.get_rs_logZ!(rs2);

# Analysis region info
CpelNano.get_nls_reg_info!(rs1,config); CpelNano.get_nls_reg_info!(rs2, config);

# Get them all
CpelNano.get_rs_log_gs!(rs1); CpelNano.get_rs_log_gs!(rs2);

# Compute E[log g(X)]
CpelNano.get_rs_exp_log_g1!(rs1); CpelNano.get_rs_exp_log_g2!(rs1);
CpelNano.get_rs_exp_log_g1!(rs2); CpelNano.get_rs_exp_log_g2!(rs2);

# Compute E[X] and E[XX]
CpelNano.get_rs_exps!(rs1); CpelNano.get_rs_exps!(rs2);

# Compute MML
CpelNano.comp_mml!(rs1); CpelNano.comp_mml!(rs2);
round.([rs1.mml rs2.mml], digits=4)

# Compute NME
CpelNano.comp_nme!(rs1); CpelNano.comp_nme!(rs2);
round.([rs1.nme rs2.nme], digits=4)

# Compute CMD
cmd = CpelNano.comp_cmd(rs1, rs2)
round.(cmd, digits=4)
