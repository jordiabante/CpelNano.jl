######################################################################################
## log(g(x)) & E[log(g(x))]
######################################################################################
using CpelNano

# Configuration
config = CpelNano.CpelNanoConfig()

# Parameter vector
a = 20.0; b = 100.0; c = 20.0; ϕ = [a,b,c];

# Other parameters
rs=CpelNano.RegStruct(); 
rs.cpg_pos = collect(100:10:3000); rs.N=length(rs.cpg_pos); 
rs.chrst = 1; rs.chrend = 3000; rs.ϕhat=ϕ;
rs.ρn=rand(0.001:0.001:0.1,rs.N); rs.dn=fill(10.0,rs.N-1); 

# Re-scale model
CpelNano.rscle_grp_mod!(rs)

# Get matrices
CpelNano.get_rs_lgtrck_mats!(rs)

# Partition function
CpelNano.get_rs_logZ!(rs)

# Analysis region info
CpelNano.get_nls_reg_info!(rs,config)

# Get them all
CpelNano.get_log_gs!(rs)

# Compute E[log g(X)]
CpelNano.get_exp_log_g1!(rs)
CpelNano.get_exp_log_g2!(rs)
rs.exps.log_g1
rs.exps.log_g2

# Compute E[X]
rs.exps.ex = CpelNano.get_E_X_log(rs.tm.log_u1,rs.tm.log_uN,rs.tm.log_Ws,rs.logZ)
rs.exps.exx = CpelNano.get_E_XX_log(rs.tm.log_u1,rs.tm.log_uN,rs.tm.log_Ws,rs.logZ)

# Compute MML
CpelNano.comp_mml!(rs)
rs.mml

# Compute NME
CpelNano.comp_nme!(rs)
rs.nme
