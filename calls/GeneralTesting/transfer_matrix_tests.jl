## Conditional Expectations E[X̄_{l}|Y;ϕ] and E[X̄_{l}X̄_{l+1}|Y;ϕ]
using CpelNano
using LinearAlgebra
using BenchmarkTools

# Parameters
L = 300; αs = fill(-1.0, L); βs = fill(0.0, L - 1);

# Data
obs = [CpelNano.MethCallCpgGrp(-20.0, 0.0) for _ in 1:L];

# Get u's
u1 = CpelNano.get_uc(αs[1], obs[1]); uN = CpelNano.get_uc(αs[end], obs[end]);
logu1 = CpelNano.get_log_uc(αs[1], obs[1]); loguL = CpelNano.get_log_uc(αs[end], obs[end]);

# Get W's
Ws = [CpelNano.get_Wc(αs[l], αs[l + 1], βs[l], obs[l], obs[l + 1]) for l = 1:(L - 1)]
logWs = [CpelNano.get_log_Wc(αs[l], αs[l + 1], βs[l], obs[l], obs[l + 1]) for l = 1:(L - 1)]

# Scaling
Ws = UniformScaling(1.0 / maximum(maximum.(Ws))) * Ws

# Get Zc
Zc = CpelNano.get_Zc(u1, uN, Ws)
logZc = CpelNano.get_log_Zc(logu1, loguL, logWs)

# Ex
CpelNano.get_Ec_X(u1,uN,Ws,Zc)
CpelNano.get_Ec_X_log(logu1,loguL,logWs,logZc)

# Exx
CpelNano.get_Ec_XX(u1,uN,Ws,Zc)
CpelNano.get_Ec_XX_log(logu1,loguL,logWs,logZc)

# E[logp(y|x)|y]
CpelNano.get_Ec_logpyx(u1,uN,Ws,Zc,obs)
CpelNano.get_Ec_logpyx_log(logu1,loguL,logWs,logZc,obs)

## Expectations E[X̄_{l};ϕ] and E[X̄_{l}X̄_{l+1};ϕ]
using CpelNano
using LinearAlgebra
using BenchmarkTools

# Parameters
L = 300; αs = fill(1.0, L); βs = fill(1.0, L - 1);

# Get u's
u1 = CpelNano.get_u(αs[1]); uL = CpelNano.get_u(αs[end]);
logu1 = CpelNano.get_log_u(αs[1]); loguL = CpelNano.get_log_u(αs[end]);

# Get W's
Ws = [CpelNano.get_W(αs[l], αs[l + 1], βs[l]) for l = 1:length(βs)];
logWs = [CpelNano.get_log_W(αs[l], αs[l + 1], βs[l]) for l = 1:length(βs)];

# Scaling
Ws = UniformScaling(1.0 / maximum(maximum.(Ws))) * Ws

# Get Z
Z = CpelNano.get_Z(u1, uL, Ws)
logZ = CpelNano.get_log_Z(logu1, loguL, logWs)

# E[X]
CpelNano.get_E_X(u1,uL,Ws,Z)
CpelNano.get_E_X_log(logu1,loguL,logWs,logZ)

# E[XX]
CpelNano.get_E_XX(u1,uL,Ws,Z)
CpelNano.get_E_XX_log(logu1,loguL,logWs,logZ)

# Marginal p(x)
CpelNano.get_marg_px(u1,uL,Ws,Z)
CpelNano.get_marg_px_log(logu1,loguL,logWs,logZ)
