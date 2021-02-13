## Estimation algorithm WGBS
using CpelNano

# Parameter vector
a = 0.25; b = -50.0; c = 5.0; ϕ = [a,b,c]; L = 200;

# Other parameters
x = CpelNano.RegStruct(); x.L = L; x.Nl = fill(1.0, x.L); 
x.ρl = rand(0.001:0.001:0.1, x.L); x.dl = fill(10.0, x.L - 1); 
α, β = CpelNano.get_αβ_from_ϕ(ϕ, x); 
config = CpelNano.CpelNanoConfig(); config.max_em_iters = 50; 
config.max_em_init = 10; config.verbose = true;

# Simulate data
M = 25; pobs = 0.25;
rs = CpelNano.cpel_samp_wgbs(M, α, β, pobs); 
rs.L = x.L; rs.Nl = x.Nl; rs.ρl = x.ρl; rs.dl = x.dl; 

# Run EM algorithm
CpelNano.get_ϕhat_wgbs!(rs,config); rs.ϕhat
α, β = CpelNano.get_αβ_from_ϕ(rs.ϕhat, rs)
CpelNano.gen_x_mc(α,β)

## EM instance
using CpelNano

# Parameter vector
a = 0.1; b = 0.2; c = 5.0; ϕ = [a,b,c];

# Other parameters
x = CpelNano.RegStruct(); x.L = 40; x.Nl = fill(5.0, x.L); x.ρl = rand(0.1:0.01:0.5, x.L); x.dl = fill(10.0, x.L - 1); 
α, β = CpelNano.get_αβ_from_ϕ(ϕ, x); M = 100; pobs = 1.0; μ_m = 40.0; σ_m = 2.0; μ_u = 80.0; σ_u = 2.0; 
config = CpelNano.CpelNanoConfig(); config.max_em_iters = 30; config.max_em_init = 1;

# Simulate data
rs = CpelNano.cpel_samp_ont(M, α, β, pobs, μ_m, σ_m, μ_u, σ_u); rs.L = x.L; rs.Nl = x.Nl; rs.ρl = x.ρl; rs.dl = x.dl; 

# Run instance
ϕhat, Q = CpelNano.em_inst(rs, zeros(3), config)

## Estimation algorithm
using CpelNano

# Parameter vector
a = 1.0; b = -20.0; c = 12.0; ϕ = [a,b,c]; L = 60;

# Other parameters
x = CpelNano.RegStruct(); x.L = L; x.Nl = fill(1.0, x.L); x.ρl = rand(0.001:0.001:0.1, x.L); x.dl = fill(10.0, x.L - 1); 
α, β = CpelNano.get_αβ_from_ϕ(ϕ, x); M = 40; pobs = 0.5; μ_m = 2.0; σ_m = 2.0; μ_u = -2.0; σ_u = 2.0; 
config = CpelNano.CpelNanoConfig(); config.max_em_iters = 50; config.max_em_init = 1;config.verbose = true;

# Simulate data
rs = CpelNano.cpel_samp_ont(M, α, β, pobs, μ_m, σ_m, μ_u, σ_u); 
rs.L = x.L; rs.Nl = x.Nl; rs.ρl = x.ρl; rs.dl = x.dl; 

# Run EM algorithm
CpelNano.get_ϕhat!(rs,config); rs.ϕhat
