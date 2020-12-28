## MML and NME

using CpelNano

# One α-subregion
𝒜s = [1:25]; ℬs = [1:24]; ϕ = [1.0,1.0];
x = CpelNano.RegStruct([],ϕ,𝒜s,ℬs);
CpelNano.get_∇logZ!(x); CpelNano.get_rs_mats!(x); CpelNano.get_Z!(x);
CpelNano.comp_mml!(x); CpelNano.comp_nme_vec!(x);
x.mml,x.nme_vec

## Output quantities
using CpelNano

# Normalized entropy
𝒜s=[1:1,2:2,3:18,19:19,20:20]; ℬs=[1:9,10:19];
ϕ=[4.0,0.0,0.0,-4.0,0.0,0.0,0.0]; x=CpelNano.RegStruct([],ϕ,𝒜s,ℬs);
CpelNano.get_∇logZ!(x); CpelNano.get_rs_mats!(x); CpelNano.get_Z!(x);
CpelNano.comp_mml!(x); CpelNano.get_log_gs!(x); CpelNano.comp_nme_vec!(x);
x.mml
x.nme_vec

## Differential quantities
using CpelNano

𝒜s=[1:1,2:2,3:18,19:19,20:20]; ℬs=[1:9,10:19];
ϕ1=[4.0,0.0,0.0,4.0,0.0,0.0,0.0];
ϕ2=[-4.0,0.0,0.0,-4.0,0.0,0.0,0.0];
x1=CpelNano.RegStruct([],ϕ1,𝒜s,ℬs); x2=CpelNano.RegStruct([],ϕ2,𝒜s,ℬs);
CpelNano.get_∇logZ!(x1); CpelNano.get_rs_mats!(x1); CpelNano.get_Z!(x1); 
CpelNano.get_log_gs!(x1); CpelNano.comp_mml!(x1); CpelNano.comp_nme_vec!(x1);
CpelNano.get_∇logZ!(x2); CpelNano.get_rs_mats!(x2); CpelNano.get_Z!(x2); 
CpelNano.get_log_gs!(x2); CpelNano.comp_mml!(x2); CpelNano.comp_nme_vec!(x2);
CpelNano.comp_gjsd_vec(x1,x1)
CpelNano.comp_gjsd_vec(x2,x2)
CpelNano.comp_gjsd_vec(x1,x2)
